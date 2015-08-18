/***************************************************************************
 *   Copyright (C) 2014, 2015 Jan Fostier (jan.fostier@intec.ugent.be)     *
 *   Copyright (C) 2014, 2015 Mahdi Heydari (mahdi.heydari@intec.ugent.be) *
 *   This file is part of Brownie                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "graph.h"
#include "kmernode.h"
#include "readfile/fastafile.h"
#include "mappedread.h"
#include "settings.h"
#include <gsl/gsl_randist.h>
#include <cmath>
#include <fstream>

using namespace std;

// sort by coverage: highest coverage goes left !
bool sortBySize(const NodeID &lhs, const NodeID &rhs) {
        return (DBGraph::graph->getSSNode(lhs).getMarginalLength() >
                DBGraph::graph->getSSNode(rhs).getMarginalLength());
}

bool sortByStd(const pair<pair<NodeID, NodeID>, pair<double, double> >& left,
               const pair<pair<NodeID, NodeID>, pair<double, double> >& right) {
        return left.second.second < right.second.second;
}

// ============================================================================
// SPECIAL FUNCTIONS
// ============================================================================

std::ostream &operator<<(std::ostream &out, const vector<NodeID>& nodeSeq)
{
        for (int i = 0; i < (int)nodeSeq.size() - 1; i++)
                out << nodeSeq[i] << ", ";
        if (nodeSeq.size() > 0)
                out << nodeSeq[nodeSeq.size() - 1];
        out << " (" << nodeSeq.size() << ")";
        return out;
}

std::ostream &operator<<(std::ostream &out, const vector<NodeSpec>& nodeSeq)
{
        for (size_t i = 0; i < nodeSeq.size(); i++)
                out << nodeSeq[i].first << ", ";
        cout << "(" << nodeSeq.size() << ")";
        return out;
}

std::ostream &operator<<(std::ostream &out, const NodePair& np) {
        out << np.first << " -> " << np.second;
        return out;
}

// ============================================================================
// DE BRUIJN GRAPH CLASS
// ============================================================================

void DBGraph::readPairedReads(const ReadLibrary& input,
                              map<NodeID, vector<MappedRead> > &output)
{
        size_t readPairCount = 0;

        ifstream ifs(input.getMappedPairedFilename().c_str(), ios::binary);
        while(ifs.good()) {

                MappedRead mRead1, mRead2;
                TString read1, read2;

                mRead1.read(ifs);
                mRead2.read(ifs);

                read1.readFromStream(ifs);
                read2.readFromStream(ifs);

                if (readPairCount++ % OUTPUT_FREQUENCY == 0) {
                        cout << "\tProcessing pair " << readPairCount << "/"
                             << input.getNumMappedPaired() << "\r";
                        cout.flush();
                }

                if (!mRead1.isMapped())
                        continue;
                if (!mRead2.isMapped())
                        continue;

                const vector<NodeID>& seq1 = mRead1.getNodeSeq();
                const vector<NodeID>& seq2 = mRead2.getNodeSeq();

                for (size_t i = 0; i < seq1.size(); i++) {
                        NodeID nodeID = seq1[i];

                        output[nodeID].push_back(mRead1);
                        output[nodeID].push_back(mRead2);
                }

                for (size_t i = 0; i < seq2.size(); i++) {
                        NodeID nodeID = seq2[i];

                        output[-nodeID].push_back(mRead2.revCompl());
                        output[-nodeID].push_back(mRead1.revCompl());
                }
        }

        ifs.close();

        cout << "\tProcessing pair " << readPairCount << "/"
             << input.getNumMappedPaired() << endl;
}

Observations DBGraph::getUnbiasedDistance(NodeID srcID, NodeID dstID,
                                     const ReadLibrary& input,
                                     const vector<MappedRead>& mReads,
                                     double biasedPos) const
{
        Observations newEst;       // this is the output

        const double insertSize = input.getFragmentSize();
        const double insertStd = input.getFragmentStd();
        const double srcLength = getSSNode(srcID).getMarginalLength();

        const double minPos = biasedPos - insertSize +
                              3.0 * insertStd + Kmer::getWordSize();
        const double maxPos = biasedPos - insertSize +
                              getSSNode(dstID).getMarginalLength() +
                              Kmer::getWordSize() - 3.0 * insertStd;

        const double LB = srcLength - 3.0 * insertStd;
        const double UB = srcLength + insertSize - Kmer::getWordSize() + 3.0 * insertStd;

        // no unbiased estimation is possible
        if (maxPos < minPos)
                return newEst;

        for (size_t k = 0; k < mReads.size(); k = k + 2) {
                const MappedRead& mRead1 = mReads[k];
                const MappedRead& mRead2 = mReads[k+1];

                // for safety
                if (!mRead1.mapsToNode(srcID))
                        continue;

                // only select reads that contribute to the unbiased estimation
                double rootPos = mRead1.getNodeReadStartPos(srcID);
                if (rootPos < minPos || rootPos > maxPos)
                        continue;

                if (!mRead2.mapsToNode(dstID))
                        continue;

                NodeDistance nodeDist = mRead1.getDistBetweenNodes(srcID, dstID,
                                                                   mRead2, insertSize);

                double dist = nodeDist.getDistance();

                if (dist > LB && dist < UB) {   // discard outliers
                        newEst.addObsOnline(dist);
                }
        }

        return newEst;
}

void DBGraph::makeLinksConsistent(map<NodePair, GaussVal>& links)
{
        typedef map<NodePair, GaussVal>::iterator AnchorLinkIt;

        for (AnchorLinkIt it = links.begin(); it != links.end(); it++) {
                NodeID srcID = it->first.first;
                NodeID dstID = it->first.second;
                const GaussVal& gv = it->second;

                double srcLen = getSSNode(srcID).getMarginalLength();
                double dstLen = getSSNode(dstID).getMarginalLength();

                NodePair revNp = (dstID > 0) ? NodePair(dstID, srcID) :
                        NodePair(-dstID, -srcID);
                GaussVal revGv = (dstID > 0) ? GaussVal(-gv.getAverage(), gv.getStd()) :
                        GaussVal(gv.getAverage() - srcLen + dstLen, gv.getStd());

                AnchorLinkIt lIt = links.find(revNp);
                if (lIt == links.end()) {
                        //cout << "Inserted: " << revNp << ": " << revGv << endl;
                        links.insert(pair<NodePair, GaussVal>(revNp, revGv));
                        continue;
                }

                if (lIt->second != revGv) {
                        //cout << "Inconsitent distances: " << endl;
                        //cout << revNp << ": " << revGv << endl;
                        //cout << lIt->first << ": " << lIt->second << endl;
                }
        }
}

bool DBGraph::compareForOverlap(NodeID nodeA, const GaussVal& posA,
                                NodeID nodeB, const GaussVal& posB)
{
        double lenA = getSSNode(nodeA).getMarginalLength();
        double lenB = getSSNode(nodeB).getMarginalLength();

        double startA = posA.getAverage();
        double startB = posB.getAverage();

        double endA = startA + lenA;
        double endB = startB + lenB;

        double OL1 = endA - startB; // shift nodeA before nodeB
        double OL2 = endB - startA; // shift nodeB before nodeA

        // they don't overlap
        if (OL1 < 0 || OL2 < 0)
                return true;

        double displacement = min<double>(OL1, OL2);
        double std = sqrt(posA.getStd() * posA.getStd() + posB.getStd() * posB.getStd());

        double misplaced = displacement/ std;

        if (misplaced > SCAFFOLD_SNAP_LINK)
                return false;

        return true;
}

void DBGraph::dbReductions()
{
        set<Scaffold*> scaffoldSet;
        buildPrimaryScaffold(scaffoldSet);

        cout << "Number of scaffolds: " << scaffoldSet.size() << endl;
        for (set<Scaffold*>::iterator it = scaffoldSet.begin(); it != scaffoldSet.end(); it++) {
                validateScaffold(*(*it));
               // cout << *(*it) << endl;
        }

        for (set<Scaffold*>::iterator it = scaffoldSet.begin(); it != scaffoldSet.end(); it++)
                fillGaps(*(*it));

        cout << "Exiting... bye!" << endl;
        exit(0);
}
