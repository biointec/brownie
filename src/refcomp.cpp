/***************************************************************************
 *   Copyright (C) 2014 - 2016 Jan Fostier (jan.fostier@intec.ugent.be)    *
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

#include "refcomp.h"
#include "graph.h"
#include "readfile/fastafile.h"

#include <iomanip>

using namespace std;

// ============================================================================
// REFERENCE COMPARISON CLASS
// ============================================================================

std::ostream &operator<<(std::ostream &out, const BreakPoint &bp)
{
        out << "RefSeq " << bp.refID
            << " [" <<  bp.begin << " - " << bp.end-1 << "]";
        return out;
}

std::ostream &operator<<(std::ostream &out, const RefSegment &rs)
{
        out << "\tRefSeq " << rs.getRefID() << " ";
        switch (rs.getType()) {
                case RefSegmentType::DELETION:
                        out << "DELETION  ";
                        break;
                case RefSegmentType::INSERTION:
                        out << "INSERTION ";
                        break;
                case RefSegmentType::PARALLEL:
                        out << "PARALLEL  ";
                        break;
                case RefSegmentType::BREAK:
                        out << "BREAK     ";
                        break;
                case RefSegmentType::CONTIG:
                        out << "CONTIG    ";
                        break;
        }

        out << "[" <<  rs.getRefBegin() << " - " << rs.getRefEnd() << "[ ";
        if (rs.getStartNodeID() == 0)
                out << "**" << " - ";
        else
                out << rs.getStartNodeID() << " (" << rs.getStartNodeBegin() << ") - ";

        if (rs.getEndNodeID() == 0)
                out << "**";
        else
                out << rs.getEndNodeID() << " (" << rs.getEndNodeEnd() << ")";

        return out;
}

void RefComp::printSegments(const std::vector<RefSegment>& segment_v) const
{
        for (const RefSegment& refSeg : segment_v)
                cout << refSeg << endl;
}

RefComp::RefComp(const std::string& refFilename)
{
        FastAFile ifs(false);
        ifs.open(refFilename.c_str());

        string read;
        while (ifs.getNextRead(read))
                refSeq_v.push_back(read);

        ifs.close();
}

size_t RefComp::getSize() const
{
        size_t ret = 0;

        for (size_t i = 0; i < refSeq_v.size(); i++)
                ret += refSeq_v[i].size();
        return ret;
}

void RefComp::alignReference(const DBGraph& dbg,
                             std::vector<RefSegment>& refSegment_v) const
{
        size_t currSize = 0, totSize = getSize();
        refSegment_v.clear();

        for (size_t refID = 0; refID < refSeq_v.size(); refID++) {
                const string& refSeq = refSeq_v[refID];

                KmerIt it(refSeq);
                NodePosPair prev = dbg.findNPP(it.getKmer());
                RefSegmentType type = (prev.isValid()) ? RefSegmentType::CONTIG : RefSegmentType::BREAK;

                refSegment_v.push_back(RefSegment(refID));
                refSegment_v.back().setBeginPos(it.getOffset(),
                                                prev.getNodeID(),
                                                prev.getPosition());
                refSegment_v.back().setType(type);

                // handle the other kmers in the reference sequence
                for (it++; it.isValid(); it++, currSize++) {
                        NodePosPair curr = dbg.findNPP(it.getKmer());

                        if (currSize % OUTPUT_FREQUENCY == 0)
                                cout << "Aligning reference to de Bruijn graph ("
                                     << setprecision(2)
                                     << Util::toPercentage(currSize, totSize) << "%)  \r";

                        // NOTE: it is possible that both if-statements below
                        // are executed consecutively. Do not change the order!

                        // encounter a breakpoint on a contig -> switch to break
                        if (type == RefSegmentType::CONTIG && !dbg.consecutiveNPP(prev, curr)) {
                                // close the contig
                                refSegment_v.back().setEndPos(it.getOffset(),
                                                              prev.getNodeID(),
                                                              prev.getPosition()+1);

                                // push a breakpoint
                                type = RefSegmentType::BREAK;
                                refSegment_v.push_back(RefSegment(refID));
                                refSegment_v.back().setBeginPos(it.getOffset(), 0, 0);
                                refSegment_v.back().setType(type);
                        }

                        // encounter a valid kmer with no open segment -> switch to contig
                        if (type == RefSegmentType::BREAK && curr.isValid()) {
                                // close the breakpoint
                                refSegment_v.back().setEndPos(it.getOffset(), 0, 0);

                                // push a contig
                                type = RefSegmentType::CONTIG;
                                refSegment_v.push_back(RefSegment(refID));
                                refSegment_v.back().setBeginPos(it.getOffset(),
                                                                curr.getNodeID(),
                                                                curr.getPosition());
                                refSegment_v.back().setType(type);
                        }

                        prev = curr;
                }

                // close open segment at the end of a sequence
                refSegment_v.back().setEndPos(refSeq.size(),
                                              prev.getNodeID(),
                                              prev.getPosition()+1);
        }

        cout << "Aligning reference to de Bruijn graph (100%)  " << endl;
}

void RefComp::annotateBreakpoints(const DBGraph& dbg,
                                  std::vector<RefSegment>& refSegment_v) const
{
        vector<size_t> dist_v(2*dbg.getNumNodes()+1, numeric_limits<size_t>::max());
        vector<bool> visited_v(2*dbg.getNumNodes()+1, false);

        for (size_t i = 0; i < refSegment_v.size(); i++) {
                RefSegment& rs = refSegment_v[i];

                cout << "Annotating breakpoints (" << setprecision(2)
                     << Util::toPercentage(i, refSegment_v.size()) << "%)  \r";

                if (rs.getType() != RefSegmentType::BREAK)
                        continue;

                size_t refID = rs.getRefID();
                bool haveLeft = ((i > 0) && (refSegment_v[i-1].getRefID() == refID));
                bool haveRight = ((i+1 < refSegment_v.size()) && (refSegment_v[i+1].getRefID() == refID));

                if (!haveLeft || !haveRight) {
                        rs.setType(RefSegmentType::DELETION);
                        continue;
                }

                NodeID srcID = refSegment_v[i-1].getEndNodeID();
                size_t srcPos = refSegment_v[i-1].getEndNodeEnd()-1;
                NodeID dstID = refSegment_v[i+1].getStartNodeID();
                size_t dstPos = refSegment_v[i+1].getStartNodeBegin();
                size_t refBegin = refSegment_v[i-1].getRefEnd();
                size_t refEnd = refSegment_v[i+1].getRefBegin();
                NodePosPair srcNPP(srcID, srcPos), dstNPP(dstID, dstPos);

                // try and find a path accross the breakpoint
                size_t maxLen = max<size_t>(100, 2*(refEnd - refBegin));
                bool pathExists = dbg.findPath(srcNPP, dstNPP, maxLen, dist_v, visited_v);
                if (!pathExists) {
                        dbg.revCompNPP(srcNPP);
                        dbg.revCompNPP(dstNPP);
                        pathExists = dbg.findPath(srcNPP, dstNPP, maxLen, dist_v, visited_v);
                }

                // no path found -> breakpoint
                if (!pathExists)
                        continue;

                if (refBegin == refEnd) {  // ref is fully present
                        rs.setType(RefSegmentType::INSERTION);
                } else {                   // ref is partially missing
                        if (srcNPP == dstNPP)
                                rs.setType(RefSegmentType::DELETION);
                        else
                                rs.setType(RefSegmentType::PARALLEL);
                }
        }

        cout << "Aligning reference to de Bruijn graph (100%)" << endl;
}

void RefComp::validateGraph(const DBGraph& dbg, size_t minContigSize)
{
        vector<RefSegment> refSegments_v;

        alignReference(dbg, refSegments_v);

        cout << "Annotating alignments" << endl;
        annotateBreakpoints(dbg, refSegments_v);
        printSegments(refSegments_v);
}

void RefComp::getTrueNodeChain(const DBGraph& dbg, vector<NodeChain>& nodeChain)
{
        // align the reference sequences to the DBG to figure out true node chains
        for (size_t refID = 0; refID < refSeq_v.size(); refID++) {
                const string& refSeq = refSeq_v[refID];

                // handle the other kmers
                NodePosPair prev;
                for (KmerIt it(refSeq); it.isValid(); it++) {
                        Kmer kmer = it.getKmer();

                        NodePosPair curr = dbg.findNPP(kmer);
                        if (!curr.isValid()) {

                                prev = curr;
                                continue;
                        }

                        if (dbg.consecutiveNPP(prev, curr)) {
                                if (curr.getNodeID() != prev.getNodeID())
                                        nodeChain.back().push_back(curr.getNodeID());
                        } else {
                                nodeChain.push_back(NodeChain());
                                nodeChain.back().setCount(1);
                                nodeChain.back().push_back(curr.getNodeID());
                        }

                        prev = curr;
                }
        }

}

void RefComp::getNodeMultiplicity(const DBGraph& dbg,
                                  vector<size_t>& multiplicity)
{
        multiplicity.clear();
        multiplicity.resize(dbg.getNumNodes()+1, 0);

        // count the number of true kmers in each node
        for (size_t refID = 0; refID < refSeq_v.size(); refID++) {

                const string& refSeq = refSeq_v[refID];

                // handle the first kmer separately
                KmerIt it(refSeq);
                NodePosPair prev = dbg.findNPP(it.getKmer());
                if (prev.isValid())
                        multiplicity[abs(prev.getNodeID())]++;

                // for all kmers in the reference sequence
                for (it++; it.isValid(); it++) {
                        Kmer kmer = it.getKmer();
                        NodePosPair curr = dbg.findNPP(kmer);

                        // update kmer counters
                        if (curr.isValid())
                                multiplicity[abs(curr.getNodeID())]++;

                        // set the true arc flag if appropriate
                        /*if (dbg.consecutiveNPP(prev, curr)) {
                                if ((prev.getPosition()+1) != curr.getPosition()) {
                                        SSNode left = dbg.getSSNode(prev.getNodeID());
                                        left.getRightArc(curr.getNodeID())->setTrueArc(true);
                                }
                        }*/

                        prev = curr;
                }
        }

        size_t numTrueNodes = 0, numErrNodes = 0;

        // compute the multiplicity for each node
        for (NodeID i = 1; i <= dbg.getNumNodes(); i++) {
                SSNode node = dbg.getSSNode(i);
                if (!node.isValid())
                        continue;

                double ML = node.getMarginalLength();
                multiplicity[i] = round((double)multiplicity[i]/ML);
                if (multiplicity[i] > 0)
                        numTrueNodes++;
                else
                        numErrNodes++;
        }


        cout << "Number of true nodes: " << numTrueNodes << "/" << dbg.getNumNodes() << " (" << Util::toPercentage(numTrueNodes, dbg.getNumNodes()) << "%)" << endl;
        cout << "Number of false nodes: " << numErrNodes << "/" << dbg.getNumNodes() << " (" << Util::toPercentage(numErrNodes, dbg.getNumNodes()) << "%)" << endl;
}


void RefComp::extractBreakpointSubgraph(const DBGraph& dbg, std::string breakpointFileName,string  tempDir){
        std::ifstream input(breakpointFileName);
        vector< pair<string, string> > breakpoints;
        if (!input.good()) {
                std::cerr << "Error opening: " << breakpointFileName << std::endl;
                return;
        }
        std::string line, id ="", DNA_sequence ="" ;
        while (std::getline(input, line).good()) {
                if (line[0] == '>') {
                        id = line.substr(1);
                        if (DNA_sequence!=""){
                                breakpoints.push_back( make_pair(id, DNA_sequence));
                        }
                        DNA_sequence.clear();
                }
                else if (line[0] != '>'){
                        DNA_sequence += line;
                }
        }
        if (DNA_sequence!=""){
                breakpoints.push_back( make_pair(id, DNA_sequence));
        }

        for (int i = 0 ;i <breakpoints.size(); i++){
                pair<string , string> breakPoint= breakpoints [i];
                set<int> nodeSet;
                extractNodeSetbySequence(dbg, breakPoint.second,nodeSet);
                dbg.writeCytoscapeGraph( tempDir+breakPoint.first,nodeSet,1);
        }
}

void RefComp::extractNodeSetbySequence(const DBGraph& dbg, std::string &refSeq, std::set<int> &nodeSet){

        for (KmerIt it(refSeq); it.isValid(); it++) {
                Kmer kmer = it.getKmer();
                NodePosPair curr = dbg.findNPP(kmer);
                if (!curr.isValid())

                        continue;
                nodeSet.insert(curr.getNodeID());
        }
}
