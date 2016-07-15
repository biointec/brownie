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

RefComp::RefComp(const std::string& refFilename)
{
        FastAFile ifs(false);
        ifs.open(refFilename.c_str());

        string read;
        while (ifs.getNextRead(read))
                reference.push_back(read);

        ifs.close();
}

void RefComp::validateGraph(const DBGraph& dbg)
{
        size_t numKmers = 0, numKmerFound = 0;
        vector<BreakPoint> breakpoint;

        // align the reference sequences to the DBG to figure out breakpoints
        for (size_t refID = 0; refID < reference.size(); refID++) {
                bool breakPointOpen = false;
                const string& refSeq = reference[refID];

                // handle the first kmer separately
                KmerIt it(refSeq);
                NodePosPair prev = dbg.findNPP(it.getKmer());
                if (prev.isValid()) {
                        numKmerFound++;
                } else {
                        breakpoint.push_back(BreakPoint(refID, 0));
                        breakPointOpen = true;
                }

                // handle the other kmers
                for (it++; it.isValid(); it++) {
                        Kmer kmer = it.getKmer();

                        NodePosPair curr = dbg.findNPP(kmer);

                        // update kmer counters
                        numKmers++;
                        if (curr.isValid())
                                numKmerFound++;

                        // if the kmer is connected to the previous one, get out
                        if (dbg.consecutiveNPP(prev, curr)) {
                                breakPointOpen = false;
                                prev = curr;
                                continue;
                        }

                        // update or create a new breakpoint
                        if (breakPointOpen) {
                                breakpoint.back().extendBreakPoint();
                        } else {
                                breakpoint.push_back(BreakPoint(refID, it.getOffset()));
                                breakPointOpen = true;
                        }

                        prev = curr;
                }
        }

        double fracFound = 100.0*(double)numKmerFound/(double) numKmers;

        cout.precision(5);
        cout << "\tValidation report: " << endl;
        cout << "\tNumber of k-mers in the DBG: " << numKmerFound << "/" << numKmers
             << " (" << fracFound << "%)" << endl;
        cout << "\tNumber of breakpoints: " << breakpoint.size() << endl;

        // figure out what kind of breakpoints we have
        for (const auto& it : breakpoint) {

                size_t begin = it.getBegin();
                size_t end = it.getEnd();
                const string& refSeq = reference[it.getRefID()];

                cout << "\t" << it << " ";

                if (begin == 0) {
                        cout << "**" << " - ";
                } else {
                        Kmer kmer(refSeq, begin-1);
                        NodePosPair npp = dbg.findNPP(kmer);
                        cout << npp.getNodeID() << " (" << npp.getPosition() << ") - ";
                }

                if (end >= (refSeq.size() - Kmer::getK() + 1)) {
                        cout << "**" << endl;
                } else {
                        Kmer kmer(refSeq, end);
                        NodePosPair npp = dbg.findNPP(kmer);
                        cout << npp.getNodeID() << " (" << npp.getPosition() << ")" << endl;
                }
        }
}

void RefComp::getTrueNodeChain(const DBGraph& dbg, vector<NodeChain>& nodeChain)
{
        // align the reference sequences to the DBG to figure out true node chains
        for (size_t refID = 0; refID < reference.size(); refID++) {
                const string& refSeq = reference[refID];

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
        for (size_t refID = 0; refID < reference.size(); refID++) {

                const string& refSeq = reference[refID];

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
                        if (dbg.consecutiveNPP(prev, curr)) {
                                if ((prev.getPosition()+1) != curr.getPosition()) {
                                        SSNode left = dbg.getSSNode(prev.getNodeID());
                                        left.getRightArc(curr.getNodeID())->setTrueArc(true);
                                }
                        }

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

        cout << "Number of true nodes in the graph: " << numTrueNodes
             << " number of false nodes: " << numErrNodes << endl;
}
