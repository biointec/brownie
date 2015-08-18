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

#include "kmernode.h"
#include "settings.h"
#include "iostream"

using namespace std;

// ============================================================================
// KMER NODE REFERENCE (PRIVATE)
// ============================================================================

const DSNode* KmerNodeRef::nodes = NULL;

// ============================================================================
// KMER NODE TABLE (PRIVATE)
// ============================================================================

KmerNodeRef KmerNodeTable::insert(const Kmer& kmer, NodeID id,
                                  PositionID pos, const DSNode &node)
{
        //ofstream kmerTablef;
        //kmerTablef.open("correctKmer.txt", std::ios_base::app);
	//kmerTablef<<"NodeID"<<","<<"kmer"<<","<<"position";

	// choose the right representative kmer
        Kmer kmerRC = kmer.getReverseComplement();
        bool reverse = (settings.isDoubleStranded()) && (kmerRC < kmer);
        const Kmer &reprKmer = (reverse) ? kmerRC : kmer;

        // translate the nodeID and position to that representative
        NodeID reprID = (reverse) ? -id : id;
        PositionID reprPos = (reverse) ? node.getLength() - pos -
                                         Kmer::getK() : pos;

        // insert value in table
        KmerNodeValue val(reprKmer, KmerNode(reprID, reprPos));
        pair<KmerNodeIt, bool> insResult = table->insert(val);
	string st=reprKmer.str();
	//kmerTablef<<reprID<<","<<st<<","<<reprPos<<endl;
	//kmerTablef.close();
        return KmerNodeRef(insResult.first, reverse);
}

// ============================================================================
// KMER NODE TABLE (PUBLIC)
// ============================================================================

KmerNodeTable::KmerNodeTable(const Settings& settings, NodeID numNodes) :
        settings(settings), numNodes(numNodes),
        table(NULL), remapInfo(NULL), timeStamp(0)
{
        // keep track of node remapping
        remapInfo = new vector<NodeEvent>[numNodes+1];
        timeStamp = 0;
}

KmerNodeTable::~KmerNodeTable()
{
        delete table;
        delete [] remapInfo;
}

void KmerNodeTable::recFindInTable(NodePosPair curr, size_t currTimeStamp,
                                   vector<NodePosPair>& npp) const
{
        if (curr.first == 0)
                return;

        const vector<NodeEvent> &cne = remapInfo[abs(curr.first)];
        for (size_t i = 0; i < cne.size(); i++) {
                const NodeEvent &e = cne[i];

                // only incorporate upstream events
                if (e.getTimeStamp() <= currTimeStamp)
                        continue;

                size_t newTarget = (curr.first > 0) ? e.getTargetID() :
                        -e.getTargetID();

                size_t newPosition = curr.second;
                newPosition += (curr.first > 0) ? e.getLeftOffset() :
                        e.getRightOffset();

                // incorporate the event
                recFindInTable(NodePosPair(newTarget, newPosition),
                               e.getTimeStamp(), npp);

                // if the event is not a duplication event, stop this branch
                if (!e.isDuplication())
                        return;
        }

        npp.push_back(curr);
}

NodePosPair KmerNodeTable::find(const Kmer& kmer) const
{
        // chose a representative kmer
        Kmer representative = settings.isDoubleStranded() ?
                kmer.getRepresentative() : kmer;

        bool reverse = (kmer != representative);

        // find the kmer in the table
        KmerNodeIt result = table->find(representative);

        // if it is not found, get out
        if (result == table->end())
                return NodePosPair(0, 0);

        // now find all occurences of the kmer in the table
        KmerNodeRef ref(result, reverse);
        return NodePosPair(ref.getNodeID(), ref.getPosition());
}

void KmerNodeTable::find(const Kmer& kmer, vector<NodePosPair>& npp) const
{
        npp.clear();

        // choose the right representative kmer
        Kmer kmerRC = kmer.getReverseComplement();
        bool reverse = (settings.isDoubleStranded()) && (kmerRC < kmer);
        const Kmer &reprKmer = (reverse) ? kmerRC : kmer;

        // find the kmer in the table
        KmerNodeIt result = table->find(reprKmer);

        // if it is not found, get out
        if (result == table->end())
                return;

        // now find all occurences of the kmer in the table
        KmerNodeRef ref(result, reverse);
        recFindInTable(NodePosPair(ref.getNodeID(), ref.getPosition()), 0, npp);
}

void KmerNodeTable::populateTable(const DSNode* nodes)
{
        //ofstream kmerTablef;
        //kmerTablef.open("correctKmer.txt");
	//kmerTablef<<"NodeID"<<","<<"kmer"<<","<<"position"<<endl;
        //kmerTablef.close();
	KmerNodeRef::setNodes(nodes);
        // count the number of k-mers in the graph
        size_t numKmers = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                const DSNode &node = nodes[id];
                if (!node.isValid())
                        continue;
                numKmers += node.getMarginalLength();
        }
        // allocate the new table
        if (table != NULL)
                delete table;
        table = new GoogleKmerNodeTable(numKmers);

        // populate the table with kmers
        for (NodeID id = 1; id <= numNodes; id++) {
                const DSNode &node = nodes[id];
                if (!node.isValid())
                        continue;
                const TString& tStr = node.getTSequence();
                Kmer kmer(tStr);
                PositionID pos = 0;
                insert(kmer, id, pos++, node);
                //comment by Madhi

		for (size_t i = Kmer::getK(); i < tStr.getLength(); i++) {
                        kmer.pushNucleotideRight(tStr[i]);
                        insert(kmer, id, pos++, node);
		        //comment by Madhi
                }
        }
}

void KmerNodeTable::mergeLeftToRight(NodeID leftID, NodeID rightID,
                                     size_t leftSize, size_t rightSize,
                                     bool duplication)
{
        timeStamp++;

        if (leftID > 0)
                remapInfo[leftID].push_back(NodeEvent(rightID, 0, rightSize,
                                                      duplication, timeStamp));
        else
                remapInfo[-leftID].push_back(NodeEvent(-rightID, rightSize, 0,
                                                       duplication, timeStamp));

        if (rightID > 0)
                remapInfo[rightID].push_back(NodeEvent(rightID, leftSize, 0,
                                                       false, timeStamp));
        else
                remapInfo[-rightID].push_back(NodeEvent(-rightID, 0, leftSize,
                                                        false, timeStamp));
}

void KmerNodeTable::sanityCheck(const DSNode* nodes, NodeID numNodes) const
{
        vector<NodePosPair> npp;

        for (KmerNodeIt it = table->begin(); it != table->end(); it++) {
                const Kmer& kmer = it->first;

                find(kmer, npp);

                for (size_t i = 0; i < npp.size(); i++) {
                        NodeID id = npp[i].first;
                        PositionID pos = npp[i].second;

                        if (nodes[abs(id)].getLength() > 5000)
                                continue;
                        string str = nodes[abs(id)].getSequence();

                        if (id < 0)
                                pos = nodes[abs(id)].getLength() - pos -
                                        Kmer::getK();

                        Kmer test(str.substr(pos, Kmer::getK()));

                        if (id < 0)
                                test.reverseComplement();

                        if (test != kmer)
                                cout << "ERROR ! " << endl;
                }
        }
}
