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
#include "settings.h"
#include "library.h"

using namespace std;

void DBGraph::parseReads(size_t thisThread,
                         vector<string>& readBuffer)
{
        for (size_t i = 0; i < readBuffer.size(); i++) {
                const string& read = readBuffer[i];

                KmerIt it(read);
                if (!it.isValid())
                        continue;

                // increase the read start coverage (only for the first valid kmer)
                NodePosPair result = table->find(it.getKmer());
                if (result.getNodeID() != 0) {
                        SSNode node = getSSNode(result.getNodeID());
                        node.setReadStartCov(node.getReadStartCov()+1);
                }

                NodeID prevID = 0;
                for (KmerIt it(read); it.isValid(); it++ ) {
                        Kmer kmer = it.getKmer();
                        NodePosPair result = table->find(kmer);
                        if (!result.isValid()) {
                                prevID = 0;
                                continue;
                        }

                        // we've found the kmer, increase the node coverage
                        NodeID thisID = result.getNodeID();
                        SSNode node = getSSNode(thisID);
                        node.incKmerCov();

                        // if the previous node was valid and different, increase the arc coverage
                        if ((prevID != 0) && (prevID != thisID)) {
                                getSSNode(prevID).getRightArc(thisID)->incReadCov();
                                getSSNode(thisID).getLeftArc(prevID)->incReadCov();
                        }

                        if (it.hasRightOverlap())
                                prevID = thisID;
                        else
                                prevID = 0;
                }
        }
}

void DBGraph::workerThread(size_t thisThread, LibraryContainer* inputs)
{
        // local storage of reads
        vector<string> myReadBuf;

        size_t blockID, recordOffset;
        while (inputs->getReadChunk(myReadBuf, blockID, recordOffset))
                parseReads(thisThread, myReadBuf);
}

void DBGraph::countNodeandArcFrequency(LibraryContainer &inputs)
{
        const unsigned int& numThreads = settings.getNumThreads();

        cout << "Populating table... ";
        cout.flush();
        populateTable();
        cout << "done" << endl;

        cout << "Number of threads: " << numThreads << endl;
        cout << "Counting node and arc frequency: " << endl;

        inputs.startIOThreads(settings.getThreadWorkSize(),
                              settings.getThreadWorkSize() * settings.getNumThreads());

        // start worker threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&DBGraph::workerThread, this, i, &inputs);

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        inputs.joinIOThreads();

        depopulateTable();
}

// ============================================================================
// PRIVATE NODE COVERAGE / MULTIPLICITY
// ============================================================================

size_t DBGraph::getLowestArcMultiplicity ( NodeID left, NodeID right )
{
        SSNode node = getSSNode ( right );
        int lowArcMult = node.getLoExpMult();

        for ( ArcIt it = node.leftBegin(); it != node.leftEnd(); it++ ) {
                if ( it->getNodeID() == left ) {
                        continue;
                }

                lowArcMult -= getSSNode ( it->getNodeID() ).getHiExpMult();
        }

        return ( lowArcMult > 0 ) ? lowArcMult : 0;
}

// ============================================================================
// PRIVATE NODE COVERAGE ROUTINES
// ============================================================================

bool sortNodeByLength(const NodeID& left, const NodeID& right)
{
        return DBGraph::graph->getDSNode ( left ).getMarginalLength() >
               DBGraph::graph->getDSNode ( right ).getMarginalLength();
}

double DBGraph::getInitialEstimateForCoverage ( const ReadLibrary& input,
        vector<size_t> &readFreq ) const
{
        // sort the nodes from big to small + compute the total nodes length
        vector<NodeID> sortedNodes;
        sortedNodes.reserve ( numNodes );
        size_t totalSize = 0;
        for ( NodeID id = 1; id <= numNodes; id++ ) {
                sortedNodes.push_back ( id );
                totalSize += getDSNode ( id ).getMarginalLength();
        }

        DBGraph::graph = this;
        sort ( sortedNodes.begin(), sortedNodes.end(), sortNodeByLength );

        size_t RL = input.getAvgReadLength();
        size_t k = Kmer::getK();

        // for the biggest nodes compute the coverage
        size_t currSize = 0, numReads = 0;
        for ( size_t i = 0; i < sortedNodes.size(); i++ ) {
                size_t MNL = getDSNode ( sortedNodes[i] ).getMarginalLength();
                currSize += MNL + RL - k;
                numReads += readFreq[sortedNodes[i]];
                if ( currSize > 0.15*totalSize ) {
                        break;
                }
        }

        return double ( numReads ) / double ( currSize );
}

double DBGraph::estimateReadStartCoverage ( const ReadLibrary &input,
        const vector<size_t> &readFreq ) const
{
        size_t nom = 0, denom = 0;
        size_t k = Kmer::getK();
        size_t RL = input.getAvgReadLength();

        for ( NodeID id = 1; id <= numNodes; id++ ) {
                DSNode& node = getDSNode ( id );

                if ( !node.isValid() ) {
                        continue;
                }
                if ( node.getRoundMult() == 0 ) {
                        continue;
                }

                size_t MNL = node.getMarginalLength();

                nom += readFreq[id];
                denom += ( MNL + RL - k )  * node.getRoundMult();
        }

        return ( double ) nom / ( double ) denom;
}
