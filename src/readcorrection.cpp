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

#include <thread>
#include <string>
#include <iomanip>

#include "library.h"
#include "readcorrection.h"

using namespace std;

// ============================================================================
// ALIGNMENT METRICS CLASS
// ============================================================================

void AlignmentMetrics::addMetrics(const AlignmentMetrics& rhs)
{
        lock_guard<mutex> lock(metricMutex);
        numReads += rhs.numReads;
        numCorrReads += rhs.numCorrReads;
        numCorrByMEM += rhs.numCorrByMEM;
        numSubstitutions += rhs.numSubstitutions;
}

void AlignmentMetrics::printStatistics() const
{
        size_t numCorrByKmer = numCorrReads - numCorrByMEM;
        size_t numUncorrected = numReads - numCorrReads;

        cout << "\nError correction report:\n";
        cout << "\tNumber of reads handled: " << numReads << endl;
        cout << "\tNumber of corrected reads: " << numCorrReads
             << fixed << setprecision(2) << " ("
             << Util::toPercentage(numCorrReads, numReads) << "%)" << endl;
        cout << "\tNumber of reads corrected by kmer seeds: " << numCorrByKmer
             << fixed << setprecision(2) << " ("
             << Util::toPercentage(numCorrByKmer, numReads) << "%)" << endl;
        cout << "\tNumber of reads corrected by MEM seeds: " << numCorrByMEM
             << fixed << setprecision(2) << " ("
             << Util::toPercentage(numCorrByMEM, numReads) << "%)" << endl;
        cout << "\tNumber of substitutions in reads: " << numSubstitutions
             << fixed << setprecision(2) << " (avg of "
             << double(numSubstitutions)/double(numReads) << " per read)" << endl;
        cout << "\tNumber of uncorrected reads: " << numUncorrected
             << fixed << setprecision(2) << " ("
             << Util::toPercentage(numUncorrected, numReads) << "%)" << endl;
}

// ============================================================================
// SEED CLASS
// ============================================================================

void Seed::createNodePosition(const DBGraph& dbg, vector<NodePosPair>& npp) const
{
        size_t currReadPos = readFirst;
        for (vector<NodeID>::const_iterator it = nodeID.begin(); it != nodeID.end(); it++) {
                SSNode node = dbg.getSSNode(*it);
                size_t nodeOffset = (it == nodeID.begin()) ? nodeFirst : 0;
                size_t nodeOL = min(node.getMarginalLength() - nodeOffset, readEnd - currReadPos);
                for (size_t i = 0; i < nodeOL; i++)
                        npp[currReadPos+i] = NodePosPair(*it, i + nodeOffset);
                currReadPos += nodeOL;
        }
}

void Seed::mergeSeeds(const vector<Seed>& seeds,
                      vector<Seed>& mergedSeeds)
{
        if (seeds.empty())
                return;

        mergedSeeds.reserve(seeds.size());
        mergedSeeds.push_back(seeds.front());

        for (size_t i = 1; i < seeds.size(); i++) {
                const Seed& left = seeds[i-1];
                const Seed& right = seeds[i];

                bool consistent = false;
                if (left.nodeID.size() == 1 && right.nodeID.size() == 1)
                        if (left.nodeID.front() == right.nodeID.front())
                                if ((right.nodeFirst - left.nodeFirst) == (right.readFirst - left.readFirst))
                                        consistent = true;

                if (consistent)
                        mergedSeeds.back().readEnd = right.readEnd;
                else
                        mergedSeeds.push_back(right);
        }
}

// ============================================================================
// READ CORRECTION CLASS
// ============================================================================

void ReadCorrection::revCompl(vector< NodePosPair >& npp)
{
        reverse(npp.begin(), npp.end());

        for (size_t i = 0; i < npp.size(); i++) {
                if (!npp[i].isValid())
                        continue;

                NodeID nodeID = npp[i].getNodeID();
                size_t pos = npp[i].getPosition();
                const SSNode node = dbg.getSSNode(nodeID);
                npp[i] = NodePosPair(-nodeID, node.getMarginalLength() - 1 - pos);
        }
}

void ReadCorrection::findNPPSlow(const string& read, vector<NodePosPair>& npp)
{
        for (KmerIt it(read); it.isValid(); it++) {
                Kmer kmer = it.getKmer();
                NodePosPair result = dbg.findNPP(kmer);
                npp[it.getOffset()] = result;
        }
}

void ReadCorrection::findNPPFast(const string& read, vector<NodePosPair>& nppv)
{
        for (KmerIt it(read); it.isValid(); it++) {
                Kmer kmer = it.getKmer();
                NodePosPair npp = dbg.findNPP(kmer);
                nppv[it.getOffset()] = npp;

                if (!npp.isValid())
                        continue;

                NodeID nodeID = npp.getNodeID();
                const SSNode node = dbg.getSSNode(nodeID);

                size_t readPos = it.getOffset() + Kmer::getK();
                size_t nodePos = npp.getPosition() + Kmer::getK();

                while ((readPos < read.size()) && (nodePos < node.getLength())) {
                        if (read[readPos] == node.getNucleotide(nodePos)) {
                                it++;
                                nppv[it.getOffset()] = NodePosPair(nodeID, nodePos - Kmer::getK() + 1);
                                nodePos++;
                                readPos++;
                        } else
                                break;

                }
        }
}

void ReadCorrection::extractSeeds(const vector<NodePosPair>& nppv,
                                     vector<Seed>& seeds)
{
        size_t prev = nppv.size();
        for (size_t i = 0; i < nppv.size(); i++) {
                if (!nppv[i].isValid())
                        continue;

                // is it the first time we encounter a valid npp?
                if (prev == nppv.size()) {
                        prev = i;
                        seeds.push_back(Seed(nppv[i].getNodeID(), nppv[i].getPosition(), i, i + 1));
                        continue;
                }

                bool consistent = false;
                const SSNode thisNode = dbg.getSSNode(nppv[i].getNodeID());
                const SSNode prevNode = dbg.getSSNode(nppv[prev].getNodeID());

                // no, check for check for consistency
                if (thisNode.getNodeID() == prevNode.getNodeID()) {
                        if ((nppv[i].getPosition() - nppv[prev].getPosition()) == (i - prev))
                                consistent = true;
                } else {                // we're in different nodes
                        if (prevNode.getRightArc(thisNode.getNodeID()) != NULL) {
                                size_t thisPos = prevNode.getMarginalLength() + nppv[i].getPosition();
                                if ((thisPos - nppv[prev].getPosition()) == (i - prev)) {
                                        seeds.back().nodeID.push_back(thisNode.getNodeID());
                                        consistent = true;
                                }
                        }
                }

                prev = i;
                if (!consistent) {
                        seeds.push_back(Seed(nppv[i].getNodeID(), nppv[i].getPosition(), i, i + 1));
                } else {
                        seeds.back().readEnd = i + 1;
                }
        }
}

void ReadCorrection::recSearch(NodeID curr, string& read, vector<NodePosPair>& npp,
                                  size_t currReadPos, size_t& counter,
                                  int currScore, int& bestScore, size_t& seedLast, bool &fullyCorrected)
{
        const SSNode node = dbg.getSSNode(curr);

        counter++;
        if (counter > (size_t)settings.getReadCorrDFSNodeLimit())
        {
                fullyCorrected = false;
                return ;
        }
        vector<DFSNode> dfsNode;

        size_t readCharLeft = getMarginalLength(read) - currReadPos;

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                NodeID nextID = it->getNodeID();
                const SSNode nextNode = dbg.getSSNode(nextID);

                size_t OLSize = min(nextNode.getMarginalLength(), readCharLeft);
                size_t nextReadPos = currReadPos + OLSize;

                string nodeOL = nextNode.substr(Kmer::getK()-1, OLSize);
                string readOL = read.substr(currReadPos + Kmer::getK() - 1, OLSize);

                int thisScore = alignment.align(readOL, nodeOL);
                int nextScore = currScore + thisScore;
                float nextRelScore = (float)thisScore / (float)OLSize;

                dfsNode.push_back(DFSNode(nextID, nextReadPos, nextScore, nextRelScore));

                // =====================
                /*string str = nextNode.substr(Kmer::getK()-1, OLSize);
                string readSubStr = read.substr(currReadPos + Kmer::getK() -1, OLSize);
                for (size_t i = 0; i < currReadPos + Kmer::getK() - 1; i++)
                        cout << " ";
                cout << str << " (Node: " << nextID << ", curr: " << nextScore << ", best: " << bestScore << ", rel. score.:" << nextRelScore << ")" << endl;*/
                // ======================
        }

        sort(dfsNode.begin(), dfsNode.end());

        for (auto& it : dfsNode) {
                NodeID nextID = it.nodeID;
                int nextScore = it.score;
                size_t nextReadPos = it.readPos;
                size_t OLSize = nextReadPos - currReadPos;

                int maxAttainScore = nextScore + getMarginalLength(read) - nextReadPos;

                // save and, if necessary, update the best score
                int prevBestScore = bestScore;
                if ((nextScore > bestScore) /*&& (currReadPos + OLSize > seedLast)*/) {
                        bestScore = nextScore;
                        seedLast = currReadPos + OLSize;
                }

                // =====================
               /* SSNode nextNode = dbg.getSSNode(nextID);
                string str = nextNode.substr(Kmer::getK()-1, OLSize);
                string readSubStr = read.substr(currReadPos + Kmer::getK() -1, OLSize);
                for (size_t i = 0; i < currReadPos + Kmer::getK() - 1; i++)
                        cout << " ";
                cout << str << " (Node: " << nextID << ", curr: " << nextScore << ", best: " << bestScore << ", max: " << maxAttainScore << ")" << endl;*/
                // ======================

                // descend in a child node only if there is chance this will improve the best score
                if (maxAttainScore > bestScore)
                        recSearch(nextID, read, npp, nextReadPos, counter, nextScore, bestScore, seedLast,fullyCorrected);

                // if the best score has been updated in this branch save the npp
                if (bestScore > prevBestScore)
                        for (size_t i = 0; i < OLSize; i++)
                                npp[currReadPos+i] = NodePosPair(nextID, i);
        }
}

bool ReadCorrection::extendSeed(string& read, vector<NodePosPair>& npp,
                                   size_t& seedFirst, size_t& seedLast)
{
        // remove at most k nucleotides from the seed
        for (size_t i = 0; i < Kmer::getK(); i++) {
                if (seedLast == seedFirst + 1)
                        break;

                npp[seedLast - 1] = NodePosPair(0, 0);
                seedLast--;
        }

        // try to extend to the seed to the right within the current node
        const SSNode node = dbg.getSSNode(npp[seedLast-1].getNodeID());
        while ((seedLast < getMarginalLength(read)) &&
               (npp[seedLast - 1].getPosition() + 1 < node.getMarginalLength())) {
                npp[seedLast] = NodePosPair(node.getNodeID(), npp[seedLast-1].getPosition() + 1);
                seedLast++;
        }

        // try to find the best right path
        //if (seedLast < getMarginalLength(read))
        //      cout << read << endl;
        size_t counter = 0; int bestScore = -(getMarginalLength(read) - seedLast);
        bool fullyCorrected = true;
        if (seedLast < getMarginalLength(read))
                 recSearch(node.getNodeID(), read, npp, seedLast, counter, 0, bestScore, seedLast, fullyCorrected);
        return fullyCorrected;
}

void ReadCorrection::applyReadCorrection(string& read,
                                         const vector<NodePosPair>& npp,
                                         size_t seedFirst, size_t seedLast,
                                         vector<NodeID>& nodeChain)
{
        // correct the seed
        size_t curr = seedFirst;
        while (curr < seedLast) {
                NodeID currID = npp[curr].getNodeID();
                nodeChain.push_back(currID);

                SSNode node = dbg.getSSNode(currID);
                size_t strLen = min(node.getLength() - npp[curr].getPosition(), read.size() - curr);
                string str = node.substr(npp[curr].getPosition(), strLen);

                read.replace(curr, strLen, str);
                curr += (strLen - Kmer::getK() + 1);
        }
}

bool ReadCorrection::correctRead(string& read, vector<NodePosPair>& npp,
                                 size_t& first, size_t& last,
                                 vector<NodeID>& nodeChain)
{
        // extend to the right
        bool rightExt= extendSeed(read, npp, first, last);

        // reverse complement all data
        Nucleotide::revCompl(read);
        first = getMarginalLength(read) - first;
        last = getMarginalLength(read) - last;
        swap(first, last);
        revCompl(npp);

        // extend to the right (which used to be left)
        bool leftExt = extendSeed(read, npp, first, last);

        Nucleotide::revCompl(read);
        first = getMarginalLength(read) - first;
        last = getMarginalLength(read) - last;
        swap(first, last);
        revCompl(npp);

        // reverse complement all data
        string original = read;
        if (rightExt && leftExt)
               applyReadCorrection(read, npp, first, last, nodeChain);
        else
                return false;

        return true;
}

void ReadCorrection::findSeedKmer(const std::string& read,
                                     vector<Seed>& mergedSeeds)
{
        vector<NodePosPair> nppv(read.length() + 1 - Kmer::getK());

        // find the node position pairs using the kmer lookup table
        findNPPFast(read, nppv);

        // transform consistent npps to seeds
        vector<Seed> seeds;
        extractSeeds(nppv, seeds);

        // sort seeds according to nodeID
        sort(seeds.begin(), seeds.end());

        // merge seeds
        Seed::mergeSeeds(seeds, mergedSeeds);

        // ----------- OUTPUT ------------
        /*cout << "NPP after kmer lookup search" << endl;
        cout << read << endl;
        for (size_t i = 0; i < Kmer::getK() - 1; i++)
                cout << " ";
        for (size_t i = 0; i < nppv.size(); i++) {
                if (nppv[i].isValid())
                        //cout << i << ":" << npp[i].getNodeID() << ":" << npp[i].getOffset() << " ";
                        cout << "|";
                else
                        cout << "*";
        }
        cout << endl;
        for (size_t i = 0; i < mergedSeeds.size(); i++) {
                for (size_t j = 0; j < mergedSeeds[i].readFirst + Kmer::getK() - 1; j++)
                        cout << " ";
                cout << "[";
                for (size_t j = mergedSeeds[i].readFirst + 1; j < mergedSeeds[i].readEnd; j++)
                        cout << "-";
                cout << "[";
                cout << " (" << mergedSeeds[i].readFirst << " to " << mergedSeeds[i].readEnd << ")\t";
                for (size_t j = 0; j < mergedSeeds[i].nodeID.size(); j++)
                        cout << mergedSeeds[i].nodeID[j] << "\t";
                cout << endl;
        }*/
        // ----------- OUTPUT ------------
}

bool sortSeedByLength(const Seed& a, const Seed& b) {
        return ((a.readEnd-a.readFirst) > (b.readEnd-b.readFirst));
}

void ReadCorrection::findSeedMEM(const string& read,
                                    vector<Seed>& mergedSeeds)
{
        vector<match_t> matches;

        int memSize = Kmer::getK() - 1;
        while (matches.size() < 100 && memSize>5) {
                matches.clear();
                //cout << "Find MEM: " << memSize << endl;
                sa.findMEM(0l, read, matches, memSize, false);
                //cout << "Number of matches for size " << memSize << ": " << matches.size() << endl;
                memSize--;
        }


        vector<Seed> seeds;
        seeds.reserve(matches.size());
        for (auto it : matches) {

                vector<long>::const_iterator e = upper_bound(startpos.begin(), startpos.end(), it.ref);
                e--;
                NodeID nodeID = distance(startpos.begin(), e) + 1;
                size_t nodeFirst = it.ref - *e;

                SSNode node = dbg.getSSNode(nodeID);
                if (nodeFirst > node.getLength()) {
                        nodeID = -nodeID;
                        nodeFirst = nodeFirst - node.getLength() - 1;
                }

                // don't select MEMs that are closer than k nucleotides to the
                // right edge of a node as this MEM is also in its right neighbor
                if (nodeFirst >= node.getMarginalLength())
                        continue;

                size_t readFirst = it.query;

                // don't select MEMs that are closer than k nucleotides to the
                // right edge of the read: NECESSARY ?!?
                if (readFirst >= getMarginalLength(read))
                        continue;

                size_t readEnd = readFirst + it.len; //(it.len > Kmer::getK() ? it.len +1 - Kmer::getK() : 1);

                seeds.push_back(Seed(nodeID, nodeFirst, readFirst, readEnd));
        }

        // sort seeds according to nodeID
        sort(seeds.begin(), seeds.end());

        // merge seeds
        Seed::mergeSeeds(seeds, mergedSeeds);

        // sort seeds according to length
        sort(mergedSeeds.begin(), mergedSeeds.end(), sortSeedByLength);

        // retain only 10 largest seeds
        if (mergedSeeds.size() > 10)
                mergedSeeds.resize(10);

        for (size_t i = 0; i < mergedSeeds.size(); i++)
                mergedSeeds[i].readEnd -= min(Kmer::getK() - 1, mergedSeeds[i].readEnd - mergedSeeds[i].readFirst - 1);


        // ----------- OUTPUT ------------
        /*cout << "NPP after kmer lookup search" << endl;
        cout << read << endl;
        for (size_t i = 0; i < mergedSeeds.size(); i++) {
                for (size_t j = 0; j < mergedSeeds[i].readFirst; j++)
                        cout << " ";
                cout << "[";
                for (size_t j = mergedSeeds[i].readFirst + 1; j < mergedSeeds[i].readEnd; j++)
                        cout << "-";
                cout << "[";
                cout << " (" << mergedSeeds[i].readFirst << " to " << mergedSeeds[i].readEnd << ")\t";
                for (size_t j = 0; j < mergedSeeds[i].nodeID.size(); j++)
                        cout << mergedSeeds[i].nodeID[j] << "\t";
                cout << endl;
        }*/
        // ----------- OUTPUT ------------
}

int ReadCorrection::correctRead(const string& read,
                                string& bestCorrectedRead,
                                const vector<Seed>& seeds,
                                vector<NodeID>& bestNodeChain)
{
        int bestScore = -read.size();

        for (auto& it : seeds) {
                size_t first = it.readFirst;
                size_t last = it.readEnd;

                // create a npp vector
                vector<NodePosPair> npp(read.length() + 1 - Kmer::getK());
                it.createNodePosition(dbg, npp);

                string correctedRead = read;
                vector<NodeID> nodeChain;
                if (!correctRead(correctedRead, npp, first, last, nodeChain))
                        continue;

                size_t correctedLength = Kmer::getK() - 1 + last - first;
                size_t uncorrectedLength = read.size() - correctedLength;

                int score = alignment.align(read, correctedRead) - uncorrectedLength;

                /*cout << "Seed: " << first << " to " << last << endl;
                alignment.printAlignment(read, correctedRead);
                cout << "Score: " << score << endl;*/

                if (score > bestScore) {
                        bestCorrectedRead = correctedRead;
                        bestScore = score;
                        bestNodeChain = nodeChain;
                }
        }
        if (read == bestCorrectedRead)
                bestScore = read.length();

        return bestScore;
}

void ReadCorrection::correctRead(ReadRecord& record,
                                 AlignmentMetrics& metrics)
{
        bool correctedByMEM = false, readCorrected = false;
        string& read = record.getRead();
        vector<NodeID>& nodeChain = record.getNodeChain();

        // if the read is too short, get out
        if (read.length() < Kmer::getK())
                return;
        vector<Seed> seeds;
        findSeedKmer(read, seeds);

        string bestCorrectedRead;
        vector<NodeID> bestNodeChain;
        int bestScore = correctRead(read, bestCorrectedRead, seeds, bestNodeChain);
        vector<Seed> seedsRC;
        string readRC = record.getRead() ;
        Nucleotide::revCompl(readRC);
        findSeedKmer(readRC, seedsRC);
        string bestCorrectedReadRC;
        vector<NodeID> bestNodeChainRC;
        int bestScoreRC = correctRead(readRC, bestCorrectedReadRC, seedsRC, bestNodeChainRC);
        if (bestScoreRC ==  -read.size() ||  bestScore ==read.size())
                bestScore = -read.size();
        else{
                if (bestScoreRC > bestScore){
                        bestScore = bestScoreRC;
                        Nucleotide::revCompl(readRC);
                        read = readRC;
                        Nucleotide::revCompl(bestCorrectedReadRC);
                        bestCorrectedRead = bestCorrectedReadRC;
                        reverse(bestNodeChainRC.begin(), bestNodeChainRC.end());
                        bestNodeChain.clear();
                        for (auto it :bestNodeChainRC)
                                bestNodeChain.push_back(-it);
                }
        }

        if (bestScore != -read.size() &&  bestScore <= ((int)read.size() / 2)) {
                correctedByMEM = true;
                findSeedMEM(read, seeds);
                bestScore = correctRead(read, bestCorrectedRead, seeds, bestNodeChain);
        }

        /*if (bestCorrectedRead!=""){
                alignment.align(read, bestCorrectedRead);
                cout << "BEST ALIGNMENT: " << bestScore << endl;
                alignment.printAlignment(read, bestCorrectedRead);
                for (auto it:bestNodeChain){
                        cout <<it <<",";
                }
                cout <<endl;
        }*/
        size_t numSubstitutions = 0;
        if (bestScore > ((int)read.size() / 2) &&dbg.validateChain(bestNodeChain)) {
                read = bestCorrectedRead;
                readCorrected = true;
                numSubstitutions = (read.length() - bestScore)/2;
                dbg.validateChain(bestNodeChain);
        }
        metrics.addObservation(readCorrected, correctedByMEM, numSubstitutions);
}

void ReadCorrection::correctChunk(vector<ReadRecord>& readChunk,
                                  AlignmentMetrics& metrics)
{
        ofstream ofs;
        for (ReadRecord& record : readChunk) {
                vector<NodeID> nodeChain;
                correctRead(record, metrics);
        }

        /*cout << readChunk.size() << endl;
        for (size_t i = 0; i < 1000; i++) {
                cout << "================== read " << i << endl;
                correctRead(readChunk[i], metrics);
        }
        cout << "Bye... " << endl;
        exit(0);*/
}

void ReadCorrectionHandler::workerThread(size_t myID, LibraryContainer& libraries,
                                         AlignmentMetrics& metrics)
{
        ReadCorrection readCorrection(dbg, settings, *sa, startpos);

        // local storage of reads
        vector<ReadRecord> myReadBuf;

        // performance counters per thread
        AlignmentMetrics threadMetrics;

        while (true) {
                size_t blockID, recordID;
                bool result = libraries.getRecordChunk(myReadBuf, blockID, recordID);

                readCorrection.correctChunk(myReadBuf, threadMetrics);

                if (result)
                        libraries.commitRecordChunk(myReadBuf, blockID, recordID);
                else
                        break;
        }

        // update the global metrics with the thread info (thread-safe)
        metrics.addMetrics(threadMetrics);
}

void ReadCorrectionHandler::initEssaMEM()
{
        size_t length = 0;
        startpos.clear();
        startpos.reserve(dbg.getNumNodes());
        for (NodeID nodeID = 1; nodeID <= dbg.getNumNodes(); nodeID++) {
                SSNode node = dbg.getSSNode(nodeID);
                assert(node.isValid());
                startpos.push_back(length);
                length += node.getLength() * 2 + 2;
        }

        reference.clear();
        reference.reserve(length);
        for (NodeID nodeID = 1; nodeID < dbg.getNumNodes(); nodeID++) {
                SSNode node = dbg.getSSNode(nodeID);
                if (!node.isValid())
                        continue;

                string thisSequence = node.getSequence();
                reference.append(thisSequence);
                reference.append(">");

                Nucleotide::revCompl(thisSequence);
                reference.append(thisSequence);
                reference.append(">");
        }

        std::vector<std::string> refdescr;
        refdescr.push_back("");

        bool printSubstring = false;
        bool printRevCompForw = false;

        sa = new sparseSA(reference,                            // reference
                          refdescr,
                          startpos,                             // start index for each string
                          false,                                // 4 column format or not
                          settings.getEssaMEMSparsenessFactor(),// ESSA sparseness factor
                          true,                                 // suffixlinks
                          true,                                 // child arrays,
                          true,                                 // kmertable
                          1,                                    // skip parameter
                          10,                                   // kmer size
                          printSubstring,
                          printRevCompForw,
                          false                                 // nucleotides only
                          );
        sa->construct();
}

void ReadCorrectionHandler::doErrorCorrection(LibraryContainer& libraries)
{
        const unsigned int& numThreads = settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        libraries.startIOThreads(settings.getThreadWorkSize(),
                                 10 * settings.getThreadWorkSize() * settings.getNumThreads(),
                                 true);

        // start worker threads
        AlignmentMetrics metrics;
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&ReadCorrectionHandler::workerThread,
                                          this, i, ref(libraries), ref(metrics));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        libraries.joinIOThreads();

        metrics.printStatistics();
}

ReadCorrectionHandler::ReadCorrectionHandler(DBGraph& g, const Settings& s) :
        dbg(g), settings(s), sa(NULL)
{
        Util::startChrono();
        cout << "Creating kmer lookup table... "; cout.flush();
        dbg.buildKmerNPPTable();
        cout << "done (" << Util::stopChronoStr() << ")" << endl;

        Util::startChrono();
        cout << "Building suffix array (sparseness factor: "
             << settings.getEssaMEMSparsenessFactor() << ")..."; cout.flush();
        initEssaMEM();
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
}

ReadCorrectionHandler::~ReadCorrectionHandler()
{
        delete sa;
        dbg.destroyKmerNPPTable();
}
