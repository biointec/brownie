/***************************************************************************
 *   Copyright (C) 2014, 2015 Jan Fostier (jan.fostier@intec.ugent.be)     *
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
#include "library.h"
#include "readcorrection.h"
#include "kmernode.h"

using namespace std;

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

void ReadCorrectionJan::revCompl(vector< NodePosPair >& npp)
{
        reverse(npp.begin(), npp.end());

        for (size_t i = 0; i < npp.size(); i++) {
                if (!npp[i].isValid())
                        continue;

                NodeID nodeID = npp[i].getNodeID();
                size_t pos = npp[i].getOffset();
                const SSNode node = dbg.getSSNode(nodeID);
                npp[i] = NodePosPair(-nodeID, node.getMarginalLength() - 1 - pos);
        }
}

void ReadCorrectionJan::findNPPSlow(const string& read, vector<NodePosPair>& npp)
{
        for (KmerIt it(read); it.isValid(); it++) {
                Kmer kmer = it.getKmer();
                NodePosPair result = dbg.getNodePosPair(kmer);
                npp[it.getOffset()] = result;
        }
}

void ReadCorrectionJan::findNPPFast(const string& read, vector<NodePosPair>& nppv)
{
        for (KmerIt it(read); it.isValid(); it++) {
                Kmer kmer = it.getKmer();
                NodePosPair npp = dbg.getNodePosPair(kmer);
                nppv[it.getOffset()] = npp;

                if (!npp.isValid())
                        continue;

                NodeID nodeID = npp.getNodeID();
                const SSNode node = dbg.getSSNode(nodeID);

                size_t readPos = it.getOffset() + Kmer::getK();
                size_t nodePos = npp.getOffset() + Kmer::getK();

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

void ReadCorrectionJan::extractSeeds(const vector<NodePosPair>& nppv,
                                     vector<Seed>& seeds)
{
        size_t prev = nppv.size();
        for (size_t i = 0; i < nppv.size(); i++) {
                if (!nppv[i].isValid())
                        continue;

                // is it the first time we encounter a valid npp?
                if (prev == nppv.size()) {
                        prev = i;
                        seeds.push_back(Seed(nppv[i].getNodeID(), nppv[i].getOffset(), i, i + 1));
                        continue;
                }

                bool consistent = false;
                const SSNode thisNode = dbg.getSSNode(nppv[i].getNodeID());
                const SSNode prevNode = dbg.getSSNode(nppv[prev].getNodeID());

                // no, check for check for consistency
                if (thisNode.getNodeID() == prevNode.getNodeID()) {
                        if ((nppv[i].getOffset() - nppv[prev].getOffset()) == (i - prev))
                                consistent = true;
                } else {                // we're in different nodes
                        if (prevNode.getRightArc(thisNode.getNodeID()) != NULL) {
                                size_t thisPos = prevNode.getMarginalLength() + nppv[i].getOffset();
                                if ((thisPos - nppv[prev].getOffset()) == (i - prev)) {
                                        seeds.back().nodeID.push_back(thisNode.getNodeID());
                                        consistent = true;
                                }
                        }
                }

                prev = i;
                if (!consistent) {
                        seeds.push_back(Seed(nppv[i].getNodeID(), nppv[i].getOffset(), i, i + 1));
                } else {
                        seeds.back().readEnd = i + 1;
                }
        }
}

void ReadCorrectionJan::recSearch(NodeID curr, string& read, vector<NodePosPair>& npp,
                                  size_t currReadPos, size_t& counter,
                                  int currScore, int& bestScore, size_t& seedLast)
{
        const SSNode node = dbg.getSSNode(curr);

        counter++;
        if (counter > settings.getMaxDepth())
                return;

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
                float nextRelScore = (float)thisScore / (float)nextNode.getMarginalLength();

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
                if (nextScore > bestScore)
                        bestScore = nextScore;

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
                        recSearch(nextID, read, npp, nextReadPos, counter, nextScore, bestScore, seedLast);

                // if the best score has been updated in this branch...
                if (bestScore <= prevBestScore)
                        continue;

                // ... save the npp
                for (size_t i = 0; i < OLSize; i++)
                        npp[currReadPos+i] = NodePosPair(nextID, i);

                if (currReadPos + OLSize > seedLast)
                        seedLast = currReadPos + OLSize;
        }
}

void ReadCorrectionJan::extendSeed(string& read, vector<NodePosPair>& npp,
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
               (npp[seedLast - 1].getOffset() + 1 < node.getMarginalLength())) {
                npp[seedLast] = NodePosPair(node.getNodeID(), npp[seedLast-1].getOffset() + 1);
                seedLast++;
        }

        // try to find the best right path
        //if (seedLast < getMarginalLength(read))
        //      cout << read << endl;
        size_t counter = 0; int bestScore = -(getMarginalLength(read) - seedLast);
        if (seedLast < getMarginalLength(read))
                recSearch(node.getNodeID(), read, npp, seedLast, counter, 0, bestScore, seedLast);
}

void ReadCorrectionJan::applyReadCorrection(string& read,
                                            const vector<NodePosPair>& npp,
                                            size_t seedFirst, size_t seedLast)
{
        // correct the seed
        size_t curr = seedFirst;
        while (curr < seedLast) {
                SSNode node = dbg.getSSNode(npp[curr].getNodeID());
                size_t strLen = min(node.getLength() - npp[curr].getOffset(), read.size() - curr);
                string str = node.substr(npp[curr].getOffset(), strLen);

                read.replace(curr, strLen, str);
                curr += (strLen - Kmer::getK() + 1);
        }
}

void ReadCorrectionJan::correctRead(string& read, vector<NodePosPair> npp,
                                    size_t& first, size_t& last)
{
        // extend to the right
        extendSeed(read, npp, first, last);

        // reverse complement all data
        Nucleotide::revCompl(read);
        first = getMarginalLength(read) - first;
        last = getMarginalLength(read) - last;
        swap(first, last);
        revCompl(npp);

        // extend to the right (which used to be left)
        extendSeed(read, npp, first, last);

        Nucleotide::revCompl(read);
        first = getMarginalLength(read) - first;
        last = getMarginalLength(read) - last;
        swap(first, last);
        revCompl(npp);

        // reverse complement all data
        string original = read;
        applyReadCorrection(read, npp, first, last);
}

void ReadCorrectionJan::findSeedKmer(const std::string& read,
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

bool sortByLength(const Seed& a, const Seed& b) {
        return ((a.readEnd-a.readFirst) > (b.readEnd-b.readFirst));
}

void ReadCorrectionJan::findSeedMEM(const string& read,
                                    vector<Seed>& mergedSeeds)
{
        vector<match_t> matches;

        int memSize = Kmer::getK() - 1;
        while (matches.size() < 100&& memSize>5) {
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
        sort(mergedSeeds.begin(), mergedSeeds.end(), sortByLength);

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

int ReadCorrectionJan::correctRead(const string& read,
                                   string& bestCorrectedRead,
                                   const vector<Seed>& seeds)
{
        int bestScore = -read.size();

        for (auto& it : seeds) {
                size_t first = it.readFirst;
                size_t last = it.readEnd;

                // create a npp vector
                vector<NodePosPair> npp(read.length() + 1 - Kmer::getK());
                it.createNodePosition(dbg, npp);

                string correctedRead = read;
                correctRead(correctedRead, npp, first, last);

                size_t correctedLength = Kmer::getK() - 1 + last - first;
                size_t uncorrectedLength = read.size() - correctedLength;

                int score = alignment.align(read, correctedRead) - uncorrectedLength;

                /*cout << "Seed: " << first << " to " << last << endl;
                alignment.printAlignment(read, correctedRead);
                cout << "Score: " << score << endl;*/

                if (score > bestScore) {
                        bestCorrectedRead = correctedRead;
                        bestScore = score;
                }
        }

        return bestScore;
}

void ReadCorrectionJan::correctRead(ReadRecord& record)
{
        string& read = record.getRead();

        // if the read is too short, get out
        if (read.length() < Kmer::getK())
                return;

        vector<Seed> seeds;
        findSeedKmer(read, seeds);

        string bestCorrectedRead;
        int bestScore = correctRead(read, bestCorrectedRead, seeds);

        if (bestScore <= ((int)read.size() / 2)) {
                findSeedMEM(read, seeds);
                bestScore = correctRead(read, bestCorrectedRead, seeds);
        }

        /*alignment.align(read, bestCorrectedRead);
        cout << "BEST ALIGNMENT: " << bestScore << endl;
        alignment.printAlignment(read, bestCorrectedRead);*/

        if (bestScore > ((int)read.size() / 2))
                read = bestCorrectedRead;
}

void ReadCorrectionJan::correctChunk(vector<ReadRecord>& readChunk)
{
        for (auto& it : readChunk)
                correctRead(it);

        /*cout << readChunk.size() << endl;
        for (size_t i = 269; i < 270; i++) {
                cout << "================== read " << i << endl;
                correctRead(readChunk[i]);
        }

        cout << "Bye... " << endl;
        exit(0);*/
}

void ReadCorrectionHandler::workerThread(size_t myID, LibraryContainer& libraries)
{
        ReadCorrectionJan readCorrection(dbg, settings, *sa, startpos);

        // local storage of reads
        vector<ReadRecord> myReadBuf;

        while (true) {
                size_t blockID, recordID;
                bool result = libraries.getRecordChunk(myReadBuf, blockID, recordID);

                readCorrection.correctChunk(myReadBuf);

                if (result)
                        libraries.commitRecordChunk(myReadBuf, blockID, recordID);
                else
                        break;
        }
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

        sa = new sparseSA(reference,             //reference
                          refdescr,             //
                          startpos,             //
                          false,                //4 column format or not
                          settings.getEssaFactor(),                    //sparseness
                          false,                //suffixlinks
                          true,                 //child arrays
                          1,                    //skip parameter
                          printSubstring,       //
                          printRevCompForw);
        sa->construct();
}

void ReadCorrectionHandler::doErrorCorrection(LibraryContainer& libraries)
{
        const unsigned int& numThreads = settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        cout << "Building suffix array..."; cout.flush();
        initEssaMEM();
        cout << "done" << endl;

        libraries.startIOThreads(settings.getThreadWorkSize(),
                                 10 * settings.getThreadWorkSize() * settings.getNumThreads(),
                                 true);

        // start worker threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&ReadCorrectionHandler::workerThread,
                                          this, i, ref(libraries));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        libraries.joinIOThreads();
}

ReadCorrectionHandler::ReadCorrectionHandler(DBGraph& g, const Settings& s) :
        dbg(g), settings(s), sa(NULL)
{
        Util::startChrono();
        cout << "Creating kmer lookup table... "; cout.flush();
        dbg.populateTable();
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
}
ReadCorrectionHandler::~ReadCorrectionHandler()
{
        delete sa;
        dbg.depopulateTable();
}
