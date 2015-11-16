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

void ReadCorrectionJan::findNPPSlow(string& read, vector<NodePosPair>& npp)
{
        for (KmerIt it(read); it.isValid(); it++) {
                Kmer kmer = it.getKmer();
                NodePosPair result = dbg.getNodePosPair(kmer);
                npp[it.getOffset()] = result;
        }
}

void ReadCorrectionJan::findNPPFast(string& read, vector<NodePosPair>& nppv)
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

void ReadCorrectionJan::checkConsistency(vector<NodePosPair>& nppv,
                                         vector<bool>& consistency)
{
        size_t prev = nppv.size();
        for (size_t i = 0; i < nppv.size(); i++) {
                if (!nppv[i].isValid())
                        continue;

                // is it the first time we encounter a valid npp?
                if (prev == nppv.size()) {
                        prev = i;
                        continue;
                }

                const SSNode thisNode = dbg.getSSNode(nppv[i].getNodeID());
                const SSNode prevNode = dbg.getSSNode(nppv[prev].getNodeID());

                // no, check for check for consistency
                if (thisNode.getNodeID() == prevNode.getNodeID()) {
                        if ((nppv[i].getOffset() - nppv[prev].getOffset()) == (i - prev)) {
                                for (size_t j = prev; j < i; j++)
                                        consistency[j] = true;
                                for (size_t j = prev+1; j < i; j++)
                                        nppv[j] = NodePosPair(prevNode.getNodeID(), nppv[prev].getOffset()+j-prev);
                        }
                } else {                // we're in different nodes
                        if (prevNode.getRightArc(thisNode.getNodeID()) != NULL) {
                                size_t thisPos = prevNode.getMarginalLength() + nppv[i].getOffset();
                                if ((thisPos - nppv[prev].getOffset()) == (i - prev)) {
                                        for (size_t j = prev; j < i; j++)
                                                consistency[j] = true;
                                        int k = 0;
                                        for (size_t j = prev+1; j < i; j++)
                                                if (nppv[prev].getOffset()+j-prev < prevNode.getMarginalLength())
                                                        nppv[j] = NodePosPair(prevNode.getNodeID(), nppv[prev].getOffset()+j-prev);
                                                else
                                                        nppv[j] = NodePosPair(thisNode.getNodeID(), k++);
                                }
                        }
                }

                prev = i;
        }
}

void ReadCorrectionJan::findSeed(const vector<NodePosPair>& nppv,
                                 const vector< bool >& consistency,
                                 vector<pair<size_t, size_t> >& seeds)
{
        bool intervalOpen = false; size_t thisFirst;
        for (size_t i = 0; i < nppv.size(); i++) {

                if (!intervalOpen) {
                        if (nppv[i].isValid()) {        // open event
                                thisFirst = i;
                                intervalOpen = true;
                        }
                } else {
                        if (!nppv[i].isValid() || !consistency[i-1]) {  // close event
                                seeds.push_back(pair<size_t, size_t>(thisFirst, i));
                                intervalOpen = false;
                        }
                }
        }

        if (intervalOpen)
                seeds.push_back(pair<size_t, size_t>(thisFirst, nppv.size()));
}

void ReadCorrectionJan::recSearch(NodeID curr, string& read, vector<NodePosPair>& npp,
                                  size_t currReadPos, size_t& counter,
                                  int currScore, int& bestScore, size_t& seedLast)
{
        const SSNode node = dbg.getSSNode(curr);

        counter++;
        if (counter > 1000)
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
        // delete unreliable right nodes: FIXME generalize !!
//         while (seedLast > seedFirst + 1) {
//                 if (npp[seedLast - 1].getOffset() == 0 || npp[seedLast - 1].getOffset() == 1 /* || npp[seedLast - 1].getOffset() == 2 || npp[seedLast - 1].getOffset() == 3*/) {
//                         npp[seedLast - 1] = NodePosPair(0, 0);
//                         seedLast--;
//                 } else
//                         break;
//         }

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
                curr += strLen;
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

bool ReadCorrectionJan::correctRead(ReadRecord& record)
{
        string& read = record.getRead();

        // if the read is too short, get out
        if (read.length() < Kmer::getK())
                return false;

        /////////////// INITAL RUN

        vector<NodePosPair> npp(read.length() + 1 - Kmer::getK());
        findNPPFast(read, npp);

        /*cout << "NPP after initial search" << endl;
        for (size_t i = 0; i < Kmer::getK() - 1; i++)
                cout << " ";
        for (size_t i = 0; i < npp.size(); i++) {
                if (npp[i].isValid())
                        cout << i << ":" << npp[i].getNodeID() << ":" << npp[i].getOffset() << " ";
                        //cout << "|";
                else
                        cout << "*";
        }
        cout << endl;*/

        /////////////// CONSISTENCY

        vector<bool> consistency(npp.size() - 1, false);
        checkConsistency(npp, consistency);

        /*for (size_t i = 0; i < Kmer::getK() - 1; i++)
                cout << " ";
        for (auto it : consistency)
                cout << it;
        cout << endl;*/

        /////////////// SEED FINDING

        vector<pair<size_t, size_t> > seeds;
        findSeed(npp, consistency, seeds);

        // if no seed is found: try the MEM approach
        if (seeds.empty())
                return false;

        int bestScore = -read.size();
        string bestCorrectedRead = read;

        for (auto& it : seeds) {
                size_t first = it.first;
                size_t last = it.second;

                //cout << "Seed: " << it.first << ", " << it.second << endl;
                string correctedRead = read;
                correctRead(correctedRead, npp, first, last);

                size_t correctedLength = Kmer::getK() - 1 + last - first;
                size_t uncorrectedLength = read.size() - correctedLength;

                int score = alignment.align(read, correctedRead) - uncorrectedLength;
                /*alignment.printAlignment(read, correctedRead);
                cout << "Score: " << score << endl;*/

                if (score > bestScore) {
                        bestCorrectedRead = correctedRead;
                        bestScore = score;
                }
        }

        if (bestScore > ((int)read.size() / 2))
                read = bestCorrectedRead;

        //cout << read << endl;

        return true;
}

void ReadCorrectionJan::correctChunk(vector<ReadRecord>& readChunk)
{
        for (auto& it : readChunk)
                correctRead(it);

        /*cout << readChunk.size() << endl;
        for (size_t i = 2123; i < 2124; i++) {
                cout << "================== read " << i << endl;
                correctRead(readChunk[i]);
        }

        cout << "Bye... " << endl;
        exit(0);*/
}

void ReadCorrectionHandler::workerThread(size_t myID, LibraryContainer& libraries)
{
        ReadCorrectionJan readCorrection(dbg, settings);

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

void ReadCorrectionHandler::doErrorCorrection(LibraryContainer& libraries)
{
        const unsigned int& numThreads = settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        libraries.startIOThreads(settings.getThreadWorkSize(),
                                 4 * settings.getThreadWorkSize() * settings.getNumThreads(),
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
        dbg(g), settings(s)
{
        Util::startChrono();
        cout << "Creating kmer lookup table... "; cout.flush();
        dbg.populateTable();
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
}
