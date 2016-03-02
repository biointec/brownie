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
#include "graph.h"
#include "readfile/fastafile.h"
extern "C" {
        #include "suffix_tree.h"
}

#include <iomanip>

void DBGraph::writeCytoscapeGraph(const std::string& filename,
                                  NodeID seedNodeID, size_t maxDepth)
{
        // a map containing nodeIDs to handle + their depth
        map<NodeID, size_t> nodeDepth;
        // set of nodes that were handled
        set<NodeID> nodesHandled;

        // if default input values are provided, write the entire graph
        if (seedNodeID == 0) {
                for (NodeID id = -numNodes; id <= numNodes; id++) {
                        if (id == 0)
                                continue;
                        if (!getSSNode(id).isValid())
                                continue;
                        nodeDepth[id] = 0;
                }
        } else {        // else check if seedNode is valid
                if ((abs(seedNodeID) > numNodes) || !getSSNode(seedNodeID).isValid()) {
                        cerr << "WARNING: trying to use an invalid node as a "
                                "seed in writeCytoscapeGraph!" << endl;
                        return;
                }
                nodeDepth[seedNodeID] = 0;
        }

        // map of all nodes in the local graph and its depth
        ofstream ofs((filename + ".arcs").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".arcs" << " for writing" << endl;
        ofs << "Source node\tTarget node\tArc coverage" << endl;

        // A) write all arcs
        while (!nodeDepth.empty()) {
                // get and erase the current node
                map<NodeID, size_t>::iterator e = nodeDepth.begin();
                NodeID thisID = e->first;
                size_t thisDepth = e->second;
                nodeDepth.erase(e);

                // if the node was already handled, skip
                if (nodesHandled.find(thisID) != nodesHandled.end())
                        continue;

                // if we're too far in the graph, stop
                if (thisDepth > maxDepth) {
                        nodesHandled.insert(thisID);
                        continue;
                }

                // write the right arcs
                SSNode thisNode = getSSNode(thisID);
                for (ArcIt it = thisNode.rightBegin(); it != thisNode.rightEnd(); it++) {
                        SSNode rNode = getSSNode(it->getNodeID());
                        if (!rNode.isValid())
                                continue;
                        if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                continue;
                        ofs << thisID << "\t" << it->getNodeID() << "\t" << it->getCoverage() << "\n";
                        if (nodeDepth.find(it->getNodeID()) != nodeDepth.end())
                                continue;
                        nodeDepth[it->getNodeID()] = thisDepth + 1;
                }

                // write the left arcs
                for (ArcIt it = thisNode.leftBegin(); it != thisNode.leftEnd(); it++) {
                        SSNode lNode = getSSNode(it->getNodeID());
                        if (!lNode.isValid())
                                continue;
                        if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                continue;
                        ofs << it->getNodeID() << "\t" << thisID << "\t" << it->getCoverage() << "\n";
                        if (nodeDepth.find(it->getNodeID()) != nodeDepth.end())
                                continue;
                        nodeDepth[it->getNodeID()] = thisDepth + 1;
                }

                // mark this node as handled
                nodesHandled.insert(thisID);

        }
        ofs.close();

        // B) write all nodes
        ofs.open((filename + ".nodes").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".nodes" << " for writing" << endl;
        ofs << "Node ID\tMarginal length\tTrue multiplicity\tEstimated multiplicity"
               "\tKmer coverage\tRead start coverage\tSequence" << "\n";
        for (set<NodeID>::iterator it = nodesHandled.begin(); it != nodesHandled.end(); it++) {
                SSNode node = getSSNode(*it);

               // cout << *it << endl;
                int thisTrueMult = (trueMult.empty()) ? 0 : trueMult[abs(*it)];

                ofs << *it << "\t" << node.getMarginalLength() << "\t"
                    << thisTrueMult << "\t" << "0" << "\t"
                    << node.getAvgKmerCov() << "\t"
                    << "\t" << double(node.getReadStartCov()/node.getMarginalLength())
                    << "\t" << node.getSequence() << "\n";
        }
        ofs.close();
}

void DBGraph::readReferenceGenome()
{
        if (!reference.empty())
                return;

        cout << "Reading reference genome..." << endl;

        FastAFile ass(false);
        ass.open("genome.fasta");
        string read;
        while (ass.getNextRead(read)) {
                read.append(read);
                reference.push_back(read);
                refST.push_back(ST_CreateTree(read.c_str(), read.size()));
                cout << "Adding a reference of length: " << read.size() / 2 << endl;
        }

        cout << "Done reading " << reference.size() << " reference contigs" << endl;
}

void DBGraph::compareToSolution(const string& filename, bool load)
{
        if (load){
                // read the reference genome (genome.fasta) from disk
                readReferenceGenome();
                // try reading the multiplicity file from disk
                trueMult.resize(numNodes + 1);
                cout << "Building a new multiplicity file" << endl;
                ofstream ofs(filename.c_str());
                for (NodeID id = 1; id <= numNodes; id++) {
                        string P = getSSNode(id).getSequence();
                        vector<vector<size_t> > pos, posRC;
                        int ntOcc = findAllTrueOccurences(P, pos, posRC);
                        trueMult[id] = ntOcc;
                        ofs << trueMult[id] << "\n";
                }
                ofs.close();
        }
        size_t sizeCorrect = 0;
        size_t sizeIncorrect = 0;
        size_t sizeTotal = 0;
        size_t numCorrect = 0, numIncorrect = 0, numTotal = 0;

        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;

                numTotal++;
                if (trueMult[id] == 0) {
                        numIncorrect++;
                        sizeIncorrect = sizeIncorrect + node.getMarginalLength();
                } else {
                        numCorrect++;
                        sizeCorrect = sizeCorrect + node.getMarginalLength();
                }

                sizeTotal = sizeTotal + node.getMarginalLength();
        }

        size_t sizeGenome = 0;
        for (size_t i = 0; i < reference.size(); i++)
                sizeGenome += reference[i].size()/2 + 1 - Kmer::getK();

        cout << "\t===== DEBUG: quality report =====" << endl;
        cout << "\tNumber of nodes that exist: " << numCorrect << "/"
             << numTotal << " (" << fixed << setprecision(2)
             << Util::toPercentage(numCorrect, numTotal) << "%) -> ("
             << Util::toPercentage(sizeCorrect, sizeTotal)
             << "% of graph sequence content)" << endl;
        cout << "\tNumber of nodes that do not exist: " << numIncorrect << "/"
             << numTotal << " (" << fixed << setprecision(2)
             << Util::toPercentage(numIncorrect, numTotal) << "%) -> ("
             << Util::toPercentage(sizeIncorrect, sizeTotal)
             << "% of graph sequence content)" << endl;
        cout << "\tThe fraction of the genome that is covered: "
             << sizeCorrect << "/" << sizeGenome << " ("
             << fixed << setprecision(2)
             << Util::toPercentage(sizeCorrect, sizeGenome) << "%)" << endl;
        cout << "\t===== DEBUG: end =====" << endl;
}

void DBGraph::compareToSolution2(const string& filename, bool load)
{
        if (load) {
                // read the reference genome (genome.fasta) from disk
                populateTable();
        }

        FastAFile ass(false);
        ass.open("genome.fasta");

        size_t numKmers = 0, numFound = 0, numBreakpoints = 0;

        string read;
        NodePosPair prev;
        while (ass.getNextRead (read)) {
                for (KmerIt it(read); it.isValid(); it++) {
                        Kmer kmer = it.getKmer();
                        numKmers++;

                        NodePosPair npp = getNodePosPair(kmer);
                        if (!npp.isValid())
                                continue;

                        numFound++;

                        if (prev.isValid()) {
                                if (prev.getNodeID() == npp.getNodeID()) {
                                        if (prev.getOffset()+1 != npp.getOffset()) {
                                                numBreakpoints++;
                                                cout << "Non-adjacent kmer in same node" << endl;
                                        }
                                } else {
                                        if (getSSNode(prev.getNodeID()).getRightArc(npp.getNodeID()) == NULL) {
                                                numBreakpoints++;
                                                cout << "Non-adjacent kmer in different node" << endl;
                                        }
                                }
                        }

                        prev = npp;
                }
        }

        ass.close();

        double fracFound = 100.0*(double)numFound/(double) numKmers;

        cout.precision(4);
        cout << "Validation report: " << endl;
        cout << "\tk-mers in table: " << numFound << "/" << numKmers
             << "(" << fracFound << "%)" << endl;
        cout << "\tnumber of breakpoints: " << numBreakpoints << endl;
}

size_t DBGraph::findAllTrueOccurences(const string& P,
                                      vector<vector<size_t> >& positions,
                                      vector<vector<size_t> >& positionsRC) const
{
        positions.clear();
        positionsRC.clear();
        size_t ntOcc = 0;

        // normal strand
        for (size_t i = 0; i < refST.size(); i++) {
                positions.push_back(vector<size_t>());
                SUFFIX_TREE *st = (SUFFIX_TREE*)refST[i];

                int *pos;
                char *c_str = const_cast<char*>(P.c_str());
                int numOcc = ST_FindAllSubstrings(st, c_str, P.size(), &pos);

                for (int j = 0; j < numOcc; j++) {
                        if ((pos[j]-1) >= (int)reference[i].size()/2)
                                continue;
                        positions.back().push_back(pos[j]-1);
                        ntOcc++;
                }

                if (numOcc > 0)
                        free(pos);
        }

        // calculate the reverse complement, if it is the same, get out now
        string RC = Nucleotide::getRevCompl(P);

        if (RC == P)
                return ntOcc;

        // opposite strand
        for (size_t i = 0; i < refST.size(); i++) {
                positionsRC.push_back(vector<size_t>());
                SUFFIX_TREE *st = (SUFFIX_TREE*)refST[i];

                int *pos;
                char *c_str = const_cast<char*>(RC.c_str());
                int numOcc = ST_FindAllSubstrings(st, c_str, RC.size(), &pos);

                for (int j = 0; j < numOcc; j++) {
                        if ((pos[j]-1) >= (int)reference[i].size()/2)
                                continue;
                        positionsRC.back().push_back(pos[j]-1);
                        ntOcc++;
                }

                if (numOcc > 0)
                        free(pos);
        }

        return ntOcc;
}
