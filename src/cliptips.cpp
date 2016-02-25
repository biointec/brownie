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

#include "graph.h"
#include "settings.h"

#include <queue>

using namespace std;

// ============================================================================
// GRAPH STATISTICS
// ============================================================================

std::ostream &operator<<(std::ostream &out, const vector<NodeID> &path)
{
        for (auto& it : path)
                out << it << " ";

        return out;
}

// ============================================================================
// CLIP TIPS
// ============================================================================

bool DBGraph::clipTips(double covCutoff, size_t maxMargLength)
{
        cout << endl << "=================== Removing tips ===================" << endl;

#ifdef DEBUG
        size_t tp=0, tn=0, fp=0,fn=0;
        size_t tps=0, tns=0, fps=0,fns=0;
        size_t tpj=0, tnj=0, fpj=0,fnj=0;
#endif

        size_t numDeleted = 0, numTotal = 0;
        cout << "Cut-off value for removing tips is: " << covCutoff << endl;

        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;
                numTotal++;

                // check for dead ends
                bool leftDE = (node.getNumLeftArcs() == 0);
                bool rightDE = (node.getNumRightArcs() == 0);
                if (!leftDE && !rightDE)
                        continue;

                SSNode startNode = (rightDE) ? getSSNode(-id) : getSSNode(id);
                bool isolated = rightDE && leftDE;
                bool joinedTip = startNode.getNumRightArcs() > 1;
                bool remove = ((startNode.getAvgKmerCov() <= covCutoff) &&
                               (startNode.getMarginalLength() <= maxMargLength));
                if (remove) {
                        removeNode(startNode.getNodeID());
                        numDeleted++;
                }

#ifdef DEBUG
                if (remove) {
                        if (trueMult.size()>0&& trueMult[id] > 0) {
                                if (isolated)
                                        fps++;
                                else if(joinedTip)
                                        fpj++;
                                else
                                        fp++;
                        } else {
                                if(isolated)
                                        tps++;
                                else if(joinedTip)
                                        tpj++;
                                else
                                        tp++;
                        }
                } else {
                        if (trueMult.size()>0&& trueMult[id] > 0) {
                                if(isolated)
                                        tns++;
                                else if(joinedTip)
                                        tnj++;
                                else
                                        tn++;
                        } else {
                                if(isolated)
                                        fns++;
                                else if(joinedTip)
                                        fnj++;
                                else
                                        fn++;
                        }
                }
#endif
        }

        cout << "Clipped " << numDeleted << "/" << numTotal << " nodes" << endl;
#ifdef DEBUG
        cout << "\t===== DEBUG: tip clipping report =====" << endl;
        cout << "\tIsolated TP: " << tps << "\tTN: "<< tns << "\tFP: " << fps << "\tFN: "<< fns << endl;
        cout << "\tSensitivity: " << 100.0 * Util::getSensitivity(tps, fns) << "%" << endl;
        cout << "\tSpecificity: " << 100.0 * Util::getSpecificity(tns, fps) << "%" << endl;
        cout << "\t****************************************" << endl;
        cout << "\tJoined TP: " << tpj << "\tTN: " << tnj << "\tFP: " << fpj << "\tFN: " << fnj << endl;
        cout << "\tSensitivity: " << 100.0 * Util::getSensitivity(tpj, fnj) << "%" << endl;
        cout << "\tSpecificity: " << 100.0 * Util::getSpecificity(tnj, fpj) << "%" << endl;
        cout << "\t****************************************" << endl;
        cout << "\tTip TP: " << tp << "\tTN: " << tn << "\tFP: " << fp << "\tFN: " << fn << endl;
        cout << "\tSensitivity: " << 100.0 * Util::getSensitivity(tp, fn) << "%" << endl;
        cout << "\tSpecificity: " << 100.0 * Util::getSpecificity(tn, fp) << "%" << endl;
        cout << "\t===== DEBUG: end =====" << endl;
#endif

        return  numDeleted > 0;
}

bool DBGraph::canonicalBubble(double covCutoff)
{
        cout << endl << "=================== Canonical Bubble ===================" << endl;

#ifdef DEBUG
        size_t tp=0, tn=0, fp=0,fn=0;
#endif

        size_t numDeleted = 0, numTotal = 0;
        cout << "Cut-off value for removing bubbles is: " << covCutoff << endl;

        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;
                numTotal++;

                // check for nodes of length k
                if (node.getMarginalLength() != Kmer::getK())
                        continue;
                if (node.getNumLeftArcs() != 1)
                        continue;
                if (node.getNumRightArcs() != 1)
                        continue;

                bool remove = (node.getAvgKmerCov() <= covCutoff);
                if (remove) {
                        removeNode(id);
                        numDeleted++;
                }

#ifdef DEBUG
                if (remove) {
                        if (trueMult.size()>0&& trueMult[id] > 0) {
                                fp++;
                        } else {
                                tp++;
                        }
                } else {
                        if (trueMult.size()>0&& trueMult[id] > 0) {
                                tn++;
                        } else {
                                fn++;
                        }
                }
#endif
        }

        cout << "Clipped " << numDeleted << "/" << numTotal << " nodes" << endl;
#ifdef DEBUG
        cout << "\t===== DEBUG: canonical bubble report =====" << endl;
        cout << "\tCanonical bubble TP: " << tp << "\tTN: " << tn << "\tFP: " << fp << "\tFN: " << fn << endl;
        cout << "\tSensitivity: " << 100.0 * Util::getSensitivity(tp, fn) << "%" << endl;
        cout << "\tSpecificity: " << 100.0 * Util::getSpecificity(tn, fp) << "%" << endl;
        cout << "\t===== DEBUG: end =====" << endl;
#endif

        return  numDeleted > 0;
}

bool DBGraph::concatenateNodes()
{
        size_t numConcatentations = 0;
        size_t numIncorrectConcatenations=0;

        for (NodeID lID = -numNodes; lID <= numNodes; lID++) {
                if (lID == 0)
                        continue;
                SSNode left = getSSNode(lID);
                if (!left.isValid())
                        continue;
                if (left.getNumRightArcs() != 1)
                        continue;

                NodeID rID = left.rightBegin()->getNodeID();

                // don't merge palindromic repeats
                if (rID == -lID)
                        continue;

                SSNode right = getSSNode (rID);
                if (right.getNumLeftArcs() != 1)
                        continue;

#ifdef DEBUG
                if (!trueMult.empty()) {
                        if (trueMult[abs(lID)] != trueMult[abs(rID)]) {
                                numIncorrectConcatenations++;
                                trueMult[abs(lID)]=0;
                                trueMult[abs(rID)]=0;
                        }
                }
#endif

                left.deleteRightArc(rID);
                right.deleteLeftArc(lID);
                left.inheritRightArcs(right);
                left.setKmerCov(left.getKmerCov()+right.getKmerCov());
                left.setReadStartCov(left.getReadStartCov()+right.getReadStartCov());
                right.invalidate();
                numConcatentations++;

                deque<SSNode> deq;
                deq.push_back(left);
                deq.push_back(right);

                string str;
                convertNodesToString(deq, str);

                left.setSequence(str);
                lID--;
        }

        if (numConcatentations > 0)
                cout << "Concatenated " << numConcatentations << " nodes" << endl;
#ifdef DEBUG
        cout << "\t===== DEBUG: tip clipping report =====" << endl;
        cout << "\t" << "Number of incorrect connections: " << numIncorrectConcatenations << endl;
        cout << "\t===== DEBUG: end =====" << endl;
#endif

        return (numConcatentations > 0);
}

class PathInfo {
public:
        NodeID nodeID;
        size_t length;

        /**
         * Default constructor
         */
        PathInfo(NodeID nodeID, size_t length) : nodeID(nodeID), length(length) {};
};

struct PathInfoComp {

        /**
         * Compare two paths (by length)
         */
        bool operator()(const PathInfo& f, const PathInfo& s) {
                return f.length >= s.length;
        }
};

vector<NodeID> DBGraph::getPath(NodeID dstID, const vector<NodeID>& prevNode) const
{
        vector<NodeID> path;
        path.push_back(dstID);

        while (true) {
                dstID = prevNode[dstID + numNodes];
                if (dstID == 0)
                        break;
                path.push_back(dstID);
        }

        reverse(path.begin(), path.end());
        return path;
}

void DBGraph::getUniquePath(const std::vector<NodeID>& path,
                           size_t& first, size_t& last)
{
        // sanity check is necessary
        if (path.size() <= 2) {
                first = 1; last = 0;
                return;
        }

        // find the first node that can be deleted
        first = 1;
        for (size_t i = path.size() - 2; i >= 1; i--) {
                if (getSSNode(path[i]).getNumRightArcs() == 1)
                        continue;
                first = i+1;
                break;
        }

        // find the last node that can be deleted
        last = path.size() - 2;
        for (size_t i = 1; i <= path.size() - 2; i++) {
                if (getSSNode(path[i]).getNumLeftArcs() == 1)
                        continue;
                last = i-1;
                break;
        }
}

double DBGraph::getPathAvgKmerCov(const vector<NodeID>& path)
{
        double totKmer = 0, totLen = 0;
        for (auto it : path) {
                totKmer += getSSNode(it).getKmerCov();
                totLen += getSSNode(it).getMarginalLength();
        }

        return totKmer / totLen;
}

void DBGraph::removePath(const vector<NodeID>& pathA)
{
        for (auto it : pathA) {
                //cout << "\tDeleted node " << it << endl;
                if (trueMult[abs(it)] > 0) {
                        cout << "Incorrectly DELETED NODE: " << it << endl;
                        //writeCytoscapeGraph(settings.getTempDirectory() + "cytIncorrect", it, 10);
                        //exit(0);
                }
                removeNode(it);
        }
}

bool DBGraph::handleParallelPaths(const vector<NodeID>& pathA,
                                  const vector<NodeID>& pathB,
                                  double covCutoff)
{
        size_t firstA, lastA;
        getUniquePath(pathA, firstA, lastA);
        vector<NodeID> subPathA;
        if (firstA <= lastA)
                subPathA = vector<NodeID>(pathA.begin() + firstA,
                                          pathA.begin() + lastA + 1);
        double covSubPathA = subPathA.empty() ? covCutoff + 1 :
                getPathAvgKmerCov(subPathA);

        size_t firstB, lastB;
        getUniquePath(pathB, firstB, lastB);
        vector<NodeID> subPathB;
        if (firstB <= lastB)
                subPathB = vector<NodeID>(pathB.begin() + firstB,
                                          pathB.begin() + lastB + 1);
        double covSubPathB = subPathB.empty() ? covCutoff + 1 :
                getPathAvgKmerCov(subPathB);

        const vector<NodeID>& lowCovPath = (covSubPathA < covSubPathB) ?
                subPathA : subPathB;
        const double& lowCov = (covSubPathA < covSubPathB) ?
                covSubPathA : covSubPathB;

       /* cout << "Path A" << endl;
        cout << pathA << endl;
        cout << "First: " << firstA << ", last: " << lastA << endl;
        cout << subPathA << endl;
        cout << covSubPathA << endl;

        cout << " -- " << endl;

        cout << "Path B" << endl;
        cout << pathB << endl;
        cout << "First: " << firstB << ", last: " << lastB << endl;
        cout << subPathB << endl;
        cout << covSubPathB << endl;*/

        if ((lowCov <= covCutoff) && !lowCovPath.empty()) {
                removePath(lowCovPath);
                return true;
        }

        // Remove final arc?
        double covFinalArcA = ((pathA.size() >= 2) && (lastA == pathA.size() - 2)) ?
                getSSNode(pathA[lastA]).getRightArc(pathA[lastA+1])->getCoverage() : covCutoff + 1;
        double covFinalArcB = ((pathB.size() >= 2) && (lastB == pathB.size() - 2)) ?
                getSSNode(pathB[lastB]).getRightArc(pathB[lastB+1])->getCoverage() : covCutoff + 1;
        if (min(covFinalArcA, covFinalArcB) <= covCutoff) {
                cout << "Removing final arcs" << endl;
                if (covFinalArcA < covFinalArcB) {
                        getSSNode(pathA[lastA]).deleteRightArc(pathA[lastA+1]);
                        getSSNode(pathA[lastA+1]).deleteLeftArc(pathA[lastA]);
                } else {
                        getSSNode(pathB[lastB]).deleteRightArc(pathB[lastB+1]);
                        getSSNode(pathB[lastB+1]).deleteLeftArc(pathB[lastB]);
                }
                return true;
        }

        // Remove first arc?
        double covFirstArcA = ((pathA.size() >= 2) && (firstA == 1)) ?
                getSSNode(pathA[0]).getRightArc(pathA[1])->getCoverage() : covCutoff + 1;
        double covFirstArcB = ((pathB.size() >= 2) && (firstB == 1)) ?
                getSSNode(pathB[0]).getRightArc(pathB[1])->getCoverage() : covCutoff + 1;
        if (min(covFirstArcA, covFirstArcB) <= covCutoff) {
                cout << "Removing first arcs" << endl;
                if (covFirstArcA < covFirstArcB) {
                        getSSNode(pathA[0]).deleteRightArc(pathA[1]);
                        getSSNode(pathA[1]).deleteLeftArc(pathA[0]);
                } else {
                        getSSNode(pathB[0]).deleteRightArc(pathB[1]);
                        getSSNode(pathB[1]).deleteLeftArc(pathB[0]);
                }
                return true;
        }

        return false;
}

bool DBGraph::bubbleDetection(NodeID srcID, vector<NodeID>& visited,
                              vector<NodeID>& prevNode,
                              vector<NodeID>& nodeColor,
                              double covCutoff,
                              size_t maxMargLength, size_t maxNodesVisited)
{
        priority_queue<PathInfo, vector<PathInfo>, PathInfoComp> heap;
        heap.push(PathInfo(srcID, 0));

        bool returnValue = false;

        while(!heap.empty()) {
                PathInfo currTop = heap.top();
                heap.pop();
                NodeID currID = currTop.nodeID;
                size_t currLength = currTop.length;
                SSNode curr = getSSNode(currID);
                if (!curr.isValid())
                        continue;
                for (ArcIt it = curr.rightBegin(); it != curr.rightEnd(); it++ ) {
                        NodeID nextID = it->getNodeID();
                        SSNode next = getSSNode(nextID);

                        // do we encounter a node previously encountered?
                        if ((prevNode[nextID + numNodes] != 0) || (nextID == srcID)) {
                                if (nodeColor[nextID + numNodes] == nodeColor[currID + numNodes])
                                        continue;

                                vector<NodeID> pathA = getPath(currID, prevNode);
                                pathA.push_back(nextID);
                                vector<NodeID> pathB = getPath(nextID, prevNode);

                                // if at least one node was deleted: get out of here!
                                if (handleParallelPaths(pathA, pathB, covCutoff)) {
                                        returnValue = true;
                                        goto exitRoutine;
                                }
                        } else {
                                prevNode[nextID + numNodes] = currID;
                                nodeColor[nextID + numNodes] = (currID == srcID) ?
                                        nextID : nodeColor[currID + numNodes];
                                visited.push_back(nextID);

                                size_t nextLength = currLength + next.getMarginalLength();
                                if (nextLength > maxMargLength)
                                        continue;
                                if (visited.size() > maxNodesVisited)
                                        continue;
                                PathInfo nextTop(nextID, nextLength);
                                heap.push(nextTop);
                        }
                }
        }

        // label definition to break out of nested loops
        exitRoutine:

        for (auto it : visited) {
                prevNode[it + numNodes] = 0;
                nodeColor[it + numNodes] = 0;
        }
        visited.clear();

        return returnValue;
}

bool DBGraph::bubbleDetection(double covCutoff, size_t maxMargLength)
{
        cout << endl << "=================== Removing bubbles ===================" << endl;

        vector<NodeID> visited;
        vector<NodeID> prevNode(2*numNodes+1, 0);
        vector<NodeID> nodeColor(2*numNodes+1, 0);

        bool returnValue = false;

        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                if (!getSSNode(id).isValid())
                        continue;
                if (getSSNode(id).getNumRightArcs() < 2)
                        continue;

                while (bubbleDetection(id, visited, prevNode, nodeColor,
                                       covCutoff, maxMargLength,
                                       settings.getBubbleDFSNodeLimit()))
                        returnValue = true;
        }

        return returnValue;
}

bool DBGraph::bubbleDetection(NodeID id, double covCutoff, size_t maxMargLength)
{
        vector<NodeID> visited;
        vector<NodeID> prevNode(2*numNodes+1, 0);
        vector<NodeID> nodeColor(2*numNodes+1, 0);

        bool returnValue = false;
        while (bubbleDetection(id, visited, prevNode, nodeColor,
                               covCutoff, maxMargLength,
                               settings.getBubbleDFSNodeLimit()))
                returnValue = true;

        return returnValue;
}
