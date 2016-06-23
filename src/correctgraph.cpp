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
#include "correctgraph.h"

#include <queue>
#include <map>

using namespace std;

// ============================================================================
// PRIVATE CORRECTGRAPH.CPP (STAGE 4 ROUTINES)
// ============================================================================

void DBGraph::removeNode(NodeID nodeID)
{
#ifdef DEBUG
        if (trueMult[abs(nodeID)] > 0)
                cout << "\tERROR removing node " << nodeID << endl;
#endif

        SSNode node = getSSNode(nodeID);
        for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
                SSNode leftNode = getSSNode(it->getNodeID());
                if (leftNode.getNodeID()==-node.getNodeID())
                        continue;
                bool result = leftNode.deleteRightArc (node.getNodeID());
                assert(result);
        }

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                SSNode rightNode = getSSNode ( it->getNodeID() );
                if (rightNode.getNodeID()==-node.getNodeID())
                        continue;
                bool result = rightNode.deleteLeftArc(node.getNodeID());
                assert(result);
        }

        node.deleteAllRightArcs();
        node.deleteAllLeftArcs();
        node.invalidate();
}

void DBGraph::detachNode(NodeID leftID, NodeID rightID)
{
#ifdef DEBUG
        if (getSSNode(leftID).getRightArc(rightID)->getTrueArc()) {
                cout << "\tERROR detaching nodes " << leftID << " and " << rightID << endl;
                exit(EXIT_SUCCESS);
        }
#endif

        getSSNode(leftID).deleteRightArc(rightID);
        getSSNode(rightID).deleteLeftArc(leftID);
}

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
        for (auto it : pathA)
                removeNode(it);
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
                if (covFinalArcA < covFinalArcB)
                        detachNode(pathA[lastA], pathA[lastA+1]);
                else
                        detachNode(pathB[lastB], pathB[lastB+1]);
                return true;
        }

        // Remove first arc?
        double covFirstArcA = ((pathA.size() >= 2) && (firstA == 1)) ?
                getSSNode(pathA[0]).getRightArc(pathA[1])->getCoverage() : covCutoff + 1;
        double covFirstArcB = ((pathB.size() >= 2) && (firstB == 1)) ?
                getSSNode(pathB[0]).getRightArc(pathB[1])->getCoverage() : covCutoff + 1;
        if (min(covFirstArcA, covFirstArcB) <= covCutoff) {
                if (covFirstArcA < covFirstArcB)
                        detachNode(pathA[0], pathA[1]);
                else
                        detachNode(pathB[0], pathB[1]);
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
        priority_queue<PathDFS, vector<PathDFS>, PathDFSComp> heap;
        heap.push(PathDFS(srcID, 0));

        bool returnValue = false;

        while(!heap.empty()) {
                PathDFS currTop = heap.top();
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
                                PathDFS nextTop(nextID, nextLength);
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

bool DBGraph::flowCorrection(NodeID nodeID, double covCutoff)
{
        SSNode node = getSSNode(nodeID);
        if (!node.isValid())
                return false;

        int expNodeMult = getExpMult(node.getAvgKmerCov());
        if (expNodeMult == 0)
                return false;

        cout << "Multiplicity for node " << nodeID << ": " << expNodeMult << endl;

        int sumArcMult = 0;
        bool candidateRemoval = false;
        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                int expArcMult = getExpMult(it->getCoverage());
                if (expArcMult == 0) {
                        candidateRemoval = true;
                        expArcMult++;
                }
                sumArcMult += expArcMult;
        }

        cout << "Sum of the right arc multiplicities: " << sumArcMult << endl;

        // we will not detach arcs in this step
        if (sumArcMult <= expNodeMult && !candidateRemoval)
                return false;

        cout << "Sum of arcs is higher than expected multiplicity" << endl;

        // a) First assume that the topology is CORRECT
        double totCorrProb = getObsProb(node.getAvgKmerCov(), node.getMarginalLength(), sumArcMult);
        cout << "Log prob of node " << nodeID << " with coverage: " << node.getAvgKmerCov() << "  having multiplicity: " << sumArcMult << ": " << totCorrProb << endl;

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                int expArcMult = getExpMult(it->getCoverage());
                if (expArcMult == 0)    // bring this to one as we assume topology to be correct
                        expArcMult++;
                double arcProb = getObsProb(it->getCoverage(), 1, expArcMult);
                cout << "Log prob of arc with coverage " << it->getCoverage() << " having multiplicity: " << expArcMult << ": " << arcProb << endl;
                totCorrProb += arcProb;
        }

        cout << "TOTAL log prob assuming topology is correct: " << totCorrProb << endl;

        // b) Now assume that the topology is INCORRECT
        double totWrongProb = getObsProb(node.getAvgKmerCov(), node.getMarginalLength(), expNodeMult);
        cout << "Log prob of node " << nodeID << " with coverage: " << node.getAvgKmerCov() << "  having multiplicity: " << expNodeMult << ": " << totWrongProb << endl;

        vector<NodeID> toDetach;
        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                int expArcMult = getExpMult(it->getCoverage());
                if (expArcMult == 0)
                        toDetach.push_back(it->getNodeID());
                double arcProb = getObsProb(it->getCoverage(), 1, expArcMult);
                cout << "Log prob of arc with coverage " << it->getCoverage() << " having multiplicity: " << expArcMult << ": " << arcProb << endl;
                totWrongProb += arcProb;
        }

        cout << "TOTAL log prob assuming topology is WRONG: " << totWrongProb << endl;

        if (totWrongProb - totCorrProb < 5.0)
                return false;

        for (auto it : toDetach) {
                detachNode(nodeID, it);
                if ((getSSNode(it).getNumLeftArcs() == 0) &&
                    (getSSNode(it).getNumRightArcs() == 0))
                        getSSNode(it).invalidate();
        }

        return !toDetach.empty();
}

// ============================================================================
// PUBLIC CORRECTGRAPH.CPP (STAGE 4 ROUTINES)
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
        /*cout << "\t===== DEBUG: tip clipping report =====" << endl;
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
        cout << "\t===== DEBUG: end =====" << endl;*/
#endif

        return  numDeleted > 0;
}

bool DBGraph::concatenateNodes()
{
        size_t numConcatenations = 0, numIncorrectConcatenations = 0;

        for (NodeID seedID = 1; seedID <= numNodes; seedID++) {

                SSNode seed = getSSNode(seedID);
                if (!seed.isValid())
                        continue;

                deque<NodeID> nodeList;
                nodeList.push_back(seedID);
                seed.setFlag(true);

                // find linear paths to the right
                SSNode curr = seed;
                while (curr.getNumRightArcs() == 1) {
                        NodeID rightID = curr.rightBegin()->getNodeID();
                        SSNode right = getSSNode(rightID);
                        // don't merge palindromic repeats / loops
                        if (right.getFlag())
                                break;
                        if (right.getNumLeftArcs() != 1)
                                break;
                        nodeList.push_back(rightID);
                        right.setFlag(true);
                        curr = right;
                }

                // find linear paths to the left
                curr = seed;
                while (curr.getNumLeftArcs() == 1) {
                        NodeID leftID = curr.leftBegin()->getNodeID();
                        SSNode left = getSSNode(leftID);
                        // don't merge palindromic repeats / loops
                        if (left.getFlag())
                                break;
                        if (left.getNumRightArcs() != 1)
                                break;
                        nodeList.push_front(leftID);
                        left.setFlag(true);
                        curr = left;
                }

                // reset the flags to false
                for (const auto& it : nodeList)
                        getSSNode(it).setFlag(false);

                // if no linear path was found, continue
                if (nodeList.size() == 1)
                        continue;

                // concatenate the path
                NodeID frontID = nodeList.front();
                SSNode front = getSSNode(frontID);
                NodeID backID = nodeList.back();
                SSNode back = getSSNode(backID);

                front.deleteAllRightArcs();
                front.inheritRightArcs(back);

                size_t newKmerCov = front.getKmerCov();
                size_t newReadStartCov = front.getReadStartCov();
                for (size_t i = 1; i < nodeList.size(); i++) {
                        newKmerCov += getSSNode(nodeList[i]).getKmerCov();
                        newReadStartCov += getSSNode(nodeList[i]).getReadStartCov();
                        getSSNode(nodeList[i]).deleteAllLeftArcs();
                        getSSNode(nodeList[i]).deleteAllRightArcs();
                        getSSNode(nodeList[i]).invalidate();
                }

                front.setKmerCov(newKmerCov);
                front.setReadStartCov(newReadStartCov);

                vector<NodeID> nodeListv;
                copy(nodeList.begin(), nodeList.end(), std::back_inserter(nodeListv));

                string str;
                convertNodesToString(nodeListv, str);

                front.setSequence(str);

                numConcatenations += nodeList.size() - 1;

#ifdef DEBUG
                if (!trueMult.empty()) {
                        for (size_t i = 1; i < nodeListv.size(); i++) {
                                NodeID lID = nodeListv[i-1];
                                NodeID rID = nodeListv[i];
                                SSNode left = getSSNode(lID);
                                SSNode right = getSSNode(rID);
                                size_t lMult = trueMult[abs(lID)];
                                size_t rMult = trueMult[abs(rID)];
                                if (lMult != rMult) {
                                        cout << "\tConcatenating " << lID << " (cov = "
                                             << left.getAvgKmerCov() << ", " << lMult
                                             << ") and " << rID << " (cov = "
                                             << right.getAvgKmerCov() << ", " << rMult
                                             << ")" << endl;
                                        numIncorrectConcatenations++;
                                        trueMult[abs(lID)] = max(lMult, rMult);
                                        trueMult[abs(rID)] = max(lMult, rMult);
                                }
                        }
                }
#endif
        }

        cout << "Concatenated " << numConcatenations << " nodes" << endl;

#ifdef DEBUG
        if (numIncorrectConcatenations > 0)
                cout << "\t" << "Number of incorrect connections: "
                     << numIncorrectConcatenations << endl;
#endif

        return (numConcatenations > 0);
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

bool DBGraph::flowCorrection()
{
        cout << endl << "=================== Flow correction ===================" << endl;

        bool returnValue = false;
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                if (!getSSNode(id).isValid())
                        continue;
                if (getSSNode(id).getNumRightArcs() < 2)
                        continue;

                if (flowCorrection(id, getCovCutoff()))
                        returnValue = true;
        }

        return returnValue;
}
