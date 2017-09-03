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
#include "alignment.h"
#include <queue>
#include <map>
#include "component.h"
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
        numValidNodes--;
}

void DBGraph::flagNode(NodeID nodeID)
{
        getSSNode(nodeID).setFlag(true);
}

void DBGraph::removeArc(NodeID leftID, NodeID rightID)
{
#ifdef DEBUG
        if (getSSNode(leftID).getRightArc(rightID)->getTrueArc()) {
                cout << "\tERROR detaching nodes " << leftID << " and " << rightID << endl;
                //exit(EXIT_SUCCESS);
        }
#endif

        getSSNode(leftID).deleteRightArc(rightID);
        getSSNode(rightID).deleteLeftArc(leftID);
}

void DBGraph::flagArc(NodeID leftID, NodeID rightID)
{
        getSSNode(leftID).getRightArc(rightID)->setFlag(true);
        getSSNode(rightID).getLeftArc(leftID)->setFlag(true);
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
                            size_t& first, size_t& last) const
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
string DBGraph::getPathSeq(const vector<NodeID>& path){
        string seq = "";
        for (int i=1;i<path.size()-1;i++) {
                seq = seq + getSSNode(path[i]).substr(Kmer::getK()-1, getSSNode(path[i]).getMarginalLength());
        }
        return seq;

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

void DBGraph::flagPath(const vector<NodeID>& pathA)
{
        for (NodeID nodeID : pathA)
                flagNode(nodeID);
}

bool DBGraph::flowCorrection(NodeID nodeID, double covCutoff, size_t maxMargLength)
{
        SSNode node = getSSNode(nodeID);
        if (!node.isValid())
                return false;

        int expNodeMult = getExpMult(node.getAvgKmerCov());
        if (expNodeMult == 0)
                return false;

        //cout << "Multiplicity for node " << nodeID << ": " << expNodeMult << endl;

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

        //cout << "Sum of the right arc multiplicities: " << sumArcMult << endl;

        // we will not detach arcs in this step
        if (sumArcMult <= expNodeMult && !candidateRemoval)
                return false;

        //cout << "Sum of arcs is higher than expected multiplicity" << endl;

        // a) First assume that the topology is CORRECT
        double totCorrProb = getObsProbLog(node.getAvgKmerCov(), node.getMarginalLength(), sumArcMult);
        //cout << "Log prob of node " << nodeID << " with coverage: " << node.getAvgKmerCov() << "  having multiplicity: " << sumArcMult << ": " << totCorrProb << endl;

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                int expArcMult = getExpMult(it->getCoverage());
                if (expArcMult == 0)    // bring this to one as we assume topology to be correct
                        expArcMult++;
                double arcProb = getObsProbLog(it->getCoverage(), 1, expArcMult);
                //cout << "Log prob of arc with coverage " << it->getCoverage() << " having multiplicity: " << expArcMult << ": " << arcProb << endl;
                totCorrProb += arcProb;
        }

        //cout << "TOTAL log prob assuming topology is correct: " << totCorrProb << endl;

        // b) Now assume that the topology is INCORRECT
        double totWrongProb = getObsProbLog(node.getAvgKmerCov(), node.getMarginalLength(), expNodeMult);
        //cout << "Log prob of node " << nodeID << " with coverage: " << node.getAvgKmerCov() << "  having multiplicity: " << expNodeMult << ": " << totWrongProb << endl;

        vector<NodeID> toDetach;
        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                int expArcMult = getExpMult(it->getCoverage());
                if (expArcMult == 0)
                        toDetach.push_back(it->getNodeID());
                double arcProb = getObsProbLog(it->getCoverage(), 1, expArcMult);
                //cout << "Log prob of arc with coverage " << it->getCoverage() << " having multiplicity: " << expArcMult << ": " << arcProb << endl;
                totWrongProb += arcProb;
        }

        //cout << "TOTAL log prob assuming topology is WRONG: " << totWrongProb << endl;

        if (totWrongProb - totCorrProb < 5)
                return false;
        bool changed = false;
        for (auto it : toDetach) {
                if (getSSNode(it).getAvgKmerCov() > covCutoff || getSSNode(it).getMarginalLength() > maxMargLength)
                        continue;
                removeArc(nodeID, it);
                if ((getSSNode(it).getNumLeftArcs() == 0) &&
                    (getSSNode(it).getNumRightArcs() == 0)) {
                        getSSNode(it).invalidate();
                        changed = true;
                        numValidNodes--;
                }
        }
        return changed ;
        //return !toDetach.empty();
}

// ============================================================================
// PUBLIC CORRECTGRAPH.CPP (STAGE 4 ROUTINES)
// ============================================================================

bool DBGraph::clipNormalTip(SSNode startNode ){
        SSNode alternative;
        SSNode currNode = getSSNode(startNode.rightBegin()->getNodeID());
        if (currNode.getNumLeftArcs()<=1)
                return false;
        ArcIt it = currNode.leftBegin();
        do{
                alternative = getSSNode(it->getNodeID());
                if (startNode.getNodeID()!=alternative.getNodeID() ){ // alternative.getAvgKmerCov()>=startNode.getAvgKmerCov()
                        string currStr="", altStr="";
                        currStr = startNode.getSequence();
                        altStr = alternative.getSequence();
                        // this is in favor or longer node to be saved.
                        //For the small ones the last k-1 base is exactly the same, more likely will be deleted
                        currStr = currStr.substr(currStr.length() - min (currStr.length(),altStr.length() ),min (currStr.length(),altStr.length() ));
                        altStr = altStr.substr(altStr.length()- min( currStr.length(),altStr.length() ),min (currStr.length(),altStr.length() ));

                        altStr = altStr.substr(0,altStr.length()-settings.getK()+1);
                        currStr = currStr.substr(0,currStr.length()-settings.getK()+1);
                        if ( currStr.length() <settings.getK()||
                                alignment.align(currStr,altStr)> ( (int) max( currStr.length(),altStr.length() ) / 3)){
                                return true;
                        }

                }
                it++;
        }
        while (it != currNode.leftEnd());
        return false;
}




bool DBGraph::clipJoinedTip(double covCutoff,size_t maxMargLength ,SSNode startNode ){


        bool remove = true;
        ArcIt it = startNode.rightBegin();
        while (it != startNode.rightEnd()){
                double arcCov = startNode.getRightArc(getSSNode(it->getNodeID()).getNodeID())->getCoverage();
                if (arcCov > covCutoff)
                        remove = false;
                it++;
        }
        if (startNode.getAvgKmerCov() > covCutoff || startNode.getMarginalLength() > maxMargLength )
                remove = false;
        if ( !remove){
                ArcIt it = startNode.rightBegin();
                SSNode nodeBefore;
                double arcCovMin = startNode.getRightArc(getSSNode(it->getNodeID()).getNodeID())->getCoverage();
                while (it != startNode.rightEnd()){
                        double arcCov = startNode.getRightArc(getSSNode(it->getNodeID()).getNodeID())->getCoverage();
                        if (arcCov <= arcCovMin){
                                nodeBefore = getSSNode(it->getNodeID());
                                arcCovMin = arcCov;
                        }
                        it++;
                }
                if (arcCovMin <= covCutoff){
                        removeArc(startNode.getNodeID(),nodeBefore.getNodeID());
                }
        }
        //we don't want to remove any join tip
        return remove;
}


bool DBGraph::clipTips(double covCutoff, size_t maxMargLength)
{

#ifdef DEBUG
        size_t tp=0, tn=0, fp=0,fn=0;
        size_t tps=0, tns=0, fps=0,fns=0;
        size_t tpj=0, tnj=0, fpj=0,fnj=0;
#endif

        size_t numDeleted = 0;
        //remove single loops
        for (NodeID id = 1; id <= numNodes; id++) {

                SSNode node = getSSNode(id);

                if (node.getAvgKmerCov() > (covCutoff) )
                        continue;
                if (node.getMarginalLength() > maxMargLength )
                        continue;

                if (node.getNumRightArcs() ==1 && node.getNumLeftArcs()==1&& node.getLeftArc(node.getNodeID()) != NULL){
                        removeNode(node.getNodeID());
                        numDeleted ++ ;
                }
        }
        //remove palindromic tips
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id ==0)
                        continue;
                SSNode first  = getSSNode( id);
                if (first.getNumRightArcs() != 1)
                        continue;
                if (first.getNumLeftArcs()  != 1)
                        continue;
                SSNode second  = getSSNode( first.rightBegin()->getNodeID());
                if (second.getRightArc(first.getNodeID()) ==NULL)
                        continue;
                if (first.getAvgKmerCov() > (covCutoff) )
                        continue;
                if (first.getMarginalLength() > maxMargLength )
                        continue;

                removeNode(first.getNodeID());
                numDeleted ++ ;
        }
        for (NodeID id = 1; id <= numNodes; id++) {

                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;
                // check for dead ends
                bool leftDE = (node.getNumLeftArcs() == 0);
                bool rightDE = (node.getNumRightArcs() == 0);
                if (!leftDE && !rightDE)
                        continue;
                SSNode startNode = (rightDE) ? getSSNode(-id) : getSSNode(id);
                bool isolated = false;
                bool joinedTip = false;
                bool remove = false;
                //isolated tips
                if (startNode.getNumRightArcs() == 0 ){
                        isolated = true;
                        remove = true;
                }

                if (startNode.getNumRightArcs() == 1){
                        remove = clipNormalTip( startNode );
                }

                if (startNode.getNumRightArcs()>1){
                        joinedTip = true;
                        remove = clipJoinedTip(covCutoff, maxMargLength, startNode);
                }
                if (startNode.getAvgKmerCov() > (covCutoff) )
                        continue;
                if (startNode.getMarginalLength() > maxMargLength && !isolated)
                        continue;
                if (startNode.getMarginalLength() > maxMargLength*2 && isolated)
                        continue;
                if (remove){
                        removeNode(startNode.getNodeID());
                        numDeleted ++;
                }
                #ifdef DEBUG
                if (remove) {
                        //cout << "node " <<id << " detected as a tip and removed" <<endl;
                        if (trueMult.size()>0&& trueMult[abs(id)] > 0) {
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
                        if (trueMult.size()>0&& trueMult[abs(id)] > 0) {
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
        cout << "\tClipped " << numDeleted << " nodes" << endl;

#ifdef DEBUG
        cout << "\t===== DEBUG: tip clipping report =====" << endl;
        cout << "\tIsolated TP: " << tps << "\tTN: "<< tns << "\tFP: " << fps << "\tFN: "<< fns << endl;
        cout << "\tSensitivity: " << 100.0 * Util::getSensitivity(tps, fns) << "%" << endl;
        cout << "\tSpecificity: " << 100.0 * Util::getSpecificity(tns, fps) << "%" << endl;
        cout << "\t****************************************" << endl;
        cout << "\t****************************************" << endl;
        cout << "\tNormal TP: " << tp << "\tTN: " << tn << "\tFP: " << fp << "\tFN: " << fn << endl;
        cout << "\tSensitivity: " << 100.0 * Util::getSensitivity(tp, fn) << "%" << endl;
        cout << "\tSpecificity: " << 100.0 * Util::getSpecificity(tn, fp) << "%" << endl;
        cout << "\t===== DEBUG: end =====" << endl;
#endif

        return  numDeleted > 0 ;
}

void DBGraph::concatenateAroundNode(NodeID seedID, vector<NodeID>& nodeListv)
{
        nodeListv.clear();

        SSNode seed = getSSNode(seedID);
        if (!seed.isValid())
                return;

        deque<NodeID> nodeListq;
        nodeListq.push_back(seedID);
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
                nodeListq.push_back(rightID);
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
                nodeListq.push_front(leftID);
                left.setFlag(true);
                curr = left;
        }

        // reset the flags to false
        for (const auto& it : nodeListq)
                getSSNode(it).setFlag(false);

        // if no linear path was found, continue
        if (nodeListq.size() == 1)
                return;

        // concatenate the path
        NodeID frontID = nodeListq.front();
        SSNode front = getSSNode(frontID);
        NodeID backID = nodeListq.back();
        SSNode back = getSSNode(backID);

        front.deleteAllRightArcs();
        front.inheritRightArcs(back);

        size_t newKmerCov = front.getKmerCov();
        size_t newReadStartCov = front.getReadStartCov();
        for (size_t i = 1; i < nodeListq.size(); i++) {
                newKmerCov += getSSNode(nodeListq[i]).getKmerCov();
                newReadStartCov += getSSNode(nodeListq[i]).getReadStartCov();
                getSSNode(nodeListq[i]).deleteAllLeftArcs();
                getSSNode(nodeListq[i]).deleteAllRightArcs();
                getSSNode(nodeListq[i]).invalidate();
        }

        front.setKmerCov(newKmerCov);
        front.setReadStartCov(newReadStartCov);

        copy(nodeListq.begin(), nodeListq.end(), std::back_inserter(nodeListv));

        string str;
        convertNodesToString(nodeListv, str);

        front.setSequence(str);
        numValidNodes -= nodeListq.size() - 1;
}

bool DBGraph::concatenateNodes()
{
        size_t numConcatenations = 0;
#ifdef DEBUG
        size_t numIncorrectConcatenations = 0;
#endif

        for (NodeID seedID = 1; seedID <= numNodes; seedID++) {
                vector<NodeID> concatenation;
                concatenateAroundNode(seedID, concatenation);

                if (!concatenation.empty())
                        numConcatenations += concatenation.size() - 1;

#ifdef DEBUG

                if (trueMult.empty())
                        continue;

                for (size_t i = 1; i < concatenation.size(); i++) {
                        NodeID lID = concatenation[i-1];
                        NodeID rID = concatenation[i];
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
#endif
        }

        cout << "\tConcatenated " << numConcatenations << " nodes" << endl;

#ifdef DEBUG
        if (numIncorrectConcatenations > 0)
                cout << "\t" << "Number of incorrect connections: "
                     << numIncorrectConcatenations << endl;
#endif

        size_t countTotal = 0;
        for (NodeID seedID = 1; seedID <= numNodes; seedID++)
                if (getSSNode(seedID).isValid())
                        countTotal++;

        return (numConcatenations > 0);
}


bool DBGraph::flowCorrection(double covCutoff, size_t maxMargLength)
{
        bool returnValue = false;
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                if (!getSSNode(id).isValid())
                        continue;
                if (getSSNode(id).getNumRightArcs() < 2)
                        continue;

                if (abs(id) % OUTPUT_FREQUENCY == 0) {
                        cout << "\tProcessing node " << id << "/" << numNodes << "\r";
                        cout.flush();
                }

                if (flowCorrection(id, covCutoff, maxMargLength))
                        returnValue = true;
        }

        cout << "\tProcessing node " << numNodes << "/" << numNodes << endl;

        return returnValue;
}


bool DBGraph::handleParallelPaths(const vector<NodeID>& pathA,
                                  const vector<NodeID>& pathB,
                                  double covCutoff, size_t maxMargLength)
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
       bool remove = true;
       if (!lowCovPath.empty()) {
               if (lowCov == covSubPathA){
                       if (covSubPathB < covSubPathA * 2)
                               remove = false;
               }else{
                       if (covSubPathA < covSubPathB * 2)
                               remove = false;
               }
       }
       AlignmentJan ali(250, 2, 1, -1, -3);

       string pathAstr = getPathSeq(pathA);
       string pathBstr = getPathSeq(pathB);
       if (abs( pathAstr.length()- pathBstr.length())> 0 )
                       remove = false;
       if ( pathAstr.length() >settings.getK() && ali.align(pathAstr,pathBstr)<((int)min( pathAstr.length(),pathBstr.length() ) / 3))
                       remove =false;
       if (remove && (lowCov <= covCutoff) && !lowCovPath.empty()) {
                //removePath(lowCovPath);
                flagPath(lowCovPath);
                return true;
        }

        // Remove final arc?
        double covFinalArcA = ((pathA.size() >= 2) && (lastA == pathA.size() - 2)) ?
                getSSNode(pathA[lastA]).getRightArc(pathA[lastA+1])->getCoverage() : covCutoff + 1;
        double covFinalArcB = ((pathB.size() >= 2) && (lastB == pathB.size() - 2)) ?
                getSSNode(pathB[lastB]).getRightArc(pathB[lastB+1])->getCoverage() : covCutoff + 1;
        if (min(covFinalArcA, covFinalArcB) <= covCutoff) {
                if (covFinalArcA < covFinalArcB)
                        //removeArc(pathA[lastA], pathA[lastA+1]);
                        if (covFinalArcA *2 < covFinalArcB)
                        flagArc(pathA[lastA], pathA[lastA+1]);
                else
                        //removeArc(pathB[lastB], pathB[lastB+1]);
                        if (covFinalArcB *2 <covFinalArcA )
                        flagArc(pathB[lastB], pathB[lastB+1]);
                return true;
        }

        // Remove first arc?
        double covFirstArcA = ((pathA.size() >= 2) && (firstA == 1)) ?
                getSSNode(pathA[0]).getRightArc(pathA[1])->getCoverage() : covCutoff + 1;
        double covFirstArcB = ((pathB.size() >= 2) && (firstB == 1)) ?
                getSSNode(pathB[0]).getRightArc(pathB[1])->getCoverage() : covCutoff + 1;
        if (min(covFirstArcA, covFirstArcB) <= covCutoff) {
                if (covFirstArcA < covFirstArcB)
                        //removeArc(pathA[0], pathA[1]);
                        if (covFirstArcA*2 < covFirstArcB)
                        flagArc(pathA[0], pathA[1]);
                else
                        //removeArc(pathB[0], pathB[1]);
                        if (covFirstArcB*2 <covFirstArcA )
                        flagArc(pathB[0], pathB[1]);
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
        if (getSSNode(srcID).getNumRightArcs() < 2)
                return false;
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
                                if (handleParallelPaths(pathA, pathB, covCutoff, maxMargLength)) {
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

void DBGraph::bubbleDetectionThread(size_t threadID, ParGraph& wlb,
                                    double covCutoff, size_t maxMargLength)
{
        vector<NodeID> visited;
        vector<NodeID> prevNode(2*numNodes+1, 0);
        vector<NodeID> nodeColor(2*numNodes+1, 0);

        while (true) {
                size_t firstNode, numNodes;
                wlb.getNodeChunk(firstNode, numNodes);

                if (numNodes == 0)
                        break;

                //cout << "Work from " << threadID << " " << firstNode << " to " << firstNode + numNodes << endl;
                for (size_t id = firstNode; id < firstNode + numNodes; id++) {
                        SSNode node = getSSNode(id);
                        if (!node.isValid())
                                continue;
                        // handle the positive node
                        bubbleDetection(id, visited, prevNode, nodeColor,
                                        covCutoff, maxMargLength,
                                        settings.getBubbleDFSNodeLimit());

                        // handle the negative node
                        bubbleDetection(-id, visited, prevNode, nodeColor,
                                        covCutoff, maxMargLength,
                                        settings.getBubbleDFSNodeLimit());
                }
        }
}

bool DBGraph::bubbleDetection(double covCutoff, size_t maxMargLength)
{
        const unsigned int& numThreads = settings.getNumThreads();
        ParGraph wlb(numNodes, settings.getThreadBubbleWorkSize());

        for (NodeID id = 1; id <= numNodes; id++)
                getSSNode(id).setFlag(false);

        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&DBGraph::bubbleDetectionThread, this,
                                          i, ref(wlb), covCutoff, maxMargLength);

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        cout << "\tProcessing graph (100%) " << endl;
        bool returnValue = false; size_t numNodesRemoved = 0;

        for (NodeID id = 1; id <= numNodes; id++) {
                if (getSSNode(id).getFlag()) {
                        removeNode(id);
                        //cout << "node " <<id << " detected as a bubble and removed" <<endl;

                        returnValue = true;
                        numNodesRemoved++;
                }
        }
        cout << "\tRemoved " << numNodesRemoved << " nodes" << endl;

        size_t numArcsRemoved = 0;
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;

                vector<pair<NodeID, NodeID> > toRemove;

                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        if (node.getRightArc(it->getNodeID())->getFlag())
                                toRemove.push_back(pair<NodeID, NodeID>(id, it->getNodeID()));
                }

                for (auto it : toRemove) {
                        removeArc(it.first, it.second);
                        returnValue = true;
                        numArcsRemoved++;
                }
        }
        cout << "\tRemoved " << numArcsRemoved << " arcs" << endl;

        return numNodesRemoved > 20;
}

size_t DBGraph::findbreakpoints(std::string breakpointFileName )
{
        std::ifstream input(breakpointFileName);
        vector< pair<string, string> > breakpoints;
        if (!input.good()) {
                std::cerr << "Error opening: " << breakpointFileName << std::endl;
                return -1;
        }
        std::string line, id ="", DNA_sequence ="" ;
        while (std::getline(input, line).good()) {
                if (line[0] == '>') {
                        if (DNA_sequence!= "")
                                breakpoints.push_back( make_pair(id, DNA_sequence));
                        id = line.substr(1);
                        DNA_sequence.clear();
                }
                else if (line[0] != '>')
                        DNA_sequence += line;
        }
        breakpoints.push_back( make_pair(id, DNA_sequence));

        size_t numberOfBreakpoints = 0 ;
        for (int i = 0 ;i <breakpoints.size(); i++){
                pair<string , string> breakPoint= breakpoints [i];
                string refSeq = breakPoint.second;
                //cout << " we are looking at ref ID : " << breakPoint.first <<endl;
                // handle the other kmers
                NodePosPair prev;

                for (KmerIt it(refSeq); it.isValid(); it++) {
                        Kmer kmer = it.getKmer();
                        NodePosPair curr = findNPP(kmer);
                        if (!curr.isValid())
                                continue;
                        if (!prev.isValid()){
                                prev = curr;
                                continue;
                        }
                        if (!consecutiveNPP(prev, curr) && prev.getNodeID() != curr.getNodeID())
                        {
                                cout << "new breakpoint happend in " <<  breakPoint.first <<endl;
                                cout << "the connection between node " << prev.getNodeID() << " and " <<curr.getNodeID()  << " is lost." <<endl;
                                numberOfBreakpoints ++;
                                break;

                        }
                        prev = curr;
                }
        }
        cout << "There are " << numberOfBreakpoints << " breakpoins in the input file. "<<endl;
        return numberOfBreakpoints;
}

