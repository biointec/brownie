#include "alignment.h"
#include "graph.h"
#include "kmernode.h"
#include "readfile/fastafile.h"
#include "settings.h"
#include <gsl/gsl_randist.h>
#include <fstream>
#include <list>

#include <iomanip>
#include <string.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <vector>
#include "essaMEM-master/sparseSA.hpp"
#include <sstream>
#include <stack>
#include <algorithm>
#include <ctime>
#include <cstdio>
#include "library.h"


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

using namespace std;
class PathInfo {
public:
        NodeID nodeID;
        size_t length;

        PathInfo(NodeID nodeID, size_t length) : nodeID(nodeID), length(length) {};
};
struct pathStruct {
        SSNode rootNode;
        SSNode lastNode;
        SSNode firstNode;
        int pathLenght;
        int numOfNode;
        queue<SSNode> pathQueueNodes;
        pathStruct(SSNode root) {
                this->rootNode=root;
                this->lastNode=root;
                pathQueueNodes.push(root);
                //to avoid devision by zero, it should be zero by the way
                this->pathLenght=1;
                numOfNode=1;
        }
        pathStruct() {}
        void addNodeToPaht(SSNode node) {
                pathQueueNodes.push(node);
                pathLenght=pathLenght+node.getMarginalLength();
                this->lastNode=node;
                numOfNode++;
                if(numOfNode==2) {
                        firstNode=node;
                }

        }
        SSNode getLastNode() {
                return lastNode;
        }
        string getSequence() {
                queue<SSNode> tempQueue=this->pathQueueNodes;
                tempQueue.pop();
                string str;
                while(!tempQueue.empty()) {
                        SSNode node=tempQueue.front();
                        str=str+node.getSequence();
                        tempQueue.pop();
                }
                return str;
        }
        double getSingleNodePahtCoverage() {

                queue<SSNode> tempQueue=this->pathQueueNodes;
                double nodeCov=0;
                int i=0;
                //pop root
                tempQueue.pop();
                while(!tempQueue.empty()) {
                        SSNode node=tempQueue.front();
                        if(node.getNumLeftArcs()==1&&node.getNumRightArcs()==1) {
                                tempQueue.pop();
                                nodeCov=nodeCov+(node.getNodeKmerCov());
                                i++;
                        } else {
                                break;
                        }
                }
                if(i!=0) {
                        nodeCov=(nodeCov/i);
                        return nodeCov;
                }
                return-1;

        }
        double getSingleArcPathCoverage() {
                queue<SSNode> tempQueue=this->pathQueueNodes;
                double arcCov=0;
                if(tempQueue.size()>=2) {
                        SSNode firstNode=tempQueue.front();
                        tempQueue.pop();
                        int i=0;
                        while(!tempQueue.empty()) {
                                SSNode secondNode=tempQueue.front();
                                if(secondNode.getNumRightArcs()==1&& secondNode.getNumLeftArcs()==1) {
                                        tempQueue.pop();
                                        arcCov=arcCov+firstNode.getRightArc(secondNode.getNodeID())->getCoverage();
                                        firstNode=secondNode;
                                        i++;
                                }
                                else {
                                        break;
                                }

                        }
                        //the initial value is -1, so it should be added again
                        if(i!=0) {
                                arcCov=(arcCov)/(i);
                                return arcCov;

                        }

                }
                return -1;
        }
        int removeSingleNodes() {
                int remove=0;
                if(pathQueueNodes.size()>=3) {
                        SSNode rNode=pathQueueNodes.front();
                        pathQueueNodes.pop();
                        SSNode rrNode =pathQueueNodes.front();
                        pathQueueNodes.pop();
                        SSNode rrrNode=pathQueueNodes.front();
                        pathQueueNodes.pop();
                        if(rrNode.getNumLeftArcs()==1 && rrNode.getNumRightArcs()==1 ) {
                                bool result= rrrNode.deleteLeftArc(rrNode.getNodeID());
                                assert(result);
                                rrNode.deleteRightArc(rrrNode.getNodeID());
                                rrNode.deleteLeftArc(rNode.getNodeID());
                                rNode.deleteRightArc(rrNode.getNodeID());
                                rrNode.deleteAllLeftArcs();
                                rrNode.deleteAllRightArcs();
                                rrNode.invalidate();
                                remove++;

                        }
                        while(rrrNode.getNumLeftArcs()==0&& rrrNode.getNumRightArcs()==1&& !pathQueueNodes.empty()) {
                                rrNode=rrrNode;
                                rrrNode=pathQueueNodes.front();
                                pathQueueNodes.pop();
                                bool result= rrrNode.deleteLeftArc(rrNode.getNodeID());
                                assert(result);
                                rrNode.deleteRightArc(rrrNode.getNodeID());
                                rrNode.deleteAllLeftArcs();
                                rrNode.deleteAllRightArcs();
                                rrNode.invalidate();
                                remove++;
                        }
                }

                return remove;

        }


        bool operator<( const pathStruct & d ) const {
                return pathLenght <= d.pathLenght;
        }

};
struct comparator {
        bool operator()(const pathStruct& first, const pathStruct& second)
        {

                return first.pathLenght>=second.pathLenght;
        }
         bool operator()(const PathInfo& f, const PathInfo& s)
        {
                return f.length >= s.length;
        }
};

/*bool DBGraph::bubbleDetection(int round) {
        cout<<"*********************<<Bubble Detection starts>>......................................... "<<endl;
        priority_queue<pathStruct,vector<pathStruct>,comparator > MinHeap;
        int maxLength=maxNodeSizeToDel;//settings.getK()*2;
        size_t numOfDel=0;
        bool remove=false;
        size_t TP=0,TN=0,FP=0,FN=0;
        for ( NodeID lID = -numNodes; lID <= numNodes; lID++ ) {
                if ( lID == 0 )
                        continue;
                SSNode leftNode = getSSNode ( lID );
                if ( !leftNode.isValid() )
                        continue;

                if(leftNode.getNumRightArcs()<2)
                        continue;
                pathStruct rootPath(leftNode);
                MinHeap.push(rootPath);
                std::set<NodeID> visitedNodes;
                std::set<Arc *>visitedArc;
                visitedNodes.insert(leftNode.getNodeID());
                map<NodeID,pathStruct> pathDic;
                pathDic[rootPath.getLastNode().getNodeID()]=rootPath;
                while(!MinHeap.empty()) {
                        pathStruct leftPath =MinHeap.top();
                        MinHeap.pop();
                        leftNode=leftPath.getLastNode();
                        for ( ArcIt it = leftNode.rightBegin(); it != leftNode.rightEnd(); it++ ) {
                                //if(it->isValid()) {
                                SSNode rightNode=getSSNode(it->getNodeID());
                                Arc* p= leftNode.getRightArc(rightNode.getNodeID());
                                if(!(visitedArc.find(p)!=visitedArc.end())) {
                                        visitedArc.insert(leftNode.getRightArc(rightNode.getNodeID()));
                                        pathStruct extendedPath;
                                        extendedPath=leftPath;
                                        //+rightNode.getMarginalLength()
                                        if(extendedPath.pathLenght>=maxLength || rightNode.getNumRightArcs()==0)
                                                continue;
                                        extendedPath.addNodeToPaht(rightNode);
                                        MinHeap.push(extendedPath);
                                        if (visitedNodes.find(rightNode.getNodeID()) != visitedNodes.end()) {
                                                pathStruct prevPath=pathDic.at(rightNode.getNodeID());
                                                if (prevPath.firstNode.getNodeID()!=extendedPath.firstNode.getNodeID()) {
                                                        double lengthpro=0;
                                                        if(prevPath.pathLenght!=0)
                                                                lengthpro=extendedPath.pathLenght/prevPath.pathLenght;
                                                        if (lengthpro<.8 || lengthpro>1.2)
                                                                continue;
                                                        if (removeBubble(prevPath.firstNode ,extendedPath.firstNode,TP,TN,FP,FN, numOfDel))
                                                                break;
                                                }
                                        } else {
                                                visitedNodes.insert(rightNode.getNodeID());
                                                pathDic[rightNode.getNodeID()]=extendedPath;
                                        }
                                }
                        }
                }

        }
        #ifdef DEBUG
        cout<< "TP:     "<<TP<<"        TN:     "<<TN<<"        FP:     "<<FP<<"        FN:     "<<FN<<endl;
        cout << "Sensitivity: ("<<100*((double)TP/(double)(TP+FN))<<"%)"<<endl;
        cout<<"Specificity: ("<<100*((double)TN/(double)(TN+FP))<<"%)"<<endl;
        #endif
        cout<<"number of  deleted nodes based on bubble detection:            "<<numOfDel<<endl;
        if (numOfDel !=0)
                return true;
        return false;
}*/
bool DBGraph::bubbleDetection() {
        cout<<"*********************<<Bubble Detection starts>>......................................... "<<endl;
        size_t numOfDel=0;
        size_t TP=0,TN=0,FP=0,FN=0;

        std::set<NodeID> visitedNodes;
        std::set<Arc *>visitedArc;
        for ( NodeID lID = -numNodes; lID <= numNodes; lID++ ) {
                if ( lID == 0 )
                        continue;
                SSNode node = getSSNode ( lID );
                if ( !node.isValid() )
                        continue;
                if(node.getNumRightArcs()<2)
                        continue;
                SSNode firstNodeInUpPath;
                SSNode firstNodeInDoPath;
                if (lID++ % OUTPUT_FREQUENCY == 0)
                        (cout << "Extracting node -" <<numNodes<< "/ "<<lID<<" /"<<numNodes
                        << " from graph.\r").flush();
                vector<pair<SSNode, SSNode> > parallelPath =ExtractBubbles(node,visitedNodes ,visitedArc );
                size_t i=0;
                while (i<parallelPath.size()){
                        pair<SSNode, SSNode>  p=parallelPath[i];
                        if (p.first.isValid() && p.second.isValid())
                                removeBubble(p.first,p.second,TP,TN,FP,FN, numOfDel);
                        i++;
                }

        }
        #ifdef DEBUG
        cout<< "TP:     "<<TP<<"        TN:     "<<TN<<"        FP:     "<<FP<<"        FN:     "<<FN<<endl;
        cout << "Sensitivity: ("<<100*((double)TP/(double)(TP+FN))<<"%)"<<endl;
        cout<<"Specificity: ("<<100*((double)TN/(double)(TN+FP))<<"%)"<<endl;
        #endif
        cout<<"number of  deleted nodes based on bubble detection:            "<<numOfDel<<endl;
        if (numOfDel !=0)
                return true;
        return false;
}

vector<pair<SSNode, SSNode> >  DBGraph:: ExtractBubbles( SSNode rootNode,std::set<NodeID>& visitedNodes , std::set<Arc *>&visitedArc){
        priority_queue<pathStruct,vector<pathStruct>,comparator > MinHeap;
        int maxLength=maxNodeSizeToDel;//settings.getK()*2;
        pathStruct rootPath(rootNode);
        pathStruct extendedPath,prevPath;
        MinHeap.push(rootPath);
        visitedNodes.insert(rootNode.getNodeID());
        map<NodeID,pathStruct> pathDic;
        pathDic[rootPath.getLastNode().getNodeID()]=rootPath;
        vector<pair<SSNode, SSNode> > parallelNodes;
        while(!MinHeap.empty()) {
                pathStruct leftPath =MinHeap.top();
                MinHeap.pop();
                rootNode=leftPath.getLastNode();
                for ( ArcIt it = rootNode.rightBegin(); it != rootNode.rightEnd(); it++ ) {
                        //if(it->isValid()) {
                        SSNode rightNode=getSSNode(it->getNodeID());
                        Arc* p= rootNode.getRightArc(rightNode.getNodeID());
                        if(!(visitedArc.find(p)!=visitedArc.end())) {
                                visitedArc.insert(rootNode.getRightArc(rightNode.getNodeID()));
                                extendedPath=leftPath;
                                //+rightNode.getMarginalLength()
                                if(extendedPath.pathLenght>=maxLength || rightNode.getNumRightArcs()==0)
                                        continue;
                                extendedPath.addNodeToPaht(rightNode);
                                MinHeap.push(extendedPath);
                                if (visitedNodes.find(rightNode.getNodeID()) != visitedNodes.end()) {
                                        prevPath=pathDic.at(rightNode.getNodeID());
                                        if (prevPath.firstNode.getNodeID()!=extendedPath.firstNode.getNodeID()) {
                                                double lengthpro=0;
                                                if(prevPath.pathLenght!=0)
                                                        lengthpro=extendedPath.pathLenght/prevPath.pathLenght;
                                                if (lengthpro<.8 || lengthpro>1.2)
                                                        continue;
                                                parallelNodes.push_back(make_pair(prevPath.firstNode,extendedPath.firstNode));
                                        }
                                } else {
                                        visitedNodes.insert(rightNode.getNodeID());
                                        pathDic[rightNode.getNodeID()]=extendedPath;
                                }
                        }
                }
        }
        visitedArc.clear();
        visitedNodes.clear();
        return parallelNodes;

}

bool DBGraph::removeBubble(SSNode &prevFirstNode ,SSNode& extendFirstNode,size_t &TP,size_t &TN,size_t &FP,size_t &FN,size_t & numOfDel ){
        bool preIsSingle=true;
        bool exteIsSingle=true;
        if (prevFirstNode.getNumLeftArcs()>1 ||prevFirstNode.getNumRightArcs()>1)
                preIsSingle=false;
        if(extendFirstNode.getNumLeftArcs()>1||extendFirstNode.getNumRightArcs()>1)
                exteIsSingle=false;
        double preCov= prevFirstNode.getNodeKmerCov();
        double extCov=extendFirstNode.getNodeKmerCov();
        if(preIsSingle && exteIsSingle) {
                bool removePre=preCov<=extCov && preCov<this->cutOffvalue ?true:false;
                bool removeExt=extCov<preCov  && extCov<this->cutOffvalue ?true:false;
                if (!removeExt&&!removePre){
                        #ifdef DEBUG
                        if (trueMult[abs( prevFirstNode.getNodeID())]>0)
                                TN++;
                        else
                                FN++;
                        if (trueMult[abs(extendFirstNode.getNodeID())]>0)
                                TN++;
                        else
                                FN++;
                        #endif
                        return false;
                }
                if (removePre){
                         #ifdef DEBUG
                        if (trueMult[abs( prevFirstNode.getNodeID())]>0)
                                FP++;
                        else
                                TP++;
                        #endif
                        if(removeNode(prevFirstNode)){
                                numOfDel++;
                                return true;
                        }
                }
                if (removeExt){
                        #ifdef DEBUG
                        if(trueMult[abs( extendFirstNode.getNodeID())]>0)
                                FP++;
                        else
                                TP++;
                        #endif
                        if (removeNode(extendFirstNode)){
                                numOfDel++;
                                return true;

                        }
                }


        }
        if(preIsSingle && !exteIsSingle) {
                if (preCov<this->cutOffvalue && preCov<=extCov) {
                        #ifdef DEBUG
                        if (trueMult[abs( prevFirstNode.getNodeID())]>0)
                                FP++;
                        else
                                TP++;
                        #endif
                        if (removeNode(prevFirstNode)){
                                numOfDel++;
                                return true;

                        }
                }
                else{
                        #ifdef DEBUG
                        if (trueMult[abs( prevFirstNode.getNodeID())]>0)
                                TN++;
                        else
                                FN++;
                        #endif
                }
        }
        if(exteIsSingle&&!preIsSingle ) {
                if ((extCov<this->cutOffvalue && extCov<=preCov)) {
                        #ifdef DEBUG
                        if(trueMult[abs( extendFirstNode.getNodeID())]>0)
                                FP++;
                        else
                                TP++;
                        #endif
                        if(removeNode(extendFirstNode)){
                                numOfDel++;
                                return true;

                        }
                }else{
                        #ifdef DEBUG
                        if (trueMult[abs( extendFirstNode.getNodeID())]>0){
                                TN++;
                        }else
                                FN++;
                        #endif
                }
        }
        //when the starting nodes of both parallel paths has more than one ingoing or arcgoing arcs
        if (!exteIsSingle&&!preIsSingle)
                return removeNotSingleBubbles(prevFirstNode,extendFirstNode, TP,TN,FP,FN, numOfDel);
        return false;
}
bool DBGraph:: removeNotSingleBubbles( SSNode &prevFirstNode ,SSNode& extendFirstNode, size_t &TP,size_t &TN,size_t &FP,size_t &FN,size_t & numOfDel){
        double preCov=prevFirstNode.getNodeKmerCov();// prevFirstNode.getExpMult()/prevFirstNode.getMarginalLength();
        double extCov=extendFirstNode.getNodeKmerCov();//extendFirstNode.getExpMult()/extendFirstNode.getMarginalLength();
        if (preCov<extCov && preCov<this->cutOffvalue){
                #ifdef DEBUG
                if (trueMult[abs( prevFirstNode.getNodeID())]>0)
                        FP++;
                else
                        TP++;
                #endif
                if (removeNode(prevFirstNode)){
                        numOfDel++;
                        return true;}
        }
        else{
                #ifdef DEBUG
                if (trueMult[abs( prevFirstNode.getNodeID())]>0)
                        TN++;
                else
                        FN++;
                #endif
        }
        if (extCov<preCov && extCov<this->cutOffvalue){
                #ifdef DEBUG
                if(trueMult[abs( extendFirstNode.getNodeID())]>0)
                        FP++;
                else
                        TP++;
                #endif
                if (removeNode(extendFirstNode)){
                        numOfDel++;
                        return true;}
        }else{
                #ifdef DEBUG
                if (trueMult[abs( extendFirstNode.getNodeID())]>0){
                        TN++;
                }else
                        FN++;
                #endif
        }
        return false;
}
bool DBGraph::whichOneIsbubble(SSNode rootNode, bool &upIsBubble, SSNode &prevFirstNode ,SSNode& extendFirstNode, bool onlySingle, double threshold){
        if (whichOneIsbubble(rootNode, upIsBubble,prevFirstNode,extendFirstNode,onlySingle)) {
                if (upIsBubble){
                        if (prevFirstNode.getNodeKmerCov()<threshold) {
                                return true;
                        } else {
                                return false;
                        }
                }else{
                        if (extendFirstNode.getNodeKmerCov()<threshold) {
                                return true;
                        } else {
                                return false;
                        }
                }
        }
        return false;
}
bool DBGraph::whichOneIsbubble(SSNode rootNode,bool &first, SSNode &prevFirstNode ,SSNode& extendFirstNode, bool onlySingle){
        double preCov=prevFirstNode.getNodeKmerCov();
        double extCov=extendFirstNode.getNodeKmerCov();
        bool preIsSingle=true;
        bool exteIsSingle=true;
        if (prevFirstNode.getNumLeftArcs()>1 ||prevFirstNode.getNumRightArcs()>1)
                preIsSingle=false;
        if(extendFirstNode.getNumLeftArcs()>1||extendFirstNode.getNumRightArcs()>1)
                exteIsSingle=false;
        if(onlySingle){
                if(preIsSingle && exteIsSingle) {
                        bool removePre=preCov<=extCov ?true:false;//&&rootNode.getNodeKmerCov()/prevFirstNode.getNodeKmerCov()>3
                        bool removeExt=extCov<preCov  ?true:false;//&&rootNode.getNodeKmerCov()/extendFirstNode.getNodeKmerCov()>3
                        if (!removeExt&&!removePre)
                                return false;
                        if (removePre){
                                first=true;
                                return true;
                        }
                        if (removeExt){
                                first=false;
                                return true;
                        }
                }
                if(preIsSingle && !exteIsSingle ) {
                        if ( preCov<=extCov) {
                                first=true;
                                return true;
                        }
                }
                if(exteIsSingle&&!preIsSingle ) {
                        if ( extCov<=preCov) {
                                first=false;
                                return true;
                        }
                }
        }else{

                if (preCov<=extCov){//&&rootNode.getNodeKmerCov()/prevFirstNode.getNodeKmerCov()>4
                        first=true;
                        return true;
                }
                if (extCov<preCov ){//&&rootNode.getNodeKmerCov()/extendFirstNode.getNodeKmerCov()>4
                        first=false;
                        return true;
                }
        }

        return false;
}

void DBGraph::extractPath(NodeID currID, const vector<NodeID>& prevNode) const
{
        vector<NodeID> path;
        path.push_back(currID);

        while (true) {
                currID = prevNode[currID + numNodes];
                if (currID == 0)
                        break;
                path.push_back(currID);
        }

        reverse(path.begin(), path.end());
        for (auto it : path)
                cout << it << " ";
}
vector<NodeID> DBGraph::getPath(NodeID currID, const vector<NodeID>& prevNode) const
{
        vector<NodeID> path;
        path.push_back(currID);

        while (true) {
                currID = prevNode[currID + numNodes];
                if (currID == 0)
                        break;
                path.push_back(currID);
        }

        reverse(path.begin(), path.end());
        return path;
}

vector<pair<vector<NodeID>, vector<NodeID> > >  DBGraph::searchForParallelNodes(SSNode node,vector<NodeID> &visited, vector<NodeID> &prevNode,vector<NodeID> &nodeColor, int depth){
        size_t maxLength = depth;
        NodeID lID=node.getNodeID();
        priority_queue<PathInfo, vector<PathInfo>, comparator> heap;
        heap.push(PathInfo(lID, 0));
        vector<pair<vector<NodeID>, vector<NodeID> > > parallelNodes;
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
                        if ((prevNode[nextID + numNodes] != 0) || (nextID == lID)) {
                                if (nodeColor[nextID + numNodes] == nodeColor[currID + numNodes])
                                        continue;

                                NodeID upNodeID=nodeColor[currID + numNodes];
                                NodeID downNodeID=nodeColor[nextID + numNodes];
                                if (upNodeID!=0&&downNodeID!=0){
                                        vector<NodeID> upPathElements=getPath(currID, prevNode);
                                        vector<NodeID> downPathElements=getPath(nextID, prevNode);
                                        parallelNodes.push_back(make_pair(upPathElements,downPathElements));
                                }

                        } else {
                                prevNode[nextID + numNodes] = currID;
                                nodeColor[nextID + numNodes] = (currID == lID) ? nextID : nodeColor[currID + numNodes];
                                visited.push_back(nextID);

                                size_t nextLength = currLength + next.getMarginalLength();
                                if (nextLength > maxLength)
                                        continue;
                                PathInfo nextTop(nextID, nextLength);
                                heap.push(nextTop);
                        }
                }
        }

       for (auto it : visited) {
                prevNode[it + numNodes] = 0;
                nodeColor[it + numNodes] = 0;
        }
        visited.clear();

        return parallelNodes;
}
vector<pair<vector<NodeID>, vector<NodeID>> >  DBGraph::searchForParallelNodes(SSNode node, int depth){
        vector<NodeID> prevNode(2*numNodes+1, 0);
        vector<NodeID> nodeColor(2*numNodes+1, 0);
        vector<NodeID> visited;
        return (searchForParallelNodes(node, visited,nodeColor,prevNode, depth));
}
bool DBGraph::bubbleDetection(int depth) {

        size_t numOfDel=0;
        size_t TP=0,TN=0,FP=0,FN=0;
        vector<NodeID> prevNode(2*numNodes+1, 0);
        vector<NodeID> nodeColor(2*numNodes+1, 0);
        vector<pair<NodeID, NodeID> > parallelNodes;
        vector<NodeID> visited;
        for (NodeID lID = -numNodes; lID <numNodes; lID++) {
                if (lID % OUTPUT_FREQUENCY == 0)
                        (cout << "Extracting node -" <<numNodes<< "/ "<<lID<<" /"<<numNodes
                        << " from graph.\r").flush();
                if ( lID == 0 )
                        continue;
                SSNode node = getSSNode(lID);
                if (!node.isValid())
                        continue;
                // only consider nodes that branch
                if (node.getNumRightArcs() < 2)
                        continue;
                if (!hasLowCovNode(node))
                        continue;
                vector<pair<vector<NodeID>, vector<NodeID>> > parallelPath=searchForParallelNodes(node, visited,nodeColor,prevNode, depth);
                for (auto it : parallelPath){

                        vector<NodeID> upPath=it.first;
                        vector<NodeID> downPath=it.second;

                        SSNode upLast=getSSNode(it.first[upPath.size()-1]);
                        SSNode downLast=getSSNode(it.second[downPath.size()-1]);
                        SSNode up=getSSNode(it.first[1]);
                        SSNode down=getSSNode(it.second[1]);
                        if(up.isValid()&&down.isValid())
                        {
                                bool upIsBubble=true;
                                bool bubbleDeleted=false;
                                if (whichOneIsbubble(node,upIsBubble,up, down,false,this->cutOffvalue)){
                                        if (upIsBubble){
                                                #ifdef DEBUG
                                                size_t mul=trueMult[abs( up.getNodeID())];
                                                #endif
                                                if (removeNode(up)){
                                                        #ifdef DEBUG
                                                        if (mul>0)
                                                                FP++;
                                                        else
                                                                TP++;
                                                        #endif
                                                        bubbleDeleted=true;
                                                        numOfDel++;
                                                }
                                                if (upLast.isValid()&&upLast.getNodeKmerCov()<cutOffvalue){//&& node.getNodeKmerCov()/upLast.getNodeKmerCov()>3){
                                                        #ifdef DEBUG
                                                        mul=trueMult[abs( upLast.getNodeID())];
                                                        #endif
                                                        if (removeNode(upLast)){
                                                                #ifdef DEBUG
                                                                if (mul>0)
                                                                        FP++;
                                                                else
                                                                        TP++;
                                                                #endif
                                                                bubbleDeleted=true;
                                                                numOfDel++;
                                                        }
                                                }
                                        }
                                        else{
                                                #ifdef DEBUG
                                                size_t mul=trueMult[abs( down.getNodeID())];
                                                #endif
                                                if( removeNode(down)){
                                                        #ifdef DEBUG
                                                        if (mul>0)
                                                                FP++;
                                                        else
                                                                TP++;
                                                        #endif
                                                        bubbleDeleted=true;
                                                        numOfDel++;
                                                }
                                                if (downLast.isValid()&&downLast.getNodeKmerCov()<cutOffvalue){// && node.getNodeKmerCov()/downLast.getNodeKmerCov()>3){
                                                        #ifdef DEBUG
                                                        mul=trueMult[abs( downLast.getNodeID())];
                                                        #endif
                                                        if( removeNode(downLast)){
                                                                #ifdef DEBUG
                                                                if (mul>0)
                                                                        FP++;
                                                                else
                                                                        TP++;
                                                                #endif
                                                                bubbleDeleted=true;
                                                                numOfDel++;
                                                        }
                                                }


                                        }
                                }
                                if(!bubbleDeleted){
                                        #ifdef DEBUG
                                        if (trueMult[abs( up.getNodeID())]>0)
                                                TN++;
                                        else
                                                FN++;
                                        if (trueMult[abs( down.getNodeID())]>0)
                                                TN++;
                                        else
                                                FN++;
                                        #endif
                                }

                        }
                }
        }
        cout<<endl;
        #ifdef DEBUG
        cout<<endl<< "TP:     "<<TP<<"        TN:     "<<TN<<"        FP:     "<<FP<<"        FN:     "<<FN<<endl;
        cout << "Sensitivity: ("<<100*((double)TP/(double)(TP+FN))<<"%)"<<endl;
        cout<<"Specificity: ("<<100*((double)TN/(double)(TN+FP))<<"%)"<<endl;
        #endif
        if (numOfDel>0)
        cout<<"number of  deleted nodes based on bubble detection: "<<numOfDel<<endl;
        if (numOfDel !=0)
                return true;
        return false;
}
double DBGraph::getMeanPathCov( vector<NodeID> &path){
        SSNode last=getSSNode(path[path.size()-1]);
        SSNode first=getSSNode(path[1]);
        return ((last.getNodeKmerCov()+first.getNodeKmerCov())/2);
}
bool DBGraph::detectFalsePath( vector<NodeID> &upPath, vector<NodeID> &downPath){
        SSNode up=getSSNode(upPath[1]);
        SSNode down=getSSNode(downPath[1]);
        float upCov= up.getNodeKmerCov(); //getMeanPathCov(upPath);
        float downCov=down.getNodeKmerCov(); ///getMeanPathCov(downPath);
        if (downCov<=upCov && downCov<this->cutOffvalue)
                removePath(downPath,upPath);
        if(upCov<=downCov&&upCov<this->cutOffvalue)
                removePath(upPath,downPath);
        return true;
}
bool DBGraph::removePath(vector<NodeID> &delPath, vector<NodeID> &correctPath){
        SSNode lastToDel=getSSNode(delPath[delPath.size()-1]);
        SSNode firstToDel=getSSNode(delPath[1]);
        if (firstToDel.getNodeKmerCov()<cutOffvalue){
                removeNode(firstToDel);
        }
        if(lastToDel.getNodeID()!=firstToDel.getNodeID()){
                if (lastToDel.getNodeKmerCov()<cutOffvalue)
                removeNode(lastToDel);
        }
        return true;
}
bool DBGraph::hasLowCovNode(SSNode root){

        ArcIt it = root.rightBegin();
        for (int i=0;i<root.getNumRightArcs();i++){
                SSNode node =getSSNode(it->getNodeID());
                if (node.getNodeKmerCov()<this->cutOffvalue)
                        return true;
                it++;
        }
        return false;

}
/*
bool DBGraph::bubbleDetection(int depth) {
        cout << " ======== Bubble Detection ======== " << endl;


        #ifdef DEBUG
        size_t numOfDel=0;
        size_t TP=0,TN=0,FP=0,FN=0;
        size_t initialTrue=0, finalTrue=0;
        size_t initialFalse=0, finalFalse=0;

        for(NodeID lID = 1; lID <numNodes; lID++){
                SSNode node = getSSNode(lID);
                if (!node.isValid())
                        continue;
                if (trueMult[abs(node.getNodeID())]>0)
                        initialTrue++;
                else
                        initialFalse++;
        }

        #endif

        vector<NodeID> prevNode(2*numNodes+1, 0);
        vector<NodeID> nodeColor(2*numNodes+1, 0);
        vector<pair<NodeID, NodeID> > parallelNodes;
        vector<NodeID> visited;
        for (NodeID lID = -numNodes; lID <numNodes; lID++) {
                if (lID % OUTPUT_FREQUENCY == 0)
                        (cout << "Extracting node -" <<numNodes<< "/ "<<lID<<" /"<<numNodes
                        << " from graph.\r").flush();
                if ( lID == 0 )
                        continue;
                SSNode node = getSSNode(lID);
                if (!node.isValid())
                        continue;
                // only consider nodes that branch
                if (node.getNumRightArcs() < 2)
                        continue;
                if (!hasLowCovNode(node))
                        continue;
                vector<pair<vector<NodeID>, vector<NodeID>> > parallelPath=searchForParallelNodes(node, visited,nodeColor,prevNode, depth);
                for (auto it : parallelPath)
                        detectFalsePath(it.first,it.second);

        }
        #ifdef DEBUG

        for(NodeID lID = 1; lID <numNodes; lID++){
                SSNode node = getSSNode(lID);
                if (!node.isValid())
                        continue;
                if (trueMult[abs(node.getNodeID())]>0)
                        finalTrue++;
                else
                        finalFalse++;
        }
        TP=initialFalse-finalFalse;
        TN=finalTrue;
        FN=finalFalse;
        FP=initialTrue-finalTrue;
        cout<< "TP:     "<<TP<<"        TN:     "<<TN<<"        FP:     "<<FP<<"        FN:     "<<FN<<endl;
        cout << "Precision: ("<<100*((double)TP/(double)(TP+FP))<<"%)"<<endl;

        //cout<<"Specificity: ("<<100*((double)TN/(double)(TN+FP))<<"%)"<<endl;
        #endif
        numOfDel=(initialTrue+initialFalse)-(finalFalse+finalTrue);
        cout<<"number of  deleted nodes based on bubble detection:            "<<numOfDel<<endl;
        if ( numOfDel!=0)
                return true;
        return false;
}*/



