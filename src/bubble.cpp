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
#include "settings.h"
#include "library.h"

#include <queue>

using namespace std;
class PathInfo {
public:
        NodeID nodeID;
        size_t length;

        PathInfo(NodeID nodeID, size_t length) : nodeID(nodeID), length(length) {};
};

struct comparator {

         bool operator()(const PathInfo& f, const PathInfo& s)
        {
                return f.length >= s.length;
        }
};

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
        size_t visitedNodesLimit = settings.getBubbleDFSNodeLimit();
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
                                if (visited.size()>visitedNodesLimit)
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
                        << " from graph          \r").flush();
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
                                                bool nodeIsSingle=true;
                                                if (up.getNumLeftArcs()>1 ||up.getNumRightArcs()>1)
                                                        nodeIsSingle=false;
                                                if (nodeIsSingle){
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

                                                }else{
                                                        removeLinks(node, up);
                                                }
                                                if (upLast.isValid()&&upLast.getNodeKmerCov()<cutOffvalue){//&& node.getNodeKmerCov()/upLast.getNodeKmerCov()>3){
                                                        #ifdef DEBUG
                                                        mul=trueMult[abs( upLast.getNodeID())];
                                                        #endif
                                                        bool nodeIsSingle=true;
                                                        if (upLast.getNumLeftArcs()>1 ||upLast.getNumRightArcs()>1)
                                                                nodeIsSingle=false;
                                                        if (nodeIsSingle &&removeNode(upLast)){
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
                                                 bool noedIsSingle=true;
                                                if (down.getNumLeftArcs()>1 ||down.getNumRightArcs()>1)
                                                        noedIsSingle=false;
                                                if (noedIsSingle){
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
                                                }
                                                else{
                                                        removeLinks(node, down);
                                                }
                                                if (downLast.isValid()&&downLast.getNodeKmerCov()<cutOffvalue){// && node.getNodeKmerCov()/downLast.getNodeKmerCov()>3){
                                                        #ifdef DEBUG
                                                        mul=trueMult[abs( downLast.getNodeID())];
                                                        #endif
                                                        bool noedIsSingle=true;
                                                        if (downLast.getNumLeftArcs()>1 ||downLast.getNumRightArcs()>1)
                                                                noedIsSingle=false;
                                                        if( noedIsSingle&&removeNode(downLast)){
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
        cout << "Number of deleted nodes based on bubble detection: " << numOfDel << endl;
        if (numOfDel !=0)
                return true;
        return false;
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



