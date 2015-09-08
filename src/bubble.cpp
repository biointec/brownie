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
#include "sparseSA.hpp"
#include "fasta.hpp"
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
};

bool DBGraph::bubbleDetection(int round) {
        cout<<"*********************<<Bubble Detection starts>>......................................... "<<endl;
        priority_queue<pathStruct,vector<pathStruct>,comparator > MinHeap;
        int maxLength=settings.getK()*2;
        size_t numOfDel=0;
        bool remove=false;
        size_t TP=0,TN=0,FP=0,FN=0;
        for ( NodeID lID = -numNodes; lID <= numNodes; lID++ ) {
                if ( lID == 0 ) {
                        continue;
                }
                SSNode leftNode = getSSNode ( lID );
                if ( !leftNode.isValid() ) {
                        continue;
                }
                if(leftNode.getNumRightArcs()<2) {
                        continue;
                }
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
}
bool DBGraph::removeBubble(SSNode &prevFirstNode ,SSNode& extendFirstNode,size_t &TP,size_t &TN,size_t &FP,size_t &FN,size_t & numOfDel ){
        bool preIsSingle=true;
        bool exteIsSingle=true;
        if (prevFirstNode.getNumLeftArcs()>1 ||prevFirstNode.getNumRightArcs()>1)
                preIsSingle=false;
        if(extendFirstNode.getNumLeftArcs()>1||extendFirstNode.getNumRightArcs()>1)
                exteIsSingle=false;
        double preCov= prevFirstNode.getNodeKmerCov();// prevFirstNode.getExpMult()/prevFirstNode.getMarginalLength();
        double extCov=extendFirstNode.getNodeKmerCov();//extendFirstNode.getExpMult()/extendFirstNode.getMarginalLength();
        if(preIsSingle && exteIsSingle) {
                if ( preCov<this->redLineValueCov && preCov<extCov)  {
                        #ifdef DEBUG
                        if (trueMult[abs( prevFirstNode.getNodeID())]>0)
                                FP++;
                        else
                                TP++;
                        #endif
                        if(removeNode(prevFirstNode)){
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
                if (extCov<this->redLineValueCov && extCov<preCov ) {
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
        }
        if(preIsSingle && !exteIsSingle) {
                if (preCov<this->safeValueCov ) {
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
        }
        if(exteIsSingle&&!preIsSingle ) {
                if ((extCov<this->safeValueCov)) {
                        #ifdef DEBUG
                        if(trueMult[abs( extendFirstNode.getNodeID())]>0)
                                FP++;
                        else
                                TP++;
                        #endif
                        if(removeNode(extendFirstNode)){
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
        }
        //when the starting nodes of both parallel paths has more than one ingoing or arcgoing arcs
        if (!exteIsSingle&&!preIsSingle)
                return removeNotSingleBublles(prevFirstNode,extendFirstNode, TP,TN,FP,FN, numOfDel);
        return false;
}
bool DBGraph:: removeNotSingleBublles( SSNode &prevFirstNode ,SSNode& extendFirstNode, size_t &TP,size_t &TN,size_t &FP,size_t &FN,size_t & numOfDel){
        double preCov=prevFirstNode.getNodeKmerCov();// prevFirstNode.getExpMult()/prevFirstNode.getMarginalLength();
        double extCov=extendFirstNode.getNodeKmerCov();//extendFirstNode.getExpMult()/extendFirstNode.getMarginalLength();
        if (preCov<extCov && preCov<this->certainVlueCov){
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
        if (extCov<preCov && extCov<this->certainVlueCov){
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
