#include "alignment.h"
#include "graph.h"
#include "kmernode.h"
#include "readfile/fastafile.h"
#include "settings.h"
#include <gsl/gsl_randist.h>
#include <fstream>
#include <list>
#include <queue>
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
    double getNodePahtCoverage() {

        queue<SSNode> tempQueue=this->pathQueueNodes;
        double nodeCov;
        //pop root
        tempQueue.pop();
        while(!tempQueue.empty()) {
            SSNode node=tempQueue.front();
            tempQueue.pop();
            nodeCov=nodeCov+(node.getExpMult()/node.getMarginalLength());
        }
        nodeCov=nodeCov/numOfNode;
        return nodeCov;

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
    double getSingleExpStReadCov() {
        queue<SSNode> tempQueue=this->pathQueueNodes;
        double stReadoCov=0;
        int i=0;
        int length;
        //pop root
        tempQueue.pop();
        while(!tempQueue.empty()) {
            SSNode node=tempQueue.front();
            if(node.getNumLeftArcs()==1&&node.getNumRightArcs()==1) {
                tempQueue.pop();
                stReadoCov=stReadoCov+(node.getReadStartCov()/node.getMarginalLength());
                length=length+node.getMarginalLength();
                i++;
            } else {
                break;
            }
        }
        if(i!=0||length<31) {
            stReadoCov=(stReadoCov/i);
            return stReadoCov;
        }
        return-1;
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


                nodeCov=nodeCov+(node.getExpMult()/node.getMarginalLength());
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
    double getSingleNodeMulCertanity(const std::map<int, std::pair<int, std::pair<double, double> > >& nodesExpMult) {

        queue<SSNode> tempQueue=this->pathQueueNodes;
        double Certanity=0;
        if(tempQueue.size()>=2) {

            tempQueue.pop();
            int i=0;
            while(!tempQueue.empty()) {
                SSNode secondNode=tempQueue.front();
                if(secondNode.getNumLeftArcs()==1) {
                    tempQueue.pop();
                    pair<int, pair<double,double> > Result= nodesExpMult.at(abs( secondNode.getNodeID()));
                    Certanity=Certanity+Result.second.first;
                    i++;

                }
                else {
                    break;
                }

            }
            //the initial value is -1, so it should be added again
            if(i!=0) {
                Certanity=(Certanity)/(i);
                return Certanity;

            }

        }
        return -1;
    }


    double getSingleNodeIncorrectness(const std::map<int, std::pair<int, std::pair<double, double> > > &temp) {

        queue<SSNode> tempQueue=this->pathQueueNodes;
        double incorrectness=0;
        if(tempQueue.size()>=2) {

            tempQueue.pop();
            int i=0;
            while(!tempQueue.empty()) {
                SSNode secondNode=tempQueue.front();
                if(secondNode.getNumLeftArcs()==1) {
                    tempQueue.pop();
                    pair<int, pair<double,double> > Result=temp.at((int)abs(secondNode.getNodeID()));
                    //pair<int, pair<double,double> > Result= temp.at(abs( firstNode.getNodeID());
                    //temp.find()
                    incorrectness=incorrectness+Result.second.second;

                    i++;

                }
                else {
                    break;
                }

            }
            //the initial value is -1, so it should be added again
            if(i!=0) {
                incorrectness=(incorrectness)/(i);
                return incorrectness;

            }

        }
        return -1;
    }
    double getArcPathCoverage() {

        queue<SSNode> tempQueue=this->pathQueueNodes;
        double arcCov;
        if(tempQueue.size()>=2) {
            SSNode firstNode=tempQueue.front();
            tempQueue.pop();
            while(!tempQueue.empty()) {
                SSNode secondNode=tempQueue.front();
                tempQueue.pop();
                arcCov=arcCov+firstNode.getRightArc(secondNode.getNodeID())->getCoverage();
                firstNode=secondNode;

            }
            arcCov=arcCov/(numOfNode-1);
            return arcCov;
        } else {
            return 0;
        }


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

    priority_queue<pathStruct,vector<pathStruct>,comparator > MinHeap;
    int maxLength=settings.getK()*2;
    int numOfIncDel=0;
    int numOfCorDel=0;
    bool remove=false;



    updateCutOffValue(round);
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
                    //why I cannot do it?
                    //visitedArc.insert(leftNode.getRightArc())
                    pathStruct extendedPath;
                    extendedPath=leftPath;
                    //+rightNode.getMarginalLength()
                    if(extendedPath.pathLenght>=maxLength || rightNode.getNumRightArcs()==0)
                        continue;
                    extendedPath.addNodeToPaht(rightNode);
                    MinHeap.push(extendedPath);
                    if (visitedNodes.find(rightNode.getNodeID()) != visitedNodes.end()) {
                        pathStruct prevPath=pathDic.at(rightNode.getNodeID());
                        // this condition can be removed later, it check to have a parallel path which starts from same node.
                        //to avoid below example (ab is common in both path, we need only one common node)
                        //abcdef
                        //abmnktgf
                        if (prevPath.firstNode.getNodeID()!=extendedPath.firstNode.getNodeID()) {
                            double lengthpro=0;
                            if(prevPath.pathLenght!=0)
                                lengthpro=extendedPath.pathLenght/prevPath.pathLenght;
                            if (lengthpro>.8 && lengthpro<1.2) {

                                double prevNodeCov=prevPath.getSingleNodePahtCoverage();
                                double prevAcrcCov=prevPath.getSingleArcPathCoverage();
                                double exteNodeCov=extendedPath.getSingleNodePahtCoverage();
                                double exteArcCov=extendedPath.getSingleArcPathCoverage();


                                double  prevProbalility=  gsl_cdf_gaussian_P((prevNodeCov-estimatedKmerCoverage)/estimatedMKmerCoverageSTD,1);
                                double  exteProbalility=  gsl_cdf_gaussian_P((exteNodeCov-estimatedKmerCoverage)/estimatedMKmerCoverageSTD,1);
                                double prevProbArcCov=gsl_cdf_gaussian_P((prevAcrcCov-estimatedArcCoverageMean)/estimatedArcCoverageSTD,1);
                                double exteProbArcCov=gsl_cdf_gaussian_P((exteArcCov-estimatedArcCoverageMean)/estimatedArcCoverageSTD,1);

                                bool preReliable=false;
                                bool find=false;
                                queue<SSNode> tempQueue;

                                if(prevNodeCov!=-1 && prevAcrcCov!=-1 && exteNodeCov!=-1&& exteArcCov!=-1) {
                                    if(prevProbalility<.1||exteProbalility<.1) {

                                        if (((prevProbArcCov<.1&&prevProbalility<.1&&prevNodeCov<this->redLineValueCov) || prevNodeCov<this->certainVlueCov)&&(exteProbArcCov>.1&&exteProbalility>.1) ) {
                                            preReliable=false;
                                            find=true;
                                            tempQueue=prevPath.pathQueueNodes;
                                        }
                                        if (((exteProbArcCov<.1&&exteProbalility<.1&&exteNodeCov<this->redLineValueCov)|| exteNodeCov<this->certainVlueCov)&& (prevProbArcCov>.1&&prevProbalility>.1)  ) {
                                            preReliable=true;
                                            find=true;
                                            tempQueue=extendedPath.pathQueueNodes;
                                        }
                                    }
                                } else {
                                    if(prevNodeCov!=-1&&prevAcrcCov!=-1) {
                                        if ((prevProbArcCov<.1&&prevProbalility<.1&&prevNodeCov<this->redLineValueCov )&&(prevNodeCov<this->certainVlueCov)) {
                                            preReliable=false;
                                            find=true;
                                            tempQueue=prevPath.pathQueueNodes;
                                        }

                                    } else {
                                        if(exteNodeCov!=-1&& exteArcCov!=-1) {
                                            if ((exteProbArcCov<.1&&exteProbalility<.1&&exteNodeCov<this->redLineValueCov)&&(exteNodeCov<this->certainVlueCov)) {
                                                preReliable=true;
                                                find=true;
                                                tempQueue=extendedPath.pathQueueNodes;
                                            }

                                        }
                                    }
                                }
                                if(nodesExpMult.size()>0 && !find) {
                                    pair<int, pair<double,double> > resultR=nodesExpMult[abs( extendedPath.rootNode.getNodeID())];
                                    if(resultR.first==1 && (resultR.second.second/resultR.second.first)<1) {
                                        pair<int, pair<double,double> > resultE=nodesExpMult[abs( extendedPath.firstNode.getNodeID())];
                                        pair<int, pair<double,double> > resultP=nodesExpMult[abs( prevPath.firstNode.getNodeID())];
                                        find =true;
                                        if ((resultE.second.second/resultE.second.first)<resultP.second.second/resultP.second.first) {
                                            preReliable=false;
                                            tempQueue=prevPath.pathQueueNodes;

                                        } else {
                                            preReliable=true;
                                            tempQueue=extendedPath.pathQueueNodes;
                                        }
                                    }

                                }


                                if(find) {
                                    if ( totalRAM<getMemoryUsedByProc())
                                        totalRAM=getMemoryUsedByProc();
                                    if (!tempQueue.empty()) {
                                        SSNode tempNode=tempQueue.front();
                                        tempQueue.pop();
                                        tempNode=tempQueue.front();
                                        while(!tempQueue.empty()&& tempNode.getNumLeftArcs()==1 && tempNode.getNumRightArcs()==1) {
                                            tempQueue.pop();
                                            if(trueMult[abs(tempNode.getNodeID())]>0)
                                                numOfIncDel++;
                                            tempNode=tempQueue.front();
                                        }
                                        if(preReliable) {
                                            numOfCorDel=numOfCorDel+extendedPath.removeSingleNodes();


                                            remove=true;
                                            break;

                                        } else {
                                            numOfCorDel=numOfCorDel+ prevPath.removeSingleNodes();


                                            remove=true;
                                            break;

                                        }

                                    }
                                }

                            }
                        }

                    } else {
                        visitedNodes.insert(rightNode.getNodeID());

                        pathDic[rightNode.getNodeID()]=extendedPath;


                    }


                }
                //}
            }
        }


    }
    if(numOfCorDel>0||numOfIncDel>0) {
        cout<<"number of correctly deleted nodes based on bubble detection: "<<numOfCorDel<<endl;
        cout<<"number of incorrectly deleted nodes based on bubble detection: "<<numOfIncDel<<endl;
    }
    if (numOfCorDel ==0 &&numOfIncDel==0)
        remove= false;

    return remove;

}
