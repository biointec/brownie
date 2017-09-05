
#include "graph.h"
#include "settings.h"
#include "correctgraph.h"
#include "alignment.h"
#include <queue>
#include <map>
#include <limits>
#include <math.h>

using namespace std;


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
bool DBGraph:: removeErroneousComponents(vector<NodeID> untrustableComponentsNodes){
        cout << untrustableComponentsNodes.size() << " number of nodes are deleted in isolated components " <<endl;
        for (NodeID n:untrustableComponentsNodes)
                      removeNode(n);
}
