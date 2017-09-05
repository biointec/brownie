
#include "graph.h"
#include "settings.h"
#include "correctgraph.h"
#include "alignment.h"
#include <queue>
#include <map>
#include <limits>
#include <math.h>
using namespace std;

bool DBGraph::removeChimericLinksByFlow(double covCutoff, size_t maxMargLength){


        map<NodeID, size_t> nodesExpMult;
        extractStatistic (nodesExpMult, maxMargLength);

        size_t numOfDel = 0;
        for ( NodeID lID = -numNodes ; lID <= numNodes; lID++ ) {
                if (lID ==0)
                        continue;
                SSNode node = getSSNode ( lID );
                if(!node.isValid())
                        continue;
                if (!checkNodeIsReliable(node, covCutoff, maxMargLength,nodesExpMult ))
                        continue;
                if (nodesExpMult[abs( node.getNodeID())] == 0)
                        continue;

                size_t nodeMultiplicity =  nodesExpMult[abs( node.getNodeID())];
                if (nodeMultiplicity ==0)
                        continue;
                if (node.getNumRightArcs() < 2)
                        continue;
                if (node.getNumRightArcs() < nodeMultiplicity )
                        continue;
                //it tries to find two consecutive nodes with the same node Multiplicity.
                //So before connecting them there should be a Chimeric Link which needs to be removed
                ArcIt it = node.rightBegin();
                bool found = false;
                SSNode nextReliableNode;
                while(it != node.rightEnd()) {
                        nextReliableNode = getSSNode(it->getNodeID());
                        if (checkNodeIsReliable(nextReliableNode , covCutoff, maxMargLength, nodesExpMult))
                        {
                                size_t reliableNodeMultiplicity = nodesExpMult[abs( nextReliableNode.getNodeID())];
                                if (nodeMultiplicity == reliableNodeMultiplicity){
                                        found =  true;
                                        break;
                                }
                        }
                        it++;
                }
                if (!found)
                        continue;

                it = node.rightBegin();
                while(it != node.rightEnd()) {
                        SSNode victim = getSSNode(it->getNodeID());
                        if(victim.getNodeID() != nextReliableNode.getNodeID()){
                                double arcCov = node.getRightArc(victim.getNodeID())->getCoverage();
                                if (arcCov <covCutoff){
                                        removeArc(node.getNodeID(),victim.getNodeID());
                                        //cout << "The link between "<< node.getNodeID() << " : " <<victim.getNodeID() << " removed, the trustable node is "<<nextReliableNode.getNodeID() <<endl;
                                        numOfDel ++;
                                        break;
                                }
                        }
                        it++;
                }

        }
        cout << "Number of deleted arcs in flow correction: " << numOfDel << endl;
        return (numOfDel > 0);
}

bool DBGraph::checkNodeIsReliable(SSNode node, double covCutoff, size_t maxMargLength , map<NodeID, size_t> &nodesExpMult){

        if ( node.getAvgKmerCov() < covCutoff)
                return false;
        size_t nodeMultiplicity = nodesExpMult[abs( node.getNodeID())];
        if (node.getMarginalLength() < maxMargLength)
                return false;
        if ( nodeMultiplicity > 4 || nodeMultiplicity ==0 )
                return false;

        return true;

}


bool cmssn(const SSNode& first, const SSNode & second ) {
    return first.getMarginalLength() > second.getMarginalLength();
}

void DBGraph ::getReadStartCovAvg( double &avg, double &variance){
        vector <SSNode> nodeArray;
        int percentage=5;
        double sumOfReadStcov=0;
        size_t totalLength=0;

        for ( NodeID lID = 1; lID <= numNodes; lID++ ) {
                SSNode node = getSSNode ( lID );
                if(node.isValid() && node.getMarginalLength()) {
                        totalLength = totalLength + node.getMarginalLength()+settings.getK()-1;
                        nodeArray.push_back(node);
                }
        }
        sort(nodeArray.begin(), nodeArray.end(),cmssn );
        double sumOfMarginalLenght = 0;
        double sizeLimit = 0;
        sizeLimit = ( totalLength * percentage)/100;
        size_t numOfVisitedNodeLimit = 100;
        size_t minNumOfsamples = 30;
        size_t num  = 0;
        while(sumOfMarginalLenght < sizeLimit  && numOfVisitedNodeLimit > num || num < minNumOfsamples) {
                num++;
                SSNode tempNode = nodeArray[num-1];
                sumOfMarginalLenght = sumOfMarginalLenght + tempNode.getMarginalLength() + settings.getK()-1;
                sumOfReadStcov = sumOfReadStcov + tempNode.getReadStartCov();
        }
        avg = sumOfReadStcov/sumOfMarginalLenght;
        size_t i = 0;
        while(i < num) {

                SSNode tempNode = nodeArray[i];
                double len =tempNode.getMarginalLength() + settings.getK()-1;
                variance = variance +(tempNode.getReadStartCov()/len - avg )* (tempNode.getReadStartCov()/len - avg );
                sumOfReadStcov = sumOfReadStcov + tempNode.getReadStartCov();
                i ++;
        }
        variance  = variance /(num -1);
        cout.precision(3);
        cout << "The avg of read start cov is "<< avg <<" whcih is calculated based on the " <<num <<" biggest nodes." <<endl;
}
void DBGraph::extractStatistic(map<NodeID, size_t> &nodesExpMult, size_t maxMargLength ){
        double avg =0 , variance = 0 ;
        getReadStartCovAvg(avg,variance);
        variance = round(variance * 10000000) /10000000.0;
        for ( NodeID lID = 1; lID <= numNodes; lID++ ) {
                SSNode node = getSSNode ( lID );
                if (!node.isValid())
                        continue;
                size_t len = node.getMarginalLength() + settings.getK()-1;
                double readStarCov = node.getReadStartCov();
                if ( len < maxMargLength || readStarCov ==0){
                        nodesExpMult[abs(node.getNodeID())] = 0;
                        continue;
                }
                double minValue = 3.5;
                size_t nodeMultiplicity = 0;
                double sigma = sqrt(variance )*len;
                size_t i = 0;
                int numOfassociateCurves = 0;
                while ( i < 10 ) {
                        i++;
                        double expectedRSTCov = avg * i * len ;
                        //double probabilityOfCurrMul = Util::poissonPDF(node.getReadStartCov(), expectedRSTCov);
                        double zvalue = (node.getReadStartCov() - (expectedRSTCov ))/(sigma);
                        if ( abs(zvalue) < minValue ){
                                minValue = abs(zvalue);
                                nodeMultiplicity =i;
                                numOfassociateCurves ++;
                        }
                }
                if (nodeMultiplicity ==0 || numOfassociateCurves > 1){
                        nodesExpMult[abs(node.getNodeID())] = 0;
                        continue;
                }
                if (abs(minValue) <3)
                        nodesExpMult[abs(node.getNodeID())] = nodeMultiplicity;
        }
}

size_t DBGraph:: findComponentsInGraph(size_t minSize)
{
        cout << "Finding disjoint components in the graph ..." <<endl;
        int srcID=1;
        int i=0;
        size_t numberOfcomponents =0;
        set<NodeID> nodesHandled;
        for (i =srcID ; i <=getNumNodes(); i++){
                set<NodeID> currentSetNodes;
                size_t componentSize =0;
                if (!getSSNode(i).isValid())
                        continue;
                if (nodesHandled.find(i) != nodesHandled.end())
                        continue;
                if (nodesHandled.find(i) != nodesHandled.end())
                        continue;
                srcID=i;
                multimap<size_t, NodeID> nodeDepth;     // map of all nodes in the local graph and its depth
                nodeDepth.insert(pair<size_t, NodeID>(0, srcID));
                // nodes that were already handled

                while (!nodeDepth.empty()) {
                        // get and erase the current node
                        multimap<size_t, NodeID>::iterator
                        e = nodeDepth.begin();
                        size_t thisDepth = e->first;
                        NodeID thisID = e->second;
                        nodeDepth.erase(e);
                        // if the node was already handled, skip
                        if (nodesHandled.find(thisID) != nodesHandled.end())
                                continue;
                        if (nodesHandled.find(-thisID) != nodesHandled.end())
                                continue;
                        // mark this node as handled
                        nodesHandled.insert(thisID);
                        currentSetNodes.insert(thisID);
                        SSNode thisNode = getSSNode(thisID);
                        componentSize = componentSize + thisNode.getMarginalLength();
                        for (ArcIt it = thisNode.rightBegin(); it != thisNode.rightEnd(); it++) {
                                SSNode rNode = getSSNode(it->getNodeID());
                                if (!rNode.isValid())
                                        continue;
                                if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                        continue;
                                nodeDepth.insert(pair<size_t, NodeID>(thisDepth + thisNode.getMarginalLength(), it->getNodeID()));
                        }
                        for (ArcIt it = thisNode.leftBegin(); it != thisNode.leftEnd(); it++) {
                                SSNode lNode = getSSNode(it->getNodeID());
                                if (!lNode.isValid())
                                        continue;
                                if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                        continue;
                                nodeDepth.insert(pair<size_t, NodeID>(thisDepth + thisNode.getMarginalLength(), it->getNodeID()));
                        }
                }
                if ( componentSize >minSize && currentSetNodes.size() >1)
                        numberOfcomponents ++;

                currentSetNodes.clear();
        }
        cout << "Number of disjoint components in the graph with more than 1 nodes and larger than : " <<minSize << " is " <<numberOfcomponents <<endl;
        return numberOfcomponents;
}




