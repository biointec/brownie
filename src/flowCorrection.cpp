
#include "graph.h"
#include "settings.h"
#include "correctgraph.h"
#include "alignment.h"
#include <queue>
#include <map>
#include <limits>
using namespace std;

bool DBGraph::removeChimericLinksByFlow(double covCutoff, size_t maxMargLength){
        size_t numOfDel = 0;
        for ( NodeID lID = -numNodes ; lID <= numNodes; lID++ ) {
                if (lID ==0)
                        continue;
                SSNode node = getSSNode ( lID );
                if(!node.isValid())
                        continue;
                if (!checkNodeIsReliable(node, covCutoff))
                        continue;

                pair<int, pair<double,double> > result=nodesExpMult[abs( node.getNodeID())];
                size_t nodeMultiplicity = result.first;
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
                        SSNode nextReliableNode = getSSNode(it->getNodeID());
                        if (checkNodeIsReliable(nextReliableNode , covCutoff))
                        {
                                pair<int, pair<double,double> > result=nodesExpMult[abs( nextReliableNode.getNodeID())];
                                size_t reliableNodeMultiplicity = result.first;
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

bool DBGraph::checkNodeIsReliable(SSNode node, double covCutoff){

        if ( node.getAvgKmerCov() < covCutoff)
                return false;
        pair<int, pair<double,double> > result=nodesExpMult[abs( node.getNodeID())];
        double confidenceRatio = result.second.first;
        double inCorrctnessRatio = result.second.second;
        size_t nodeMultiplicity = result.first;
        if ( nodeMultiplicity > 4 )
                return false;
        if (confidenceRatio <1000 || inCorrctnessRatio > 2 )
                return false;
        return true;

}

bool cmssn(const SSNode& first, const SSNode & second ) {
    return first.getMarginalLength() > second.getMarginalLength();
}

double DBGraph ::getStartReadAvg(){
        vector <SSNode> nodeArray;
        int percentage=5;
        double sumOfReadStcov=0;
        size_t totalLength=0;
        double avg = 0;
        for ( NodeID lID = 1; lID <= numNodes; lID++ ) {
                SSNode node = getSSNode ( lID );
                if(node.isValid() && node.getMarginalLength()) {
                        totalLength = totalLength + node.getMarginalLength()+settings.getK()-1;
                        nodeArray.push_back(node);
                }

        }
        sort(nodeArray.begin(), nodeArray.end(),cmssn );
        double sumOfMarginalLenght=0;
        double sizeLimit=0;
        sizeLimit= ( totalLength * percentage)/100;
        size_t numOfVisitedNodeLimit = 100;
        size_t num  = 0;
        while(sumOfMarginalLenght < sizeLimit  && numOfVisitedNodeLimit > num) {
                SSNode tempNode = nodeArray[num];
                sumOfMarginalLenght = sumOfMarginalLenght + tempNode.getMarginalLength() + settings.getK()-1;
                sumOfReadStcov = sumOfReadStcov + tempNode.getReadStartCov();
                num++;
        }
        cout.precision(3);
        avg = sumOfReadStcov/sumOfMarginalLenght;
        cout << "The avg of read start cov is "<< avg <<" whcih is calculated based on the " <<num <<" biggest nodes." <<endl;
        return avg;
}
void DBGraph::extractStatistic(){
        double avg = getStartReadAvg();
        for ( NodeID lID = 1; lID <= numNodes; lID++ ) {
                SSNode node = getSSNode ( lID );
                if (!node.isValid())
                        continue;
                double maxProb = .0001;
                size_t nodeMultiplicity = 0;
                double confidenceRatio = 0;
                double inCorrctnessRatio = std::numeric_limits<double>::max();
                size_t len = node.getMarginalLength() + settings.getK()-1;

                size_t i = 0;
                while ( i < 10 ) {
                        i++;
                        double expectedRSTCov = avg * i * len ;
                        double probabilityOfCurrMul = Util::poissonPDF(node.getReadStartCov(), expectedRSTCov);
                        if ( probabilityOfCurrMul > maxProb ){
                                maxProb = probabilityOfCurrMul;
                                nodeMultiplicity =i;
                        }
                }
                if (nodeMultiplicity ==0 ){
                        nodesExpMult[abs(node.getNodeID())] = make_pair(nodeMultiplicity , make_pair( confidenceRatio , inCorrctnessRatio));
                        continue;
                }

                double denominator = 0;
                double expectedRSTCov = avg * nodeMultiplicity * len ;
                double observedprob = Util::poissonPDF(node.getReadStartCov(),expectedRSTCov);
                i = 1;
                while( i < 10 ){
                        if (i != nodeMultiplicity){
                                double expectedRSTCov = avg * i * len;
                                double newProbability = Util::poissonPDF(node.getReadStartCov(),expectedRSTCov);
                                denominator = denominator + newProbability;
                        }
                        i = i+1;
                }

                if (denominator > 0 )
                        confidenceRatio = maxProb / denominator;
                expectedRSTCov = avg * nodeMultiplicity * len ;
                observedprob = Util::poissonPDF(node.getReadStartCov(),expectedRSTCov);
                inCorrctnessRatio = Util::poissonPDF(expectedRSTCov,expectedRSTCov)/observedprob;
                nodesExpMult[abs(node.getNodeID())] = make_pair(nodeMultiplicity , make_pair( confidenceRatio,inCorrctnessRatio));
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




