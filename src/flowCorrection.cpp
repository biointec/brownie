
#include "graph.h"
#include "settings.h"
#include "correctgraph.h"
#include "alignment.h"
#include <queue>
#include <map>
using namespace std;

bool DBGraph::removeChimericLinksByFlow(double covCutoff, size_t maxMargLength){
        size_t numOfDel = 0;
        for ( NodeID lID = -numNodes ; lID <= numNodes; lID++ ) {
                if (lID ==0)
                        continue;
                SSNode node = getSSNode ( lID );
                if(!node.isValid())
                        continue;
                if (!checkNodeIsReliable(node, covCutoff, maxMargLength))
                        continue;
                if (node.getNumRightArcs() < 2)
                        continue;
                if (node.getNumRightArcs() < node.getExpMult(node.getAvgKmerCov()) )
                        continue;
                ArcIt it = node.rightBegin();
                bool found = false;
                SSNode nextReliableNode;
                while(it != node.rightEnd()) {
                        SSNode nextReliableNode = getSSNode(it->getNodeID());
                        if (checkNodeIsReliable(nextReliableNode , covCutoff, maxMargLength))
                        {
                                found=true;
                                break;
                        }
                        it++;
                }
                if (!found)
                        continue;

                ArcIt it = node.rightBegin();
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

bool DBGraph::checkNodeIsReliable(SSNode node, double covCutoff, size_t maxMargLength ){

        if (node.getMarginalLength() < maxMargLength || node.getAvgKmerCov() < covCutoff)
                return false;
        pair<int, pair<double,double> > result=nodesExpMult[abs( node.getNodeID())];
        double confidenceRatio = result.second.first;
        double inCorrctnessRatio = result.second.second;
        if (node.getExpMult(node.getAvgKmerCov())==0 )
                return false;
        if (node.getExpMult(node.getAvgKmerCov()) >4)
                return false;
        if (1/confidenceRatio>.001 ||1/inCorrctnessRatio<.001)
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
        double avg=0;
        for ( NodeID lID = 1; lID <= numNodes; lID++ ) {
                SSNode node = getSSNode ( lID );
                if(node.isValid() && node.getMarginalLength()) {
                        totalLength=totalLength+node.getMarginalLength();
                        nodeArray.push_back(node);
                }

        }
        sort(nodeArray.begin(), nodeArray.end(),cmssn );
        double sumOfCoverage=0;
        double sumOfMarginalLenght=0;
        double sizeLimit=0;
        sizeLimit= (totalLength*percentage)/100;
        size_t num  = 0;
        while(sumOfMarginalLenght < sizeLimit) {
                SSNode tempNode = nodeArray[num];
                sumOfMarginalLenght = sumOfMarginalLenght + tempNode.getMarginalLength();
                sumOfReadStcov = sumOfReadStcov + tempNode.getReadStartCov();
                sumOfCoverage = sumOfCoverage + tempNode.getKmerCov();  //tempNode.getExpMult();
                num++;
        }
        avg = sumOfReadStcov/sumOfMarginalLenght;
        return avg;
}
void DBGraph::extractStatistic(){
        double avg = getStartReadAvg();
        for ( NodeID lID = 1; lID <= numNodes; lID++ ) {

                SSNode node = getSSNode ( lID );
                if (!node.isValid())
                        continue;
                int nodeMultiplicity = 1;
                double maxProb = 0;
                double newValue = avg*nodeMultiplicity*node.getMarginalLength();
                double newProbability = Util::poissonPDF(node.getReadStartCov(),newValue);
                while(newProbability > maxProb) {
                        nodeMultiplicity++;
                        maxProb=newProbability;
                        newValue = avg*nodeMultiplicity*node.getMarginalLength();
                        newProbability = Util::poissonPDF(node.getReadStartCov(),newValue);
                }

                nodeMultiplicity--;
                double denominator = 0;
                double currentProb = maxProb;
                int i = 1;
                bool minus=true;
                do{
                        currentProb = newProbability;
                        double newValue = 0;
                        if (minus && nodeMultiplicity >i){
                                newValue = avg*(nodeMultiplicity-i)*node.getMarginalLength();
                                minus = false;
                        }
                        else{
                                newValue = avg*(nodeMultiplicity+i)*node.getMarginalLength();
                                minus = true;
                                i++;
                        }
                        newProbability = Util::poissonPDF(node.getReadStartCov(),newValue);
                        denominator = denominator+newProbability;
                }while(abs(newProbability-currentProb)> .000001|| i < 5);
                double confidenceRatio = maxProb/denominator;
                double expectToSee = avg*nodeMultiplicity*node.getMarginalLength();
                double observedprob = 0;
                if(node.getReadStartCov() < expectToSee)
                        observedprob=Util::poissonPDF(node.getReadStartCov(),expectToSee);
                else
                        observedprob=Util::poissonPDF(expectToSee,expectToSee);
                double inCorrctnessRatio = Util::poissonPDF(expectToSee,expectToSee)/observedprob;
                nodesExpMult[abs(node.getNodeID())] = make_pair(nodeMultiplicity , make_pair( confidenceRatio,inCorrctnessRatio));

        }


}
