
#include "graph.h"
#include "settings.h"
#include "correctgraph.h"
#include "alignment.h"
#include <queue>
#include <map>
using namespace std;



bool DBGraph::deleteUnreliableNodes(double covCutoff, size_t maxMargLength){
      //bool changeIn1 = deleteExtraAttachedNodes(covCutoff,  maxMargLength);
      bool changeIn2 = connectSameMulNodes(covCutoff,  maxMargLength);
      return ( changeIn2);
}

bool DBGraph::deleteExtraAttachedNodes(double covCutoff, size_t maxMargLength){
        double tp=0, tn=0, fp=0,fn=0;
        size_t numOfDel=0;
        for ( NodeID lID = -numNodes; lID <= numNodes; lID++ ) {
                if (lID % OUTPUT_FREQUENCY == 0)
                        (cout << "Extracting node -" <<numNodes<< "/ "<<lID<<" /"<<numNodes
                        << " from graph        \r").flush();
                if (lID==0)
                        continue;
                SSNode node = getSSNode ( lID );

                if(!node.isValid())
                        continue;
                if (!checkNodeIsReliable(node))
                        continue;
                int mul = node.getExpMult(node.getAvgKmerCov());
                size_t rightArcs = node.getNumRightArcs();
                if(node.getExpMult(node.getAvgKmerCov())> node.getNumRightArcs())
                        continue;

                ArcIt it = node.rightBegin();
                SSNode currNode = getSSNode(it->getNodeID());
                while(it != node.rightEnd()) {
                        bool tip = currNode.getNumRightArcs()==0 && currNode.getNumLeftArcs()==1;
                        if (currNode.getAvgKmerCov()<covCutoff && tip ){
                                removeNode(currNode.getNodeID());
                                #ifdef DEBUG
                                if (trueMult[abs(currNode.getNodeID())]>0 ){
                                        fp++;
                                }
                                else{
                                        tp++;
                                }
                                #endif
                                numOfDel++;
                                break;
                        }
                        else{
                                #ifdef DEBUG
                                if (trueMult[abs(currNode.getNodeID())]>0){
                                        tn++;
                                }
                                else{
                                        fn++;
                                }
                                #endif
                        }
                        it++;
                }
        }
        cout<<endl;
        if (numOfDel>0)
                cout << "Number of deleted nodes in deleteExtraAttachedNodes: " << numOfDel << endl;
        #ifdef DEBUG
        if (numOfDel >0){
                cout<< "TP:     "<<tp<<"        TN:     "<<tn<<"        FP:     "<<fp<<"        FN:     "<<fn<<endl;
                cout << "Sensitivity: ("<<100*((double)tp/(double)(tp+fn))<<"%)"<<endl;
                cout<<"Specificity: ("<<100*((double)tn/(double)(tn+fp))<<"%)"<<endl;
        }
        #endif
        return (numOfDel>0);
}
bool DBGraph::connectSameMulNodes(double covCutoff, size_t maxMargLength){
        size_t secondFP = 0;
        size_t secondTP = 0;
        size_t numOfDel = 0;
        for ( NodeID lID = 1 ; lID <= numNodes; lID++ ) {
                SSNode node = getSSNode ( lID );
                if(!node.isValid())
                        continue;
                if (!checkNodeIsReliable(node))
                        continue;
                ArcIt it = node.rightBegin();
                SSNode nextReliableNode = getSSNode(it->getNodeID());
                bool found = false;
                while(it != node.rightEnd()) {
                        if (!checkNodeIsReliable(nextReliableNode))
                                it++;
                        if (nextReliableNode.getExpMult(nextReliableNode.getAvgKmerCov()) != node.getExpMult(node.getAvgKmerCov())
                                || nextReliableNode.getExpMult(node.getAvgKmerCov()) < 2)
                                it++;
                        found=true;
                        break;
                }
                if (found)
                {
                        ArcIt it = node.rightBegin();
                        while(it != node.rightEnd()) {
                                SSNode victim = getSSNode(it->getNodeID());
                                if(victim.getNodeID() != nextReliableNode.getNodeID()){
                                        double arcCov = node.getRightArc(victim.getNodeID())->getCoverage();
                                        if (arcCov <covCutoff)
                                                removeArc(node.getNodeID(),victim.getNodeID());
                                        numOfDel ++;
                                        break;
                                }
                                it++;
                        }
                        /*it = nextReliableNode.leftBegin();
                        while(it!=nextReliableNode.leftEnd()){
                                SSNode victim=getSSNode(it->getNodeID());
                                if (!victim.isValid())
                                        it++;
                                if(victim.getNodeID()!= node.getNodeID()){
                                        double arcCov = nextReliableNode.getRightArc(victim.getNodeID())->getCoverage();
                                        if (arcCov <covCutoff)
                                                removeArc(nextReliableNode.getNodeID(),victim.getNodeID());
                                                numOfDel++;
                                        break;
                                }
                                it++;
                        }*/
                }

        }

        if (numOfDel > 0)
                cout << "Number of deleted arcs in connectSameMulNodes: " << numOfDel << endl;
        return (numOfDel>0);

}

bool DBGraph::checkNodeIsReliable(SSNode node){
        if (node.getMarginalLength()< Kmer::getK()) // smaller nodes might not be correct, these ndoes can never be deleted
                return false;
        pair<int, pair<double,double> > result=nodesExpMult[abs( node.getNodeID())];
        double confidenceRatio=result.second.first;
        double inCorrctnessRatio=result.second.second;
        if (node.getNumRightArcs()<2)
                return false;
        if (1/confidenceRatio>.001)
                return false;
        if( 1/inCorrctnessRatio<.001)
                return false;
        return true;
}

bool cmssn(const SSNode& first, const SSNode & second ) {
    return first.getMarginalLength()>second.getMarginalLength();
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
void DBGraph::extractStatistic()
{
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
