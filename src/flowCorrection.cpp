
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
        size_t maxNumOfVisitedNode = 100;
        size_t minNumOfVisitedNode = 30 < nodeArray.size() ? 30:nodeArray.size() ;
        size_t num  = 0;
        while(sumOfMarginalLenght < sizeLimit  && maxNumOfVisitedNode > num  || num < minNumOfVisitedNode ) {
                num++;
                SSNode tempNode = nodeArray[num-1];
                sumOfMarginalLenght = sumOfMarginalLenght + tempNode.getMarginalLength() + settings.getK()-1;
                sumOfReadStcov = sumOfReadStcov + tempNode.getReadStartCov();
        }
        avg = sumOfReadStcov/sumOfMarginalLenght;
        size_t i = 0;
        while(i < num) {

                SSNode tempNode = nodeArray[i];
                double len = tempNode.getMarginalLength() + settings.getK()-1;
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




