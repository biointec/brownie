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

using namespace std;

void DBGraph::parseReads(size_t thisThread,
                         vector<string>& readBuffer) const
{
        for (size_t i = 0; i < readBuffer.size(); i++) {
                const string& read = readBuffer[i];

                KmerIt it(read);
                if (!it.isValid())
                        continue;

                // increase the read start coverage (only for the first valid kmer)
                NodePosPair result = table->find(it.getKmer());
                if (result.getNodeID() != 0) {
                        SSNode node = getSSNode(result.getNodeID());
                        node.setReadStartCov(node.getReadStartCov()+1);
                }

                NodeID prevID = 0;
                for (KmerIt it(read); it.isValid(); it++ ) {
                        Kmer kmer = it.getKmer();
                        NodePosPair result = table->find(kmer);
                        if (!result.isValid()) {
                                prevID = 0;
                                continue;
                        }

                        // we've found the kmer, increase the node coverage
                        NodeID thisID = result.getNodeID();
                        SSNode node = getSSNode(thisID);
                        node.incKmerCov();

                        // if the previous node was valid and different, increase the arc coverage
                        if ((prevID != 0) && (prevID != thisID)) {
                                getSSNode(prevID).getRightArc(thisID)->incReadCov();
                                getSSNode(thisID).getLeftArc(prevID)->incReadCov();
                        }

                        if (it.hasRightOverlap())
                                prevID = thisID;
                        else
                                prevID = 0;
                }
        }
}

void DBGraph::workerThread(size_t thisThread, LibraryContainer* inputs)
{
        // local storage of reads
        vector<string> myReadBuf;

        while (inputs->getReadChunk(myReadBuf, settings.getThreadWorkSize()))
                parseReads(thisThread, myReadBuf);
}

void DBGraph::countNodeandArcFrequency(LibraryContainer &inputs)
{
        const unsigned int& numThreads = settings.getNumThreads();

        cout << "Populating table... ";
        cout.flush();
        table = new KmerNodeTable(settings, numNodes);
        table->populateTable (nodes);
        cout << "done" << endl;

        cout << "Number of threads: " << numThreads << endl;
        cout << "Counting node and arc frequency: " << endl;

        inputs.initiateReadThreading();

        // start worker threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&DBGraph::workerThread, this, i, &inputs);

        inputs.threadReads();

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        inputs.finalizeReadThreading();

        delete table;
}

// ============================================================================
// PRIVATE NODE COVERAGE / MULTIPLICITY
// ============================================================================

//comment by mahdi


bool cmssn(const SSNode& first, const SSNode & second ) {

    return first.getMarginalLength()>second.getMarginalLength();


}

void DBGraph::extractStatistic(int round) {
        vector <SSNode> nodeArray;
        int percentage=5;
        double sumOfReadStcov=0;

        size_t totalLength=0;

        double coverage=0;
        double StandardErrorMean=0;
        double avg=0;
        for ( NodeID lID = 1; lID <= numNodes; lID++ ) {

                SSNode node = getSSNode ( lID );
                if(node.isValid() && node.getMarginalLength()) {
                        totalLength=totalLength+node.getMarginalLength();
                        nodeArray.push_back(node);

                }

        }
        sort(nodeArray.begin(), nodeArray.end(),cmssn );
        size_t i=0;
        double sumOfCoverage=0;
        double sumOfMarginalLenght=0;
        double sizeLimit= ((totalLength*percentage)/100)>settings.getGenomeSize()?settings.getGenomeSize():((totalLength*percentage)/100);
        while(sumOfMarginalLenght<sizeLimit) {
                SSNode tempNode=nodeArray[i];
                sumOfMarginalLenght=sumOfMarginalLenght+tempNode.getMarginalLength();
                sumOfReadStcov=sumOfReadStcov+tempNode.getReadStartCov();
                sumOfCoverage=sumOfCoverage+  tempNode.getKmerCov();  //tempNode.getExpMult();
                i++;
        }
        avg=sumOfReadStcov/sumOfMarginalLenght;
        estimatedKmerCoverage=(sumOfCoverage/sumOfMarginalLenght);
        size_t num=i;
        double sumOfSTD=0;
        i=0;
        while(i<num) {
                SSNode tempNode=nodeArray[i];
                sumOfSTD=sumOfSTD+((tempNode.getNodeKmerCov())-estimatedKmerCoverage)*((tempNode.getNodeKmerCov())-estimatedKmerCoverage);
                i++;
        }
        if(num>0)
                StandardErrorMean=sqrt( sumOfSTD/num);
        if (round==0)
                estimatedMKmerCoverageSTD=StandardErrorMean;//          sqrt(num)*std;

        for ( NodeID lID = 1; lID <= numNodes; lID++ ) {

                SSNode node = getSSNode ( lID );
                if (!node.isValid())
                        continue;
                int nodeMultiplicity=1;
                double numerator=0;
                double newValue=avg*nodeMultiplicity*node.getMarginalLength();
                double newProbability=gsl_ran_poisson_pdf(node.getReadStartCov(),newValue);

                while(newProbability>numerator) {
                        nodeMultiplicity++;
                        numerator=newProbability;
                        newValue=avg*nodeMultiplicity*node.getMarginalLength();
                        newProbability=gsl_ran_poisson_pdf(node.getReadStartCov(),newValue);
                }

                nodeMultiplicity--;

                double denominator=0;
                double currentProb=numerator;
                double inCorrctnessRatio=0;
                int i=nodeMultiplicity-1;
                if(i>0) {
                        double newValue=avg*i*node.getMarginalLength();
                        newProbability=gsl_ran_poisson_pdf(node.getReadStartCov(),newValue);
                        denominator=denominator+newProbability;
                        i--;
                        while(i>=0&& abs(newProbability-currentProb)> .00000001) {
                                currentProb=newProbability;
                                newValue=avg*i*node.getMarginalLength();
                                newProbability=gsl_ran_poisson_pdf(node.getReadStartCov(),newValue);
                                denominator=denominator+newProbability;
                                i--;
                        }
                }
                i=nodeMultiplicity+1;
                newValue=avg*i*node.getMarginalLength();
                newProbability=gsl_ran_poisson_pdf(node.getReadStartCov(),newValue);
                denominator=denominator+newProbability;
                i++;
                while(i>0&& abs(newProbability-currentProb)> .00000001) {
                        currentProb=newProbability;
                        newValue=avg*i*node.getMarginalLength();
                        newProbability=gsl_ran_poisson_pdf(node.getReadStartCov(),newValue);
                        denominator=denominator+newProbability;
                        i++;
                }

                newValue=avg*nodeMultiplicity*node.getMarginalLength();
                if(node.getReadStartCov()<newValue)
                        newProbability=gsl_ran_poisson_pdf(node.getReadStartCov(),newValue);
                else
                        newProbability=gsl_ran_poisson_pdf(newValue,newValue);
                inCorrctnessRatio=gsl_ran_poisson_pdf(newValue,newValue)/newProbability;
                /*if (nodeMultiplicity==1) {
                        int shiftedStReadCov=leftNode.getReadStartCov()+avg*leftNode.getMarginalLength();
                        //i=2
                        newValue=2*avg*leftNode.getMarginalLength();
                        numerator=gsl_ran_poisson_pdf(shiftedStReadCov,newValue);

                        //for i=1
                        newValue=avg*leftNode.getMarginalLength();
                        denominator=gsl_ran_poisson_pdf(shiftedStReadCov,newValue);
                        i=3;
                        newValue=avg*i*leftNode.getMarginalLength();
                        newProbability=gsl_ran_poisson_pdf(shiftedStReadCov,newValue);
                        currentProb=numerator;

                        while(abs(newProbability-currentProb)> .00000001) {
                                i++;
                                denominator=denominator+newProbability;
                                currentProb=newProbability;
                                newValue=avg*i*leftNode.getMarginalLength();
                                newProbability=gsl_ran_poisson_pdf(shiftedStReadCov,newValue);
                        }
                }*/
                if (newProbability/gsl_ran_poisson_pdf(newValue,newValue)<.02)
                        node.setExpMult(0);
                else
                        node.setExpMult(nodeMultiplicity);

                double confidenceRatio=numerator/denominator;
                nodesExpMult[node.getNodeID()]=make_pair(nodeMultiplicity,make_pair( confidenceRatio,inCorrctnessRatio));

        }

}


bool DBGraph::deleteUnreliableNodes(int round){
        double tp=0, tn=0, fp=0,fn=0;
        size_t numOfDel=0;
        bool change=false;
        for ( NodeID lID =-numNodes; lID <= numNodes; lID++ ) {
                if (lID==0)
                        continue;
                SSNode leftNode = getSSNode ( lID );
                if(!leftNode.isValid())
                        continue;

                if(leftNode.getNumRightArcs()- leftNode.getExpMult()>0) {
                        ArcIt it = leftNode.rightBegin();
                        while(it != leftNode.rightEnd()) {
                                SSNode currNode = getSSNode(it->getNodeID());
                                if ((currNode.getNodeKmerCov()<cutOffvalue|| (currNode.getNumRightArcs()==0&&currNode.getNumLeftArcs()==1))&&(currNode.getExpMult()==0)){
                                        if (trueMult[abs(currNode.getNodeID())]>0){
                                                fp++;
                                        }
                                        else{
                                                tp++;
                                        }
                                        change= removeNode(currNode);
                                        numOfDel++;
                                        break;
                                }
                                else{
                                        if (trueMult[abs(currNode.getNodeID())]>0){
                                                tn++;
                                        }
                                        else{
                                                fn++;
                                        }
                                }
                                it++;
                        }

                        if(leftNode.getExpMult()==0&& leftNode.getNodeKmerCov()<cutOffvalue&& !change){
                                if (trueMult[abs(leftNode.getNodeID())]>0){
                                        fp++;
                                }
                                else{
                                        tp++;;
                                }
                                numOfDel++;
                                removeNode(leftNode);
                        }else{
                                if (trueMult[abs(leftNode.getNodeID())]>0){
                                        tn++;
                                }
                                else{
                                        fn++;
                                }

                        }

                }
        }

        cout<<"number of deleted nodes in deleteUnreliableNodes rutin:: "<<numOfDel<<endl;
        #ifdef DEBUG
        cout<< "TP:     "<<tp<<"        TN:     "<<tn<<"        FP:     "<<fp<<"        FN:     "<<fn<<endl;
        cout << "Sensitivity: ("<<100*((double)tp/(double)(tp+fn))<<"%)"<<endl;
        cout<<"Specificity: ("<<100*((double)tn/(double)(tn+fp))<<"%)"<<endl;
        #endif
        return change;
}
/*
bool DBGraph::deleteUnreliableNodes( int round) {
        cout<<"*********************<<Delete Unreliable Nodes starts>>......................................... "<<endl;
        bool modify=false;
        int numberOfDel=0;
        vector<double> tempArray1;
        vector<double> tempArray2;
        int i=0;
        for ( NodeID lID =1; lID <= numNodes; lID++ ) {
                SSNode leftNode = getSSNode ( lID );
                if(!leftNode.isValid())
                        continue;
                i++;
                pair<int, pair<double,double> > result=nodesExpMult[abs( leftNode.getNodeID())];
                if (!std::isinf( result.second.first))
                        tempArray1.push_back(result.second.second/result.second.first);
                else
                        tempArray1.push_back(0);
                tempArray2.push_back(result.second.second);
        }
        sort(tempArray1.begin(), tempArray1.end());
        sort(tempArray2.begin(), tempArray2.end());
        double max= tempArray1.at(floor(i*.80));
        double min=tempArray2.at(floor(i*.90));
        max=max<5?max:5;
        min=min>100?min:100;
        for ( NodeID lID =1; lID <= numNodes; lID++ ) {
                SSNode leftNode = getSSNode ( lID );
                if(!leftNode.isValid())
                        continue;
                pair<int, pair<double,double> > result=nodesExpMult[abs( leftNode.getNodeID())];
                double confidenceRatio=result.second.first;
                double inCorrctnessRatio=result.second.second;
                double nodeMultiplicity=result.first;
                bool change=false;
                double ratio=0;

                if (std::isinf(confidenceRatio)) {
                        ratio=0;
                }
                else {
                        ratio=inCorrctnessRatio/confidenceRatio;
                }
                if((inCorrctnessRatio>min&& leftNode.getMarginalLength()<this->maxNodeSizeToDel)) {

                        if (leftNode.getNumLeftArcs()==0)
                                continue;
                        if (leftNode.getNodeKmerCov()<this->redLineValueCov)
                                change=removeNode(leftNode);
                        if(change) {
                                modify=true;
                                numberOfDel++;
                        }

                } else {

                        if(ratio<.5&& leftNode.getMarginalLength()<this->maxNodeSizeToDel) {
                                if(leftNode.getNumRightArcs()- nodeMultiplicity>=1) {
                                        do {

                                                if (leftNode.getNodeKmerCov()<this->redLineValueCov)
                                                        change=deleteExtraRightLink(leftNode, round);

                                        } while ((change&& leftNode.getNumRightArcs()- nodeMultiplicity>=1 ));
                                }

                                if(leftNode.getNumLeftArcs()-nodeMultiplicity>=1) {
                                        do {
                                                change=deleteExtraLeftLink(leftNode, round);

                                        } while ((change&& (leftNode.getNumLeftArcs()- nodeMultiplicity)>=1 ));
                                }
                                if(change) {
                                        modify=true;
                                        numberOfDel++;
                                }
                        }

                }

        }

        cout<<"number of noeds deleted in guessNodeMultiplicity procedure is:	"<<numberOfDel<<endl;
        return modify;
}

*/

bool DBGraph::deleteExtraRightLink(SSNode rootNode, int round) {
        bool change=false;
        ArcIt it = rootNode.rightBegin();
        SSNode bestLNode = getSSNode(it->getNodeID());
        SSNode trustedNode=getSSNode(it->getNodeID());

        pair<int, pair<double,double> > bestResult=nodesExpMult[abs( bestLNode.getNodeID())];

        pair<int, pair<double,double> > trustedResult=bestResult;
        it++;
        while(it != rootNode.rightEnd()) {

                SSNode currNode = getSSNode(it->getNodeID());
                pair<int, pair<double,double> > curResult= nodesExpMult[abs( currNode.getNodeID())];

                if (bestResult.second.second<curResult.second.second) {
                        bestLNode=currNode;
                        bestResult=curResult;
                }
                if(trustedResult.second.second>curResult.second.second) {
                        trustedNode=currNode;
                        trustedResult=curResult;
                }
                it++;
        }

        // if (((bestResult.second.second/rootResult.second.second)>10&&(rootResult.second.first/ bestResult.second.first)>10 )) {
        double ratio=bestResult.second.first/bestResult.second.second;
        if (ratio<.01) {
                for ( ArcIt it3 = bestLNode.rightBegin(); it3 != bestLNode.rightEnd(); it3++ ) {
                        SSNode rrNode = getSSNode ( it3->getNodeID() );
                        if (rrNode.getNodeID()==-bestLNode.getNodeID())
                                continue;
                        bool result = rrNode.deleteLeftArc ( bestLNode.getNodeID() );

                        assert ( result );
                }
                for ( ArcIt it2 = bestLNode.leftBegin(); it2 != bestLNode.leftEnd(); it2++ ) {
                        SSNode llNode = getSSNode ( it2->getNodeID() );
                        if (llNode.getNodeID()==-bestLNode.getNodeID())
                                continue;
                        bool result = llNode.deleteRightArc ( bestLNode.getNodeID() );
                        assert ( result );
                }
                if (trueMult[abs( bestLNode.getNodeID())]>0) {
                        // cout<<"you deleted a correct node in guessNodeMultiplicity, node number is:	"<<bestLNode.getNodeID()<<"	cov:"<<bestLNode.getStReadCov()<<"	length:"<<bestLNode.getMarginalLength()<<"	"<< bestResult.first<<"	"<<bestResult.second <<endl;
                }
                bestLNode.deleteAllRightArcs();
                bestLNode.deleteAllLeftArcs();
                bestLNode.invalidate();
                change=true;

        }
        return change;
}

bool DBGraph::deleteExtraLeftLink(SSNode rootNode, int round) {
        bool change=false;
        try {


                ArcIt it = rootNode.leftBegin();
                SSNode bestLNode = getSSNode(it->getNodeID());
                pair<int, pair<double,double> > bestResult= nodesExpMult[abs (bestLNode.getNodeID())];
                it++;
                while(it != rootNode.leftEnd()) {

                        SSNode currNode = getSSNode(it->getNodeID());
                        pair<int, pair<double,double> > curResult= nodesExpMult[abs (currNode.getNodeID())];

                        if (bestResult.second.second<curResult.second.second) {
                                bestLNode=currNode;
                                bestResult=curResult;
                        }
                        it++;

                }

                double ratio=bestResult.second.first/bestResult.second.second;
                //if ( ((bestResult.second.second/rootResult.second.second)>10&&(rootResult.second.first/ bestResult.second.first>10))) {
                if (ratio<.01) {
                        for ( ArcIt it3 = bestLNode.rightBegin(); it3 != bestLNode.rightEnd(); it3++ ) {
                                SSNode rrNode = getSSNode ( it3->getNodeID() );
                                if (rrNode.getNodeID()==-bestLNode.getNodeID())
                                        continue;
                                bool result = rrNode.deleteLeftArc ( bestLNode.getNodeID() );
                                assert ( result );
                        }

                        for ( ArcIt it2 = bestLNode.leftBegin(); it2 != bestLNode.leftEnd(); it2++ ) {
                                SSNode llNode = getSSNode ( it2->getNodeID() );
                                if (llNode.getNodeID()==-bestLNode.getNodeID())
                                        continue;
                                bool result = llNode.deleteRightArc ( bestLNode.getNodeID() );
                                assert ( result );
                        }

                        if (trueMult[abs(bestLNode.getNodeID())]>0) {
                                // cout<<"you deleted a correct node in guessNodeMultiplicity, node number is:	"<<bestLNode.getNodeID()<<"	"<<bestLNode.getStReadCov()<<"	"<<bestResult.first<<"	"<<bestResult.second <<endl;
                        }
                        bestLNode.deleteAllRightArcs();
                        bestLNode.deleteAllLeftArcs();
                        bestLNode.invalidate();
                        change=true;

                }
        }  catch(exception e) {
                cout<<"Fatal error occured in deleteExtraLeftLink around node: "<<rootNode.getNodeID()<< " error Message:"<<e.what()<<endl;

        }

        return change;
}

void DBGraph::setAnchorNodes()
{
        size_t numAnchors = 0, totAnchorLen = 0, totNonAnchorLen = 0, maxAnchors = 0, maxAnchorLen = 0;

        for ( NodeID id = 1; id <= numNodes; id++ ) {
                SSNode node = getSSNode ( id );

                // first criterion: they should have expected multiplicity == 1
                bool isAnchor = ( node.getRoundMult() == 1 ) && ( !node.multIsDubious() );

                // second criterion: left nodes
                size_t numSignLeft = 0;
                for ( ArcIt it = node.leftBegin(); it != node.leftEnd(); it++ ) {
                        SSNode left = getSSNode ( it->getNodeID() );
                        if ( left.getRoundMult() >= 1 ) {
                                numSignLeft++;
                        }
                }
                if ( numSignLeft > 1 ) {
                        isAnchor = false;
                }

                // second criterion: right nodes
                size_t numSignRight = 0;
                for ( ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                        SSNode right = getSSNode ( it->getNodeID() );
                        if ( right.getRoundMult() >= 1 ) {
                                numSignRight++;
                        }
                }
                if ( numSignRight > 1 ) {
                        isAnchor = false;
                }

                node.setAnchor ( isAnchor );

                if ( isAnchor ) {
                        numAnchors++;
                        totAnchorLen += node.getMarginalLength();
                } else {
                        totNonAnchorLen += node.getMarginalLength();
                }

                if ( !isAnchor ) {
                        continue;
                }

                if ( trueMult[id] == 1 ) {
                        maxAnchors++;
                        maxAnchorLen += node.getMarginalLength();
                }

                if ( trueMult[id] == 1 ) {
                        continue;
                }

                cout << "Anchor: " << id << " (length: " << node.getMarginalLength()
                << "): " << node.getNodeKmerCov() << " +/- " << node.getReadStartCov();

                cout << ( trueMult[id] == 1 ? "" : " (ERROR)" ) << endl;
        }

        cout << "Number of anchors based on coverage statistics: " << numAnchors << "/" << maxAnchors << endl;
        cout << "The total lenght of the anchors is: " << totAnchorLen << "/" << maxAnchorLen << endl;
        cout << "The total length of the non-anchors is: " << totNonAnchorLen << endl;
}

bool DBGraph::mergeSingleNodes(bool force)
{
        size_t numDeleted = 0;
        size_t numOfIncorrectConnection=0;
        //comment by mahdi
        //the initial was lID=1
        bool change=false;//  deleteSuspiciousNodes();

        for ( NodeID lID = -numNodes; lID <= numNodes; lID++ ) {
                if ( lID == 0 ) {
                        continue;
                }
                if (lID==-404||lID==473){
                        int stop=0;
                }
                SSNode left = getSSNode ( lID );
                if ( !left.isValid() ) {
                        continue;
                }
                if ( left.getNumRightArcs() != 1 ) {
                        continue;
                }
                //the previous one
                NodeID rID = left.rightBegin()->getNodeID();

                // don't merge palindromic repeats
                if ( rID == -lID ) {
                        continue;
                }
                SSNode right = getSSNode ( rID );
                if ( right.getNumLeftArcs() != 1 ) {
                        continue;
                }
                #ifdef DEBUG
                if ( ( ( trueMult[abs ( lID )] >= 1 ) && ( trueMult[abs ( lID )] == 0 ) ) ||
                        ( ( trueMult[abs ( rID )] >= 1 ) && ( trueMult[abs ( lID )] == 0 ) ) )
                        numOfIncorrectConnection++;
                #endif
                if (!force&& (left.getNodeKmerCov()<cutOffvalue || right.getNodeKmerCov()<cutOffvalue))
                        continue;

                left.deleteRightArc ( rID );
                right.deleteLeftArc ( lID );
                left.inheritRightArcs ( right );
                //left.setExpMult ( left.getExpMult() + right.getExpMult() );
                left.setKmerCov(left.getKmerCov()+right.getKmerCov());
                //comment by mahdi
                left.setReadStartCov(left.getReadStartCov()+right.getReadStartCov());
                right.invalidate();
                numDeleted++;

                deque<SSNode> deq;
                deq.push_back ( left );
                deq.push_back ( right );

                string str;
                convertNodesToString ( deq, str );

                left.setSequence ( str );
                lID--;

        }
        cout << "Concatenated " << numDeleted << " nodes" << endl;
        #ifdef DEBUG
        cout <<numOfIncorrectConnection<< " of connections are between correct and incorrect node"<<endl;
        #endif
        return (numDeleted > 0||change);
}
//comment mahdi creating new procedure for removing diamonds
//////////////////////


bool DBGraph::deleteSuspiciousNodes() {
        bool modify=false;

        if ( !updateCutOffValue(0))
                return false;
        double TP=0;
        double TN=0;
        double FP=0;
        double FN=0;
        for ( NodeID lID = -numNodes; lID <= numNodes; lID++ ) {
                if ( lID == 0 ) {
                        continue;
                }
                SSNode leftNode = getSSNode ( lID );
                if ( !leftNode.isValid() ) {
                        continue;
                }
                if(leftNode.getNumRightArcs()!=1) {
                        continue;
                }
                SSNode rNode = getSSNode ( leftNode.rightBegin()->getNodeID() );
                rNode.getLeftKmer();
                if(rNode.getNumLeftArcs()!=1 ) {
                        continue;
                }
                //double leftCov=leftNode.getExpMult()/leftNode.getMarginalLength();
                //double rightCov=rNode.getExpMult()/rNode.getMarginalLength();
                double leftCov=leftNode.getNodeKmerCov();
                double rightCov=rNode.getNodeKmerCov();

                if (leftCov < rightCov) {

                        if (leftCov<this->safeValueCov) {
                                if (trueMult[abs( leftNode.getNodeID())]>0) {
                                        TP++;
                                }
                                else {
                                        FP++;
                                }

                                removeNode(leftNode);
                                modify=true;
                        } else {
                                if (trueMult[abs (leftNode.getNodeID())]>0) {
                                        TN++;
                                }
                                else {
                                        FN++;
                                }
                        }
                }
                else {

                        if (rightCov< this->safeValueCov) {
                                if (trueMult[abs (rNode.getNodeID())]>0) {
                                        TP++;
                                }
                                else {
                                        FP++;
                                }
                                removeNode(rNode);
                                modify=true;
                        } else {
                                if (trueMult[abs(rNode.getNodeID())]>0) {
                                        TN++;
                                }
                                else {
                                        FN++;
                                }
                        }


                }

        }

        cout<<"TP:"<< TP<<"		TN:"<< TN<<"	FP:" <<FP<<"	FN:"<< FN<<endl;
        return modify;
}


size_t DBGraph::getLowestArcMultiplicity ( NodeID left, NodeID right )
{
        SSNode node = getSSNode ( right );
        int lowArcMult = node.getLoExpMult();

        for ( ArcIt it = node.leftBegin(); it != node.leftEnd(); it++ ) {
                if ( it->getNodeID() == left ) {
                        continue;
                }

                lowArcMult -= getSSNode ( it->getNodeID() ).getHiExpMult();
        }

        return ( lowArcMult > 0 ) ? lowArcMult : 0;
}

bool DBGraph::arcFilteringBasedOnExpMultiplicity()
{
        size_t numArcsDetatched = 0, incorrDetach = 0;

        for ( NodeID id = -numNodes; id <= numNodes; id++ ) {
                if ( id == 0 ) {
                        continue;
                }

                SSNode node = getSSNode ( id );
                if ( !node.isValid() ) {
                        continue;
                }

                // at this point, there is too much flow to the right
                vector<NodeID> toDelete;
                for ( ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                        SSNode right = getSSNode ( it->getNodeID() );
                        if ( right.getRoundMult() > 0 ) {
                                continue;
                        }

                        toDelete.push_back ( it->getNodeID() );
                }

                for ( size_t i = 0; i < toDelete.size(); i++ ) {
                        node.deleteRightArc ( toDelete[i] );
                        SSNode toDel = getSSNode ( toDelete[i] );
                        toDel.deleteLeftArc ( id );
                        numArcsDetatched++;
                        if ( getSSNode ( toDelete[i] ).getFlag() == 1 ) {
                                cout << "Incorrectly removed arcs between node " << id << " (" << node.getMarginalLength() << "): " << node.getNodeKmerCov() << " +/- " << node.getReadStartCov() << " and node " << toDelete[i] << " (" << toDel.getMarginalLength() << "): " << toDel.getNodeKmerCov() << " +/- " << toDel.getReadStartCov() << endl;
                                incorrDetach++;
                        }
                }
        }

        cout << "Number of arcs detached " << numArcsDetatched << " arcs" << endl;
        cout << "Number of incorrectly detached " << incorrDetach << " arcs" << endl;

        return numArcsDetatched > 0;
}

bool DBGraph::arcFilteringBasedOnSignExpMultiplicity()
{
        size_t numArcsDetatched = 0, incorrDetach = 0;

        for ( NodeID id = -numNodes; id <= numNodes; id++ ) {
                if ( id == 0 ) {
                        continue;
                }

                SSNode node = getSSNode ( id );
                if ( !node.isValid() ) {
                        continue;
                }

                // at this point, there is too much flow to the right
                vector<NodeID> toDelete;
                for ( ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                        SSNode right = getSSNode ( it->getNodeID() );
                        if ( right.multIsDubious() ) {
                                continue;
                        }
                        if ( right.getRoundMult() > 0 ) {
                                continue;
                        }

                        toDelete.push_back ( it->getNodeID() );
                }

                for ( size_t i = 0; i < toDelete.size(); i++ ) {
                        node.deleteRightArc ( toDelete[i] );
                        SSNode toDel = getSSNode ( toDelete[i] );
                        toDel.deleteLeftArc ( id );
                        numArcsDetatched++;
                        if ( toDel.getFlag() == 1 ) {
                                cout << "Incorrectly removed arcs between node " << id << " (" << node.getMarginalLength() << "): " << node.getNodeKmerCov() << " +/- " << node.getReadStartCov() << " and node " << toDelete[i] << " (" << toDel.getMarginalLength() << "): " << toDel.getNodeKmerCov() << " +/- " << toDel.getReadStartCov() << endl;
                                incorrDetach++;
                        }
                }
        }

        cout << "Number of arcs detached " << numArcsDetatched << " arcs" << endl;
        cout << "Number of incorrectly detached " << incorrDetach << " arcs" << endl;

        return numArcsDetatched > 0;
}

bool DBGraph::flowConservationBasedArcFiltering()
{
        size_t numArcsDetatched = 0, incorrDetach = 0;

        for ( NodeID id = -numNodes; id <= numNodes; id++ ) {
                if ( id == 0 ) {
                        continue;
                }

                SSNode node = getSSNode ( id );
                if ( !node.isValid() ) {
                        continue;
                }

                // only one arc is present, continue
                if ( node.getNumRightArcs() == 1 ) {
                        continue;
                }

                // check if there are potential candidates to delete
                bool OK = true;
                for ( ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ )
                        if ( getSSNode ( it->getNodeID() ).getRoundMult() == 0 ) {
                                OK = false;
                        }
                        if ( OK ) {
                                continue;
                        }

                        int leftMult = node.getHiExpMult();  // left flow

                        int rightMult = 0;                   // right flow
                        for ( ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                                int thisRight = getLowestArcMultiplicity ( id, it->getNodeID() );
                                rightMult += ( thisRight > 0 ) ? thisRight : 1;
                        }

                        // if there is no excess on right flow, continue
                        if ( leftMult >= rightMult ) {
                                continue;
                        }

                        // at this point, there is too much flow to the right
                        multimap<double, NodeID> toDelete;
                        for ( ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                                SSNode right = getSSNode ( it->getNodeID() );
                                if ( right.getRoundMult() > 0 ) {
                                        continue;
                                }

                                toDelete.insert ( pair<double, NodeID> ( right.getNodeKmerCov(), it->getNodeID() ) );
                        }

                        int maxDelete = ( rightMult -leftMult );
                        for ( map<double, NodeID>::iterator it = toDelete.begin(); it != toDelete.end(); it++ ) {
                                if ( maxDelete == 0 ) {
                                        break;
                                }
                                SSNode nodeToDelete = getSSNode ( it->second );
                                node.deleteRightArc ( it->second );
                                nodeToDelete.deleteLeftArc ( id );
                                numArcsDetatched++;
                                maxDelete--;

                                if ( getDSNode ( abs ( it->second ) ).getFlag() == 1 ) {
                                        incorrDetach++;
                                        cout << " Wrong !!!!!!! " << endl;
                                        cout << nodeToDelete.getNodeKmerCov() << " +/- " << nodeToDelete.getReadStartCov() << endl;
                                } else {
                                        // cout << " Correct" << endl;
                                }
                        }

                        if ( maxDelete > 0 ) {
                                cout << " left arcs remaining " << endl;
                        }
        }

        //    cout << "Flow conserved on " << numOK << " nodes" << endl;
        //   cout << "Flow NOT conserved on " << numNOK << " nodes" << endl;
        cout << "Number of arcs detached " << numArcsDetatched << " arcs" << endl;
        cout << "Number of incorrectly detached " << incorrDetach << " arcs" << endl;

        return numArcsDetatched > 0;
}

bool DBGraph::findIsolatedCluster ( SSNode& curr, vector<NodeID>& cluster )
{
        curr.setVisited ( true );
        cluster.push_back ( curr.getNodeID() );

        for ( ArcIt it = curr.leftBegin(); it != curr.leftEnd(); it++ ) {
                NodeID nextID = it->getNodeID();
                SSNode next = getSSNode ( nextID );

                if ( !next.isValid() || next.isVisited() || next.getRoundMult() > 0 ) {
                        return false;
                }

                if ( !findIsolatedCluster ( next, cluster ) ) {
                        return false;
                }
        }

        for ( ArcIt it = curr.rightBegin(); it != curr.rightEnd(); it++ ) {
                NodeID nextID = it->getNodeID();
                SSNode next = getSSNode ( nextID );

                if ( !next.isValid() || next.isVisited() || next.getRoundMult() > 0 ) {
                        return false;
                }

                if ( !findIsolatedCluster ( next, cluster ) ) {
                        return false;
                }
        }

        return true;
}

void DBGraph::filterIsolatedChaff()
{
        size_t numDeleted = 0, incorrDeleted = 0;

        for ( NodeID id = 1; id <= numNodes; id++ ) {
                SSNode node = getSSNode ( id );

                if ( !node.isValid() ) {
                        continue;
                }

                if ( node.getRoundMult() > 0 ) {
                        continue;
                }

                vector<NodeID> cluster;
                bool foundCluster = findIsolatedCluster ( node, cluster );

                if ( !foundCluster ) {
                        continue;
                }

                for ( size_t i = 0; i < cluster.size(); i++ ) {
                        if ( getSSNode ( cluster[i] ).getRoundMult() > 0 ) {
                                cout << "ERROR " << endl;
                        }
                        getSSNode ( cluster[i] ).invalidate();
                        if ( getDSNode ( abs ( cluster[i] ) ).getFlag() == 1 ) {
                                incorrDeleted++;
                        }
                }

                numDeleted += cluster.size();
        }

        for ( NodeID id = 1; id <= numNodes; id++ ) {
                getSSNode ( id ).setVisited ( false );
        }

        cout << "Deleted " << numDeleted << " isolated nodes with expected coverage 0" << endl;
        cout << "Incorrectly deleted " << incorrDeleted << " nodes." << endl;
}

void DBGraph::validateKmerFrequency()
{
        /*   table = new KmerNodeTable(settings, numNodes);
         *       table->populateTable(nodes);
         *
         *       vector<size_t> oldFreq = kmerFreq;
         *
         *       kmerFreq.clear();
         *       kmerFreq.resize(numNodes+1, 0);
         *
         *       cout << "Estimating the node coverage: " << endl;
         *       for (size_t i = 0; i < settings.getInputCount(); i++) {
         *               const Input &input = settings.getInput(i);
         *
         *               cout << "Processing file " << i+1 << "/"
         *                    << settings.getInputCount() << ": "
         *                    << input.getFilename()
         *                    << ", type: " << input.getFileType()
         *                    << ", reads: " << input.getReadType() << endl;
         *
         *               countKmerFreq(input.getIsolatedFilename(),
         *                             input.getNumIsolatedReads(),
         *table);
         *               countKmerFreq(input.getPairedFilename(),
         *                             input.getNumPairedReads(),
         *table);
}

delete table;

for (size_t i = 0; i < kmerFreq.size(); i++) {
        if (kmerFreq[i] != oldFreq[i])
                cout << "ERROR IN KMER FREQUENCY !!" << endl;
}*/
}


// ============================================================================
// PRIVATE NODE COVERAGE ROUTINES
// ============================================================================

bool sortByLength ( const NodeID& left, const NodeID& right )
{
        return DBGraph::graph->getDSNode ( left ).getMarginalLength() >
        DBGraph::graph->getDSNode ( right ).getMarginalLength();
}

double DBGraph::getInitialEstimateForCoverage ( const ReadLibrary& input,
        vector<size_t> &readFreq ) const
{
        // sort the nodes from big to small + compute the total nodes length
        vector<NodeID> sortedNodes;
        sortedNodes.reserve ( numNodes );
        size_t totalSize = 0;
        for ( NodeID id = 1; id <= numNodes; id++ ) {
                sortedNodes.push_back ( id );
                totalSize += getDSNode ( id ).getMarginalLength();
        }

        DBGraph::graph = this;
        sort ( sortedNodes.begin(), sortedNodes.end(), sortByLength );

        size_t RL = input.getReadLength();
        size_t k = Kmer::getK();

        // for the biggest nodes compute the coverage
        size_t currSize = 0, numReads = 0;
        for ( size_t i = 0; i < sortedNodes.size(); i++ ) {
                size_t MNL = getDSNode ( sortedNodes[i] ).getMarginalLength();
                currSize += MNL + RL - k;
                numReads += readFreq[sortedNodes[i]];
                if ( currSize > 0.15*totalSize ) {
                        break;
                }
        }

        return double ( numReads ) / double ( currSize );
}

bool DBGraph::setNodeMultiplicity ( const vector<vector<size_t> > &readFreq )
{
        /*  bool changedST = false;
         *    //ofstream ofs ("nodeCov.txt"); // <------ Writing file
         *
         *    const size_t k = Kmer::getWordSize();
         *
         *    for ( NodeID id = 1; id <= numNodes; id++ ) {
         *        DSNode &node = getDSNode ( id );
         *        const size_t MNL = node.getMarginalLength();
         *
         *        // count the number of observed reads for this node
         *        size_t obsReads = 0;
         *        for ( size_t i = 0; i < settings.getInputCount(); i++ ) {
         *            obsReads += readFreq[i][id];
}

// determine number of expected reads for this node
double expReads = 0;
for ( size_t i = 0; i < settings.getInputCount(); i++ ) {
        size_t RL = settings.getInput ( i ).getReadLength();
        double readStartCov = settings.getInput ( i ).getReadStartCov();
        expReads += double ( MNL + RL - k ) * readStartCov;
}

// calculate the new multiplicity
double obsMult = double ( obsReads ) / expReads;
int rndMult = round ( obsMult ); // important to make this "signed"

// nominator = chance to observe "obsReads" given this multiplicity is correct
double lambdaN = double ( rndMult ) * expReads;
double nominator = gsl_ran_poisson_pdf ( obsReads, lambdaN );

// denom = chance to observe "obsReads" given this multiplicity is NOT correct
double denom = 0.0;
for ( int i = rndMult + 1; i < 2 * rndMult + 10; i++ ) {
        double lambdaD = double ( i ) * expReads;
        double newTerm = gsl_ran_poisson_pdf ( obsReads, lambdaD );
        denom += newTerm;

        // do not compute beyond 5 digits of accuracy
        if ( newTerm / denom < 1e-5 ) {
                break;
}
}

for ( int i = rndMult - 1; i > 0; i-- ) { // make sure rndMult is of "int" type
        double lambdaD = double ( i ) * expReads;
        double newTerm = gsl_ran_poisson_pdf ( obsReads, lambdaD );
        denom += newTerm;

        // do not compute beyond 5 digits of accuracy
        if ( newTerm / denom < 1e-5 ) {
                break;
}
}

// how much more likely is it that the rndMult is the actual multiplicity?
double odds = nominator / denom;

// + infinite occurs quite often when denom == 0
if ( std::isinf ( odds ) || odds > 1E10 ) {
        odds = 1E10;
}

// could be the result of 0/0 -> do not trust multiplicity
if ( std::isnan ( odds ) ) {
        odds = 0.0;
}

//cout << rndMult << "\t" << odds << "\t"
//     << nominator << "\t" << denom << endl;

//ofs << id << "\t" << MNL << "\t" << expReads << "\t" << obsReads << "\t"
//    << trueMult[id] << "\t" << rndMult << "\t" << odds << endl;

size_t oldMult = node.getRoundMult();

node.setExpMult ( obsMult );

//comment by mahdi
//node.setMultStd ( odds );

if ( node.getRoundMult() != oldMult ) {
        changedST = true;
}
}

//ofs.close();

return changedST;*/
}

double DBGraph::estimateReadStartCoverage ( const ReadLibrary &input,
        const vector<size_t> &readFreq ) const
{
        size_t nom = 0, denom = 0;
        size_t k = Kmer::getK();
        size_t RL = input.getReadLength();

        for ( NodeID id = 1; id <= numNodes; id++ ) {
                DSNode& node = getDSNode ( id );

                if ( !node.isValid() ) {
                        continue;
                }
                if ( node.getRoundMult() == 0 ) {
                        continue;
                }

                size_t MNL = node.getMarginalLength();

                nom += readFreq[id];
                denom += ( MNL + RL - k )  * node.getRoundMult();
        }

        return ( double ) nom / ( double ) denom;
}


