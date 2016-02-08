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
#include "settings.h"
#include "library.h"

using namespace std;

void DBGraph::parseReads(size_t thisThread,
                         vector<string>& readBuffer)
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

        size_t blockID, recordOffset;
        while (inputs->getReadChunk(myReadBuf, blockID, recordOffset))
                parseReads(thisThread, myReadBuf);
}

void DBGraph::countNodeandArcFrequency(LibraryContainer &inputs)
{
        const unsigned int& numThreads = settings.getNumThreads();

        cout << "Populating table... ";
        cout.flush();
        populateTable();
        cout << "done" << endl;

        cout << "Number of threads: " << numThreads << endl;
        cout << "Counting node and arc frequency: " << endl;

        inputs.startIOThreads(settings.getThreadWorkSize(),
                              settings.getThreadWorkSize() * settings.getNumThreads());

        // start worker threads
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&DBGraph::workerThread, this, i, &inputs);

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        inputs.joinIOThreads();

        depopulateTable();
}

// ============================================================================
// PRIVATE NODE COVERAGE / MULTIPLICITY
// ============================================================================

//comment by mahdi


bool cmssn(const SSNode& first, const SSNode & second ) {

    return first.getMarginalLength()>second.getMarginalLength();


}
/**
 * this routine do the following things
 * 1. calculation of avg and std of node kmer coverage
 * 2. calculation of read start coverage mean
 * 3. for every node it calculates the MULTIPLICITY
 * 4. for every node it calculates the certainity of our guess about MULTIPLICITY
 * 5. for every node it calculates the inCorrctnessRatio of our guess about MULTIPLICITY
 *
 */


void DBGraph::extractStatistic(int round) {
        vector <SSNode> nodeArray;
        int percentage=5;
        double sumOfReadStcov=0;
        size_t totalLength=0;
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

        double sizeLimit=0;
        sizeLimit= (totalLength*percentage)/100;
        /*if (settings.getGenomeSize()>0)
                sizeLimit=((totalLength*percentage)/100)>settings.getGenomeSize()?settings.getGenomeSize():((totalLength*percentage)/100);
        else
                sizeLimit= (totalLength*percentage)/100;*/
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
                double maxProb=0;
                double newValue=avg*nodeMultiplicity*node.getMarginalLength();
                double newProbability=Util::poissonPDF(node.getReadStartCov(),newValue);
                while(newProbability>maxProb) {
                        nodeMultiplicity++;
                        maxProb=newProbability;
                        newValue=avg*nodeMultiplicity*node.getMarginalLength();
                        newProbability=Util::poissonPDF(node.getReadStartCov(),newValue);
                }

                nodeMultiplicity--;

                double denominator=0;
                double currentProb=maxProb;

                int i = 1;
                bool minus=true;
                do{
                        currentProb=newProbability;
                        double newValue=0;
                        if (minus&& nodeMultiplicity>i){
                                newValue=avg*(nodeMultiplicity-i)*node.getMarginalLength();
                                minus=false;
                        }
                        else{
                                newValue=avg*(nodeMultiplicity+i)*node.getMarginalLength();
                                minus=true;
                                i++;
                        }
                        newProbability=Util::poissonPDF(node.getReadStartCov(),newValue);
                        denominator=denominator+newProbability;
                }while(abs(newProbability-currentProb)> .000001|| i<5);
                double confidenceRatio=maxProb/denominator;

                double expectToSee=avg*nodeMultiplicity*node.getMarginalLength();
                double observedprob=0;
                if(node.getReadStartCov()<expectToSee)
                        observedprob=Util::poissonPDF(node.getReadStartCov(),expectToSee);
                else
                        observedprob=Util::poissonPDF(expectToSee,expectToSee);
                double inCorrctnessRatio=Util::poissonPDF(expectToSee,expectToSee)/observedprob;
                SSNode reverseNode=getSSNode(-node.getNodeID());
                reverseNode.setExpMult(nodeMultiplicity);
                node.setExpMult(nodeMultiplicity);
                nodesExpMult[abs(node.getNodeID())]=make_pair(nodeMultiplicity,make_pair( confidenceRatio,inCorrctnessRatio));


        }

}
/**
 * this routine is used to check if the currNode can be a parallel path to some other nodes
 * the root node is also given.
 * this routine is used as a subroutin of deleteUnreliableNodes
 * @return true if the given nodes is a bubble initiate from rootNode
 *
 */

bool DBGraph::nodeIsBubble(SSNode rootNode, SSNode currNode){
        if(rootNode.getNumRightArcs()>1&& hasLowCovNode(rootNode)){
                vector<pair<vector<NodeID>, vector<NodeID>> > possibleBubbles=searchForParallelNodes(rootNode, 500);
                for (auto b : possibleBubbles){
                        vector<NodeID> upPath=b.first;
                        vector<NodeID> downPath=b.second;
                        SSNode up=getSSNode(b.first[1]);
                        SSNode down=getSSNode(b.second[1]);
                        if (up.getNodeID()!=currNode.getNodeID()&& down.getNodeID()!=currNode.getNodeID())
                                continue;
                        if(up.isValid()&&down.isValid())
                        {
                                bool upIsBubble=true;
                                if (whichOneIsbubble(rootNode,upIsBubble,up, down,false,this->estimatedKmerCoverage)){
                                        if (upIsBubble&& up.getNodeID()==currNode.getNodeID())
                                                return true;
                                        if(!upIsBubble&&down.getNodeID()==currNode.getNodeID())
                                                return true;
                                }

                        }


                }
        }
        return false;
}
/**
 * check nodes to see wether is reliable or not
 * a node is reliable if it has a high confidenceRatio and low inCorrctnessRatio
 * these terms are defined in extractStatistic routine
 * @return true if the node is reliable
 */


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
/**
 * this routine loops over all nodes, for a reliable node N with MULTIPLICITY(M)
 * it should have at most M outgoing arcs. Therefore we find extra nodes with low coverage
 * or extra nodes wich are appeared as tips or bubbles.
 * @return true if any changes happens
 * the second part of this routine looks for the adjacent nodes with same MULTIPLICITY
 * if they have some outgoing or ingoing arcs with low coverage those nodes should be
 * deleted and these two adjacent reliable nodes should be connected later.
 */
bool DBGraph::deleteUnreliableNodes(){
      bool changeIn1=deleteExtraAttachedNodes();
      bool changeIn2=connectSameMulNodes();
      return (changeIn1||changeIn2);
}
bool DBGraph::deleteExtraAttachedNodes(){
        double tp=0, tn=0, fp=0,fn=0;
        size_t numOfDel=0;
        for ( NodeID lID =-numNodes; lID <= numNodes; lID++ ) {
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
                if(node.getExpMult()<node.getNumRightArcs())
                        continue;

                ArcIt it = node.rightBegin();
                SSNode currNode = getSSNode(it->getNodeID());
                while(it != node.rightEnd()) {
                        bool tip=currNode.getNumRightArcs()==0&&currNode.getNumLeftArcs()==1;
                        SSNode firstNodeInUpPath;
                        SSNode firstNodeInDoPath;
                        bool bubble=nodeIsBubble(node,currNode);
                        double threshold=(!tip&& !bubble)? this->redLineValueCov:(this->estimatedKmerCoverage-this->estimatedMKmerCoverageSTD*2>this->redLineValueCov)?this->estimatedKmerCoverage-this->estimatedMKmerCoverageSTD*2:this->redLineValueCov;
                        if (currNode.getNodeKmerCov()<threshold &&removeNode(currNode)){
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
        cout<< "TP:     "<<tp<<"        TN:     "<<tn<<"        FP:     "<<fp<<"        FN:     "<<fn<<endl;
        cout << "Sensitivity: ("<<100*((double)tp/(double)(tp+fn))<<"%)"<<endl;
        cout<<"Specificity: ("<<100*((double)tn/(double)(tn+fp))<<"%)"<<endl;
        #endif
        return (numOfDel>0);
}
bool DBGraph::connectSameMulNodes(){
        size_t secondFP=0;
        size_t secondTP=0;
        size_t numOfDel=0;
        double threshold=this->redLineValueCov;
        for ( NodeID lID =-numNodes; lID <= numNodes; lID++ ) {
                if (lID==0)
                        continue;
                SSNode node = getSSNode ( lID );
                if(!node.isValid())
                        continue;
                if (!checkNodeIsReliable(node))
                        continue;
                ArcIt it = node.rightBegin();
                SSNode currNode = getSSNode(it->getNodeID());
                bool found=false;
                while(it != node.rightEnd()) {
                        if (!checkNodeIsReliable(currNode))
                                it++;
                        if (currNode.getExpMult()!=node.getExpMult() ||currNode.getExpMult()<2)
                                it++;
                        found=true;
                        break;
                }
                if (found)
                {
                        ArcIt it = node.rightBegin();
                        while(it != node.rightEnd()) {
                                SSNode victim=getSSNode(it->getNodeID());
                                if(victim.getNodeID()!=currNode.getNodeID()&&victim.getNodeKmerCov()<threshold&& removeNode(victim)){
                                        numOfDel++;
                                        #ifdef DEBUG
                                        if (trueMult[abs(victim.getNodeID())]>0)
                                                secondFP++;
                                        else
                                                secondTP++;
                                        #endif
                                        break;
                                }
                                it++;
                        }
                        it = currNode.leftBegin();
                        while(it!=currNode.leftEnd()){
                                SSNode victim=getSSNode(it->getNodeID());
                                if (!victim.isValid())
                                        it++;
                                if(victim.getNodeID()!=node.getNodeID()&&victim.getNodeKmerCov()<threshold &&removeNode(victim)){
                                        numOfDel++;
                                        #ifdef DEBUG
                                        if (trueMult[abs(victim.getNodeID())]>0)
                                                secondFP++;
                                        else
                                                secondTP++;
                                        #endif
                                        break;
                                }
                                it++;
                        }
                }

        }

        if (numOfDel>0)
                cout << "Number of deleted nodes in connectSameMulNodes: " << numOfDel << endl;
        #ifdef DEBUG
        cout<<endl<<"second FP: "<<secondFP<<endl;
        cout<<"second TP: "<<secondTP<<endl;
        #endif
        return (numOfDel>0);
}



/**
 * merge single nodes together if they have only one ingoing and only one outgoin Arc
 * node kmer coverage and start read coverage should be updated after deleting the node
 * @return true if they merge any nodes.
 *
 */
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

                if (!force&& (left.getNodeKmerCov()<cutOffvalue || right.getNodeKmerCov()<cutOffvalue))
                        continue;

                #ifdef DEBUG
                if (trueMult.size()>0)
                if ( ( ( trueMult[abs ( lID )] >= 1 ) && ( trueMult[abs ( rID )] == 0 ) ) ||
                        ( ( trueMult[abs ( rID )] >= 1 ) && ( trueMult[abs ( lID )] == 0 ) ) ){
                                numOfIncorrectConnection++;
                                trueMult[abs(lID)]=0;
                                trueMult[abs(rID)]=0;

                        }

                #endif

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
        if (numDeleted>0)
        cout << "Concatenated " << numDeleted << " nodes" << endl;
        #ifdef DEBUG
        cout <<numOfIncorrectConnection<< " of connections are between correct and incorrect node"<<endl;
        #endif
        return (numDeleted > 0||change);
}


/**
 * this routine is used to before merging nodes, if two adjacent nodes have a very different coverage
 * they shouldn't be merged together
 * @return true if they do any delete or modification
 *
 * this routine is not used currently
 */
bool DBGraph::deleteSuspiciousNodes() {
        bool modify=false;
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
                                #ifdef DEBUG
                                if (trueMult[abs( leftNode.getNodeID())]>0) {
                                        TP++;
                                }
                                else {
                                        FP++;
                                }
                                #endif
                                removeNode(leftNode);
                                modify=true;
                        } else {
                                #ifdef DEBUG
                                if (trueMult[abs (leftNode.getNodeID())]>0) {
                                        TN++;
                                }
                                else {
                                        FN++;
                                }
                                #endif
                        }
                }
                else {

                        if (rightCov< this->safeValueCov) {
                                #ifdef DEBUG
                                if (trueMult[abs (rNode.getNodeID())]>0) {
                                        TP++;
                                }
                                else {
                                        FP++;
                                }
                                #endif
                                removeNode(rNode);
                                modify=true;
                        } else {
                                #ifdef DEBUG
                                if (trueMult[abs(rNode.getNodeID())]>0) {
                                        TN++;
                                }
                                else {
                                        FN++;
                                }
                                #endif
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
/**
 * remvoes nodes with coverage lower that input argument
 * @return true if some nodes will be deleted
 */
bool DBGraph::filterCoverage(float cutOff)
{
        cout << endl << " ================== Coverage Filter ==================" << endl;

        cout << "Cut-off value for removing nodes: " << cutOff << endl;
        int tp=0;
        int tn=0;
        int fp=0;
        int fn=0;
        size_t numFiltered = 0;
        for (NodeID i =1; i <= numNodes; i++) {

                SSNode n = getSSNode(i);
                if(!n.isValid())
                        continue;
                if (n.getNodeKmerCov()>cutOff) {
                        #ifdef DEBUG
                        if (trueMult.size()>0&& trueMult[abs(i)] >= 1)
                                tn++;
                        else
                                fn++;
                        #endif
                        continue;
                }
                bool Del=removeNode(n);
                if (Del){
                        numFiltered++;
                        #ifdef DEBUG
                        if (trueMult.size()>0&& trueMult[abs(i)] >= 1)
                                fp++;
                        else
                                tp++;
                        #endif
                }
                else
                {
                        #ifdef DEBUG
                        if (trueMult[abs(i)] >= 1)
                                tn++;
                        else
                                fn++;
                        #endif
                }
        }
        cout << "Number of nodes deleted based on coverage: " << numFiltered<<endl;
        #ifdef DEBUG
        cout << " Gain value is ("<<100*((double)(tp-fp)/(double)(tp+fn))<< "%)"<<endl;
        cout<< "TP:     "<<tp<<"        TN:     "<<tn<<"        FP:     "<<fp<<"        FN:     "<<fn<<endl;
        cout << "Sensitivity: ("<<100*((double)tp/(double)(tp+fn))<<"%)"<<endl;
        cout<<"Specificity: ("<<100*((double)tn/(double)(tn+fp))<<"%)"<<endl;
        #endif

        return numFiltered > 0;
}

/**
 * all routines call this routine to remove nodes
 * removes the node if its size be shorter than max_node_size_ToDelete
 * @return true if node will be deleted,
 *
 */

bool DBGraph::removeNode(SSNode & rootNode) {
        if (rootNode.getMarginalLength()>maxNodeSizeToDel)
                return false;
        for ( ArcIt it2 = rootNode.leftBegin(); it2 != rootNode.leftEnd(); it2++ ) {
                SSNode llNode = getSSNode ( it2->getNodeID() );
                if (llNode.getNodeID()==-rootNode.getNodeID())
                        continue;
                bool result = llNode.deleteRightArc ( rootNode.getNodeID() );
                assert ( result );
        }
        for ( ArcIt it3 = rootNode.rightBegin(); it3 != rootNode.rightEnd(); it3++ ) {
                SSNode rrNode = getSSNode ( it3->getNodeID() );
                if (rrNode.getNodeID()==-rootNode.getNodeID())
                        continue;
                bool result = rrNode.deleteLeftArc ( rootNode.getNodeID() );
                assert ( result );
        }
        rootNode.deleteAllRightArcs();
        rootNode.deleteAllLeftArcs();
        rootNode.invalidate();
        return true;
}
