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
#include "tkmer.h"
#include "kmeroverlap.h"
#include "readfile/sequencefile.h"
#include "util.h"
#include "settings.h"
#include "nodeendstable.h"
#include "alignment.h"
#include <string>
#include <fstream>

#include <gsl/gsl_math.h>
#include "ExpMaxClustering.h"

#include "kmernode.h"

using namespace std;

DSNode* SSNode::nodes = NULL;
const DBGraph* DBGraph::graph = NULL;


DBGraph::DBGraph(const Settings& settings) : settings(settings), nodes(NULL),
    arcs(NULL), numNodes(0), numArcs(0), table(NULL),
    mapType(SHORT_MAP) {
    DBGraph::graph = this;
    //mahdi comment my
    initialize();
}




void DBGraph::initialize()
{
    //correct it later it should read this information from setting file, but setting dosn't have such infomratin right now
    readLength=100;//settings.getReadLength();
    coverage=100;//settings.getCoverage();
    minCertainVlueCov=2;
    minSafeValueCov=2.5;
    minRedLineValueCov=3;
    kmerSize= settings.getK();
    estimatedKmerCoverage=(coverage/readLength)*(readLength-kmerSize+1)/12;
    estimatedMKmerCoverageSTD=sqrt( estimatedKmerCoverage);
    certainVlueCov=minCertainVlueCov;
    safeValueCov=minSafeValueCov;
    redLineValueCov=minRedLineValueCov;
    updateCutOffValueRound=1;
    maxNodeSizeToDel=readLength*3;
}

int DBGraph::parseLine(char* line)
{
    int i = strlen(line);
    while (*line < '0' || *line > '9') line++;
    line[i-3] = '\0';
    i = atoi(line);
    return i;
}

double DBGraph::getMemoryUsedByProc()
{
    double res = 0;


    FILE* proc_self_status_v = fopen("/proc/self/status", "r");
    char line[128];

    while (fgets(line, 128, proc_self_status_v) != NULL)
    {
        if (strncmp(line, "VmSize:", 7) == 0)
        {
            res = parseLine(line);
            break;
        }
    }
    res = res/1024;
    fclose(proc_self_status_v);
    return res;
}


void DBGraph::canonicalBubble()
{
    // sanity check
    for (NodeID i = -numNodes; i <= numNodes; i++) {
        if (i == 0)     // skip node 0, doesn't exist
            continue;
        SSNode n = getSSNode(i);
        if(!n.isValid())
            continue;
        // cout << (int)n.getFlag() << endl;
        if (n.getNumRightArcs() != 2)
            continue;
        ArcIt it = n.rightBegin();
        SSNode rn1 = getSSNode(it->getNodeID());
        it++;
        SSNode rn2 = getSSNode(it->getNodeID());
        if (rn1.getNumLeftArcs() != 1)
            continue;
        if (rn1.getNumRightArcs() != 1)
            continue;
        if (rn2.getNumLeftArcs() != 1)
            continue;
        if (rn2.getNumRightArcs() != 1)
            continue;
        if (rn1.rightBegin()->getNodeID() != rn2.rightBegin()->getNodeID())
            continue;
        // SSNode rr2 = getSSNode(rn1.rightBegin()->getNodeID());
        //  cout << (int)n.getFlag() << "\t" << (int)rn1.getFlag() << "\t" << (int)rn2.getFlag() << "\t" << (int)rr2.getFlag() << endl;
        /*if (rn1.getFlag() == 0 || rn2.getFlag() == 0)
                cout << "OK" << endl;
        else
                cout << "Ouch" << endl;*/
        //   cout << "Canonical node " << i << endl;
    }
}



bool DBGraph::continueEdit(size_t& bigestN50, string& nodeFileName, string& arcFileName,string& metaDataFilename) {
    size_t genomeSize=settings.getGenomeSize();
    size_t currentN50=getN50();
    if ((sizeOfGraph/(genomeSize)<.9) && (genomeSize/1000)<bigestN50) { // && preSize<(genomeSize+genomeSize*.1)|| (preSize<genomeSize)
        if (currentN50>bigestN50) {
            bigestN50=currentN50;
        }
        return false;
    }
    if (currentN50>=bigestN50) {
        if(currentN50>bigestN50)
            writeGraph(nodeFileName, arcFileName, metaDataFilename);
        bigestN50=currentN50;
    }
    return true;
}
bool DBGraph::filterChimeric(int round)
{
        cout<<"*********************<<filter Chimeric starts>>......................................... "<<endl;
        size_t numFiltered = 0;
        size_t numOfIncorrectlyFiltered=0;
        // sanity check
        for (NodeID i = -numNodes; i <= numNodes; i++) {
                if (i==0)
                        continue;
                SSNode n = getSSNode(i);
                if(!n.isValid())
                        continue;

                if (n.getNumRightArcs() != 1)
                        continue;
                if (n.getNumLeftArcs() != 1)
                        continue;
                if (n.getMarginalLength() < settings.getK())  //????
                        continue;
                if (n.getNodeKmerCov() > this->redLineValueCov || n.getMarginalLength()>this->maxNodeSizeToDel )
                        continue;
                SSNode rrNode = getSSNode(n.rightBegin()->getNodeID());
                if (rrNode.getNodeID()==-n.getNodeID())
                        continue;

                // first, disconnect the right node from all its right neighbors

                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {

                        SSNode rrNode = getSSNode(it->getNodeID());
                        if (rrNode.getNodeID()==-n.getNodeID())
                                continue;
                        bool result = rrNode.deleteLeftArc(n.getNodeID());
                        assert(result);
                        //break;
                }

                // first, disconnect the right node from all its right neighbors
                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                        SSNode llNode = getSSNode(it->getNodeID());
                        if (llNode.getNodeID()==-n.getNodeID())
                                continue;
                        bool result = llNode.deleteRightArc(n.getNodeID());
                        assert(result);
                        //break;
                }

                n.deleteAllLeftArcs();
                n.deleteAllRightArcs();
                n.invalidate();

                numFiltered++;
                #ifdef DEBUG
                if (trueMult[i] >= 1)
                        numOfIncorrectlyFiltered++;

                #endif

        }
        cout << "\tIncorrectly deleted chimeric connection, namely node: " << numOfIncorrectlyFiltered << endl;
        cout << "Number of chimeric nodes deleted: " << numFiltered << endl;
        return numFiltered > 0;
}

bool DBGraph::updateCutOffValue(int round)
{
        char roundStr[15];
        sprintf(roundStr, "%d", round);
        string strRound =string (roundStr);
        vector<pair<double,pair<size_t , pair<bool, int> > > > allNodes;
        vector<pair< pair< int , int> , pair<double,int> > > frequencyArray;
        size_t validNodes=0;
        for (NodeID i =1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if(!n.isValid())
                        continue;
                validNodes++;
                double kmerCoverage=n.getNodeKmerCov();
                int marginalLength=n.getMarginalLength();//+ settings.getK();
                bool correctNode=false;
                double nodeMultiplicityR=1;
                //#ifdef DEBUG
                nodeMultiplicityR=trueMult[abs(n.getNodeID())];
                if (nodeMultiplicityR>0)
                        correctNode=true;
                //#endif

                allNodes.push_back( make_pair(kmerCoverage, make_pair( nodeMultiplicityR, make_pair(correctNode,marginalLength)) ));
        }
        std::sort (allNodes.begin(), allNodes.end());
        unsigned int i=0;
        double Interval=.1;
        double  St=allNodes[0].first;
        while (i<allNodes.size()-1) {
                double intervalCount=0;
                int sumOfMarginalLenght=0;
                int correctNodeMarginalLength=0;
                while (allNodes[i].first>=St && allNodes[i].first<St+Interval&& i<allNodes.size()-1) {
                        intervalCount++;
                        i++;
                        sumOfMarginalLenght=sumOfMarginalLenght+allNodes[i].second.second.second ;
                        if (allNodes[i].second.second.first)
                                correctNodeMarginalLength=correctNodeMarginalLength+allNodes[i].second.second.second;
                }
                double representative=St+Interval/2;
                if (intervalCount!=0)
                        frequencyArray.push_back(make_pair(make_pair( sumOfMarginalLenght, correctNodeMarginalLength) , make_pair(representative,intervalCount)));
                St=St+Interval;
        }
        plotCovDiagram(frequencyArray);

        cout<<"certainValue: "<<this->certainVlueCov<<endl;
        cout<<"safeValue: "<<this->safeValueCov<<endl;
        cout<<"redLineValue: "<<this->redLineValueCov<<endl;
        //writing to file to make a plot
        //#ifdef DEBUG
        //if (round>0)

        //#endif
        //return true;
}
void DBGraph::plotCovDiagram(vector<pair< pair< int , int> , pair<double,int> > >& frequencyArray){

        ofstream expcovFile;
        ofstream sexpcovFile;
        std::string roundstr="";
        if (updateCutOffValueRound<10)
                roundstr="00"+std::to_string( updateCutOffValueRound);
        else if (updateCutOffValueRound<100)
                roundstr="0"+std::to_string( updateCutOffValueRound);
        else
                roundstr=std::to_string( updateCutOffValueRound);
        string filename=settings.getTempDirectory()+ "cov/cov_"+roundstr+".dat";
        string sFileName=settings.getTempDirectory()+"cov/scov_"+roundstr+".dat";
        string outputFileName= settings.getTempDirectory()+"cov/cov_"+roundstr+".pdf";
        expcovFile.open( filename.c_str());
        sexpcovFile.open(sFileName.c_str());
        expcovFile<<"certainValue:"<<this->certainVlueCov<<endl;
        expcovFile<<"safeValue:"<<this->safeValueCov<<endl;
        expcovFile<<"redLineValue:"<<this->redLineValueCov<<endl;
        expcovFile<<"cutOffvalue:"<<this->cutOffvalue<<endl;
        expcovFile<<"estimatedKmerCoverage:"<<this->estimatedKmerCoverage<<endl;
        expcovFile<<"estimatedMKmerCoverageSTD:"<<this->estimatedMKmerCoverageSTD<<endl;
        expcovFile<<"representative_1\trepresentative_2\tsumOfMarginalLength\tsumOfcorrectLength\tsumOfIncorrectLength"<<endl;
        sexpcovFile<<"representative_1,representative_2,sumOfcorrectLength,sumOfIncorrectLength"<<endl;
        int i=0;
        double added=(frequencyArray[1].second.first-frequencyArray[0].second.first)/2;
        while(i<frequencyArray.size()) {
                pair<double,int> pre=frequencyArray[i].second;
                expcovFile<<pre.first<< "\t"<<pre.first+added<<"\t"<<frequencyArray[i].first.first<<"\t"<< frequencyArray[i].first.second<<"\t"<< frequencyArray[i].first.first-frequencyArray[i].first.second <<endl;
                sexpcovFile<<pre.first<< ","<< pre.first+added<<","<<frequencyArray[i].first.second<<","<< frequencyArray[i].first.first-frequencyArray[i].first.second <<endl;
                i++;
        }
        sexpcovFile.close();
        expcovFile.close();
        if (updateCutOffValueRound==1){
                string correctCluster=settings.getTempDirectory()+"correcNodeCluster.dta";
                string erroneousCluster=settings.getTempDirectory()+"erroneousCluster.dat";
                ExpMaxClustering exp(sFileName,erroneousCluster,correctCluster,.01,1, 50 );
                exp.doClassification();
                this->redLineValueCov= exp.intersectionPoint*.95;
                this->safeValueCov=this->redLineValueCov*.85;
                this->certainVlueCov=this->redLineValueCov*.75;

        }
        this->updateCutOffValueRound++;
        #ifdef DEBUG
        string address =" plot.dem";
        string command="gnuplot -e  \"filename='"+filename+"'\" -e \"outputfile='"+outputFileName+"'\""+address;
        system(command.c_str());
        #endif



}

struct readStructStr
{
    int intID;
    string strID;
    string corrctReadContent;
    string erroneousReadContent;
    string qualityProfile ;
    string orientation;

};


bool DBGraph::filterCoverage(float round)
{
        cout<<"*********************<<Filter Coverage starts>>......................................... "<<endl;
        int tp=0;
        int tn=0;
        int fp=0;
        int fn=0;
        size_t numFiltered = 0;
        for (NodeID i =1; i <= numNodes; i++) {

                SSNode n = getSSNode(i);
                if(!n.isValid())
                        continue;
                if (n.getNodeKmerCov()>this->certainVlueCov) {
                        if (trueMult[abs(i)] >= 1)
                                tn++;
                        else
                                fn++;
                        continue;
                }
                bool Del=removeNode(n);
                if (Del){
                        numFiltered++;
                        if (trueMult[abs(i)] >= 1)
                                fp++;
                        else
                                tp++;
                }
                else
                {
                        if (trueMult[abs(i)] >= 1)
                                tn++;
                        else
                                fn++;
                }
        }
        cout << "Number of coverage nodes deleted: " << numFiltered<<endl;
#ifdef DEBUG
        cout << " Gain value is ("<<100*((double)(tp-fp)/(double)(tp+fn))<< "%)"<<endl;
        cout<< "TP:	"<<tp<<"	TN:	"<<tn<<"	FP:	"<<fp<<"	FN:	"<<fn<<endl;
        cout << "Sensitivity: ("<<100*((double)tp/(double)(tp+fn))<<"%)"<<endl;
        cout<<"Specificity: ("<<100*((double)tn/(double)(tn+fp))<<"%)"<<endl;
#endif

        return numFiltered > 0;

}
bool DBGraph::removeNode(SSNode & rootNode) {
       //if (rootNode.getMarginalLength()>maxNodeSizeToDel)
       //       return false;

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

bool DBGraph::removeIncorrectLink()
{
        cout<<"*********************<<Remove Incorrect Link starts>>......................................... "<<endl;
        int num=0;
        bool modify=false;
        for(int i=1; i<numNodes; i++) {
                SSNode left=getSSNode(i);
                if ( !left.isValid() )
                        continue;
                pair<int, pair<double,double> > result=nodesExpMult[abs( left.getNodeID())];
                double confidenceRatio=result.second.first;
                double inCorrctnessRatio=result.second.second;
                double nodeMultiplicity=result.first;
                int j=0;
                if(nodeMultiplicity==1&& ( confidenceRatio/inCorrctnessRatio)>10 && (left.getNumLeftArcs()>1||left.getNumRightArcs()>1)) {
                        if (left.getNumLeftArcs()>1) {
                                j=0;
                                for ( ArcIt it = left.leftBegin(); it != left.leftEnd(); it++ ) {
                                        SSNode llNode = getSSNode ( it->getNodeID() );
                                        pair<int, pair<double,double> > resultL=nodesExpMult[abs( llNode.getNodeID())];
                                        double confidenceRatioL=result.second.first;
                                        double inCorrctnessRatioL=result.second.second;
                                        double nodeMultiplicityL=result.first;
                                        if(nodeMultiplicityL==1&& (confidenceRatioL/inCorrctnessRatioL)>10)
                                                j++;
                                }
                                if (j>1) {

                                        SSNode minNode=left.leftBegin();
                                        int leftCov=left.getLeftArc(minNode.getNodeID())->getCoverage();
                                        for ( ArcIt it = left.leftBegin(); it != left.leftEnd(); it++ ) {
                                                SSNode llNode = getSSNode ( it->getNodeID() );
                                                if(left.getLeftArc(llNode.getNodeID())->getCoverage()<leftCov) {
                                                        leftCov=left.getLeftArc(llNode.getNodeID())->getCoverage();
                                                        minNode=llNode;
                                                }
                                        }
                                        if (abs(left.getNodeID())!=abs(minNode.getNodeID())) {
                                                bool dresult = left.deleteLeftArc ( minNode.getNodeID() );
                                                assert ( dresult );
                                                dresult=minNode.deleteRightArc(left.getNodeID());
                                                assert(dresult);
                                                modify=true;
                                                num++;
                                        }
                                }
                        }
                        j=0;
                        if (left.getNumRightArcs()>1) {
                                for ( ArcIt it2 = left.rightBegin(); it2 != left.rightEnd(); it2++ ) {

                                        SSNode rrNode = getSSNode ( it2->getNodeID() );
                                        pair<int, pair<double,double> > resultR=nodesExpMult[abs( rrNode.getNodeID())];
                                        double confidenceRatioR=resultR.second.first;
                                        double inCorrctnessRatioR=resultR.second.second;
                                        double nodeMultiplicityR=resultR.first;
                                        if(nodeMultiplicityR==1&& ( confidenceRatioR/inCorrctnessRatioR)>10)
                                                j++;
                                }
                                if (j>1) {
                                        SSNode minNode=left.rightBegin();
                                        int rightCov=left.getRightArc(minNode.getNodeID())->getCoverage();
                                        for ( ArcIt it = left.rightBegin(); it != left.rightEnd(); it++ ) {
                                                SSNode rrNode = getSSNode (it->getNodeID());
                                                if(left.getRightArc(rrNode.getNodeID())->getCoverage()<rightCov) {
                                                        rightCov=left.getRightArc(rrNode.getNodeID())->getCoverage();
                                                        minNode=rrNode;
                                                }
                                        }
                                        if (abs(left.getNodeID())!=abs(minNode.getNodeID())) {
                                                bool dresult = left.deleteRightArc ( minNode.getNodeID() );
                                                assert ( dresult );
                                                dresult=minNode.deleteLeftArc(left.getNodeID());
                                                assert(dresult);
                                                modify=true;
                                                num++;
                                        }
                                }
                        }
                }

        }
        cout<<"The number of deleted arcs are:"<<num<<endl;
        return modify;

}

void DBGraph::sanityCheck()
{
    // sanity check
    for (NodeID i = -numNodes; i <= numNodes; i++) {
        if (i == 0)     // skip node 0, doesn't exist
            continue;

        // shortcuts
        SSNode n = getSSNode(i);

        string sequence = n.getSequence();
        if (!n.isValid()) {
            if ((n.getNumLeftArcs() != 0) || (n.getNumRightArcs() != 0))
                cerr << "\t\tNode " << n.getNodeID()
                     << " is invalid but has arcs." << endl;
        }
        // check the continuity of the kmers
        Kmer firstKmer(sequence);
        for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
            SSNode t = getSSNode(it->getNodeID());
            if (!t.isValid())
                cerr << "\t\tNode " << n.getNodeID()
                     << " has a left arc to an invalid node." << endl;

            assert(t.getRightArc(i)->getCoverage() == it->getCoverage());
            //  assert(it->getCoverage() > 0);
            Kmer tKmer = t.getRightKmer();
            tKmer.pushNucleotideRight(firstKmer.peekNucleotideRight());
            if (firstKmer != tKmer)
                cerr << "\tError in the continuity of the nodes" << endl;
            assert(firstKmer == tKmer);
        }
        Kmer lastKmer(sequence, sequence.size()-Kmer::getK());
        for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
            SSNode t = getSSNode(it->getNodeID());
            if (!t.isValid())
                cerr << "\t\tNode " << n.getNodeID()
                     << " has a right arc to an invalid node." << endl;

            assert(t.getLeftArc(i)->getCoverage() == it->getCoverage());

            Kmer tKmer = t.getLeftKmer();
            tKmer.pushNucleotideLeft(lastKmer.peekNucleotideLeft());

            if (lastKmer != tKmer)
                cerr << "\tError in the continuity of the nodes" << endl;
            assert (lastKmer == tKmer);
        }
    }
}

DBGraph::~DBGraph()
{
    delete [] arcs;
    delete [] nodes;
}

bool DBGraph::getLeftUniqueSSNode(const SSNode &node, SSNode &leftNode) const
{
    if (node.getNumLeftArcs() != 1)
        return false;

    leftNode = getSSNode(node.leftBegin()->getNodeID());

    if (leftNode.getNumRightArcs() > 1)
        return false;

    return true;
}

bool DBGraph::getRightUniqueSSNode(const SSNode &node, SSNode &rightNode) const
{
    if (node.getNumRightArcs() != 1)
        return false;

    rightNode = getSSNode(node.rightBegin()->getNodeID());

    if (rightNode.getNumLeftArcs() > 1)
        return false;

    return true;
}

bool DBGraph::simplyfyGraph()
{
    // count the number of nodes
    int numInitial = 0;
    for (NodeID i = 1; i <= numNodes; i++) {

        if (nodes[i].isValid())
            numInitial++;
    }
    while (clipTips(false));

    // count the number of clipped nodes
    int numRemaining = 0;
    for (NodeID i = 1; i <= numNodes; i++)
        if (nodes[i].isValid())
            numRemaining++;

    size_t numClipped = numInitial - numRemaining;

    cout << "Clipped " << numClipped << "/" << numInitial << " nodes." << endl;
    return numClipped > 0;
}


void DBGraph::increaseCoverage(const NodeEndRef &left, const NodeEndRef &right)
{
    SSNode lNode = getSSNode(left.getNodeID());
    SSNode rNode = getSSNode(right.getNodeID());

    // Avoid increasing the coverage of a non-existing intra-node arc.
    // This might arise when the node consists of two overlapping kmers
    // or when (the first part of) the node is a palindrome
    if (abs(left.getNodeID()) == abs(right.getNodeID()))
        if (lNode.getLength() > Kmer::getK())
            if (lNode.getLeftKmer() == left.getKmer())
                return;

    lNode.getRightArc(right.getNodeID())->incReadCov();
    rNode.getLeftArc(left.getNodeID())->incReadCov();
}

void DBGraph::threadThroughReads(const string& filename, size_t numReads,
                                 const NodeEndTable &table)
{
    if (numReads == 0)
        return;

    ifstream ifs(filename.c_str(), ios::in | ios::binary);
    TString read;

    size_t readCount = 0;
    while (ifs.good()) {
        read.readFromStream(ifs);

        // read too short ?
        if (read.getLength() < Kmer::getK()) continue;

        if (readCount++ % OUTPUT_FREQUENCY == 0) {
            cout << "Processing read " << readCount << "/"
                 << numReads << "\r";
            cout.flush();
        }

        // process the initial kmer of the string
        Kmer kmer(read);
        NodeEndRef refLeft = table.find(kmer);

        // process the rest of the kmers in the string
        for (size_t i = Kmer::getK(); i < read.getLength(); i++) {
            kmer.pushNucleotideRight(read[i]);
            NodeEndRef refRight = table.find(kmer);

            if (refLeft.first != table.end() &&
                    refRight.first != table.end())
                increaseCoverage(refLeft, refRight);

            refLeft = refRight;
        }
    }

    ifs.close();

    cout << "Processing read " << readCount << "/"
         << numReads << endl;
}

void DBGraph::threadThroughReads()
{
    NodeEndTable table(settings.isDoubleStranded(), 2*numNodes);

    for (NodeID id = 1; id <= numNodes; id++) {
        DSNode& node = getDSNode(id);
        Kmer firstKmer = node.getLeftKmer();
        Kmer finalKmer = node.getRightKmer();

        if (node.getLength() > 1) {
            table.insert(firstKmer, id);
            table.insert(finalKmer, id);
        } else
            table.insert(firstKmer, id);
    }

    /*   for (size_t i = 0; i < settings.getInputCount(); i++) {
           const Input &input = settings.getInput(i);

           cout << "Processing file " << i+1 << "/"
                << settings.getInputCount() << ": "
                << input.getFilename()
                << ", type: " << input.getFileType()
                << ", reads: " << input.getReadType() << endl;

           threadThroughReads(input.getIsolatedFilename(),
                              input.getNumIsolatedReads(),
                              table);
           threadThroughReads(input.getPairedFilename(),
                              input.getNumPairedReads(),
                              table);
       }*/
}

void DBGraph::createFromFile(const string& nodeFilename,
                             const string& arcFilename,
                             const string& metaDataFilename)
{
    // auxiliary variables
    string dS, descriptor;
    int dI, length;
    double expMult, readStartCov;

    // read the metadata
    ifstream metaDataFile(metaDataFilename.c_str());
    if (!metaDataFile)
        throw ios_base::failure("Can't open " + metaDataFilename);
    metaDataFile >> numNodes >> numArcs;
    metaDataFile.close();

    NodeEndTable table(settings.isDoubleStranded(), 2*numNodes);

    // A) create the nodes
    ifstream nodeFile(nodeFilename.c_str());
    if (!nodeFile)
        throw ios_base::failure("Can't open " + nodeFilename);

    nodes = new DSNode[numNodes+1];
    SSNode::setNodePointer(nodes);
    for (NodeID id = 1; id <= numNodes; id++) {
        // read the node info
        nodeFile >> dS >> dI >> length >> expMult >> readStartCov >> descriptor;

        DSNode& node = getDSNode(id);
        node.setSequence(descriptor);

        node.setExpMult(expMult);
        node.setKmerCov(expMult);
        //comment by mahdi
        node.setReadStartCov(readStartCov);

        Kmer firstKmer = node.getLeftKmer();
        Kmer finalKmer = node.getRightKmer();

        if (node.getMarginalLength() > 1) {
            if (!table.insert(firstKmer, id))
                cerr << "ERROR: Multiple nodes start/end with the same k-mer!" << endl;
            if (!table.insert(finalKmer, id))
                cerr << "ERROR: Multiple nodes start/end with the same k-mer!" << endl;
        } else if (!table.insert(firstKmer, id))
            cerr << "ERROR: Multiple nodes start/end with the same k-mer!" << endl;
    }

    nodeFile.close();

    // B) create the arcs
    ifstream arcFile(arcFilename.c_str());
    if (!arcFile)
        throw ios_base::failure("Can't open " + arcFilename);

    // +2 because index 0 isn't used, final index denotes 'end'.
    arcs = new Arc[numArcs+2];
    DSNode::setArcsPointer(arcs);

    ArcID arcOffset = 1;
    for (NodeID i = 1; i <= numNodes; i++) {
        int dI, bfLeft, bfRight;

        // read the arc info
        arcFile >> dI >> bfLeft >> bfRight;
        KmerOverlap overlap((bfLeft << 4) + bfRight);

        int numLeftArcs = overlap.getNumLeftOverlap();
        int numRightArcs = overlap.getNumRightOverlap();

        DSNode& node = getDSNode(i);

        node.setNumLeftArcs(numLeftArcs);
        node.setNumRightArcs(numRightArcs);

        node.setFirstLeftArcID(arcOffset);
        node.setFirstRightArcID(arcOffset+numLeftArcs);

        // connect the left arcs
        Kmer firstKmer = node.getLeftKmer();
        for (NucleotideID j = 0; j < 4; j++) {
            char n = Nucleotide::nucleotideToChar(j);
            if (!overlap.hasLeftOverlap(n))
                continue;

            Kmer kmer = firstKmer;
            kmer.pushNucleotideLeft(n);

            NodeEndRef ref = table.find(kmer);
            if (ref.first == table.end())
                throw ios_base::failure("Mismatch between nodes"
                                        "and arc file.");
            int arcCov;
            arcFile >> arcCov;
            arcs[arcOffset].setCoverage(arcCov);
            arcs[arcOffset++].setNodeID(ref.getNodeID());
        }

        // connect the right arcs
        Kmer finalKmer = node.getRightKmer();
        for (NucleotideID j = 0; j < 4; j++) {
            char n = Nucleotide::nucleotideToChar(j);
            if (!overlap.hasRightOverlap(n))
                continue;

            Kmer kmer = finalKmer;
            kmer.pushNucleotideRight(n);

            NodeEndRef ref = table.find(kmer);
            if (ref.first == table.end())
                throw ios_base::failure("Mismatch between nodes"
                                        "and arc file.");
            int arcCov;
            arcFile >> arcCov;
            arcs[arcOffset].setCoverage(arcCov);
            arcs[arcOffset++].setNodeID(ref.getNodeID());
        }
    }

    arcFile.close();

    if(arcOffset != numArcs+1)
        throw ios_base::failure("Mismatch between nodes and arc file.");
}

void DBGraph::convertNodesToString(const vector<NodeID> &nodeSeq,
                                   string &output)
{
    output.clear();

    if (nodeSeq.empty())
        return;

    size_t size = 0;
    for (size_t i = 0; i < nodeSeq.size(); i++)
        size += getSSNode(nodeSeq[i]).getLength();
    size -= (nodeSeq.size() - 1) * (Kmer::getK() - 1);

    output = getSSNode(nodeSeq[0]).getSequence();
    for (size_t i = 1; i < nodeSeq.size(); i++)
        output.append(getSSNode(nodeSeq[i]).getSequence().substr(Kmer::getK() - 1));

    assert(size == output.size());
}

void DBGraph::convertNodesToString(const deque<SSNode> &nodeSeq,
                                   string &output)
{
    output.clear();

    if (nodeSeq.empty())
        return;

    size_t size = 0;
    for (size_t i = 0; i < nodeSeq.size(); i++)
        size += nodeSeq[i].getLength();
    size -= (nodeSeq.size() - 1) * (Kmer::getK() - 1);

    output = nodeSeq[0].getSequence();
    for (size_t i = 1; i < nodeSeq.size(); i++)
        output.append(nodeSeq[i].getSequence().substr(Kmer::getK() - 1));

    assert(size == output.size());
}

void DBGraph::convertNodesToString(const vector<NodeID> &nodeSeq,
                                   int startPos,
                                   int stopPos,
                                   string &output)
{
    size_t size = stopPos + getSSNode(nodeSeq.front()).getMarginalLength() - startPos;
    for (size_t i = 1; i < nodeSeq.size() - 1; i++)
        size += getSSNode(nodeSeq[i]).getMarginalLength();

    for (size_t i = 0; i < nodeSeq.size(); i++)
        cout << nodeSeq[i] << " ";
    cout << endl;

    cout << "StartPos " << startPos << " stopPos " << stopPos << endl;

    cout << size << endl;

    output = getSSNode(nodeSeq[0]).getSequence().substr(startPos);
    for (size_t i = 1; i < nodeSeq.size()-1; i++)
        output.append(getSSNode(nodeSeq[i]).getSequence().substr(Kmer::getK() - 1));
    output.append(getSSNode(nodeSeq.back()).getSequence().substr(0, stopPos));

    cout << output << endl;
    cout << output.size() << endl;

    assert(size == output.size());
}

void DBGraph::writeGraph(const std::string& nodeFilename,
                         const std::string& arcFilename,
                         const std::string& metaDataFilename)
{
    vector<size_t> nodeLengths;
    ofstream nodeFile(nodeFilename.c_str());
    ofstream arcFile(arcFilename.c_str());

    size_t numExtractedNodes = 0, numExtractedArcs = 0;
    for (NodeID id = 1; id <= numNodes; id++) {
        SSNode node = getSSNode(id);

        if (!node.isValid())
            continue;

        numExtractedNodes++;

        //nodeFile << "NODE" << "\t" << numExtractedNodes << "\t"
        nodeFile << "NODE" << "\t" << id << "\t"
                 << node.getLength() << "\t" << node.getKmerCov()
                 << "\t" << node.getReadStartCov() << "\n"
                 << node.getSequence() << "\n";

        KmerOverlap ol;
        for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
            char c = getSSNode(it->getNodeID()).getRightKmer().peekNucleotideLeft();
            ol.markLeftOverlap(c);
        }

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
            char c = getSSNode(it->getNodeID()).getLeftKmer().peekNucleotideRight();
            ol.markRightOverlap(c);
        }

        arcFile << numExtractedNodes << "\t" << (int)ol.getLeftOverlap()
                << "\t" << (int)ol.getRightOverlap();

        for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++)
            arcFile << "\t" << it->getCoverage();

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++)
            arcFile << "\t" << it->getCoverage();

        arcFile << endl;

        numExtractedArcs += ol.getNumLeftOverlap() + ol.getNumRightOverlap();

        nodeLengths.push_back(node.getLength());
    }

    nodeFile.close();
    arcFile.close();

    ofstream metadataFile(metaDataFilename.c_str());
    metadataFile << numExtractedNodes << "\t" << numExtractedArcs << endl;
    metadataFile.close();

    cout << "Extracted " << numExtractedNodes << " nodes and "
         << numExtractedArcs << " arcs." << endl;

    sort(nodeLengths.begin(), nodeLengths.end());

    size_t totalLength = 0;
    for (size_t i = 0; i < nodeLengths.size(); i++)
        totalLength += nodeLengths[i];

    size_t currLength = 0;
    for (size_t i = 0; i < nodeLengths.size(); i++) {
        currLength += nodeLengths[i];
        if (currLength >= 0.5*totalLength) {
            cout << endl << "N50 is " << nodeLengths[i] << " (total length: " << totalLength << ")" << endl;
            cout << "This was found at node " << i << "/" << numNodes << endl;
            break;
        }
    }

    cout << "The largest node contains " << nodeLengths.back() << " basepairs." << endl;
}

size_t DBGraph::getN50()
{
    vector<size_t> nodeLengths;

    sizeOfGraph=0;
    size_t numExtractedNodes = 0, numExtractedArcs = 0;
    for (NodeID id = 1; id <= numNodes; id++) {
        SSNode node = getSSNode(id);

        if (!node.isValid())
            continue;

        numExtractedNodes++;
        //check later for adding kmerSize
        sizeOfGraph=sizeOfGraph+node.getMarginalLength()+kmerSize;
        KmerOverlap ol;
        for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
            char c = getSSNode(it->getNodeID()).getRightKmer().peekNucleotideLeft();
            ol.markLeftOverlap(c);
        }

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
            char c = getSSNode(it->getNodeID()).getLeftKmer().peekNucleotideRight();
            ol.markRightOverlap(c);
        }
        numExtractedArcs += ol.getNumLeftOverlap() + ol.getNumRightOverlap();

        nodeLengths.push_back(node.getLength());
    }

    numberOfValidNodes=numExtractedNodes;

    cout<<"size of graph: "<<sizeOfGraph<<endl;
    cout<<"number of valid Node: "<<numberOfValidNodes<<endl;

    cout << "Extracted " << numExtractedNodes << " nodes and "
         << numExtractedArcs << " arcs." << endl;

    sort(nodeLengths.begin(), nodeLengths.end());

    size_t totalLength = 0;
    for (size_t i = 0; i < nodeLengths.size(); i++)
        totalLength += nodeLengths[i];

    size_t currLength = 0;
    size_t n50;
    for (size_t i = 0; i < nodeLengths.size(); i++) {
        currLength += nodeLengths[i];
        if (currLength >= 0.5*totalLength) {
            cout << endl << "N50 is " << nodeLengths[i] << " (total length: " << totalLength << ")" << endl;
            cout << "This was found at node " << i << "/" << numNodes << endl;
            n50= nodeLengths[i];
            break;
        }
    }
    cout << "The largest node contains " << nodeLengths.back() << " basepairs." << endl;
    return n50;

}

void DBGraph::populateTable() {
        table = new KmerNodeTable(settings, numNodes);
        table->populateTable ( nodes );
}
bool DBGraph::kmerExistsInGraph(Kmer const &kmer) const {
        return getNodePosPair(kmer).isValid();
}
NodePosPair DBGraph::getNodePosPair(Kmer const &kmer) const {
        return table->find(kmer);
}
double DBGraph::getReadLength() const {
        return readLength;
}

void DBGraph::writeGraphExplicit() const
{
        if (!settings.getSkipStage5()) {
        	return;
        }
        vector<size_t> nodeLengths;
        ofstream nodeFile("DBGraph.txt");

        size_t numExtractedNodes = 0, numExtractedArcs = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);

                if (!node.isValid())
                        continue;

                numExtractedNodes++;

                nodeFile << "NODE" << "\t" << numExtractedNodes << "\t"
                         << node.getLength() << "\n"
                         << node.getSequence() << "\n";

                KmerOverlap ol;
                nodeFile << "LEFT ARCS" << "\t" << (int)node.getNumLeftArcs();
                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
                        nodeFile << "\t" << it->getNodeID();
                }

                nodeFile << "\nRIGHT ARCS" << "\t" << (int)node.getNumRightArcs();
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        nodeFile << "\t" << it->getNodeID();
                }
                nodeFile << "\n";
        }

        nodeFile.close();
}
