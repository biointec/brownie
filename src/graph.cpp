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
#include "kmeroverlap.h"
#include "settings.h"
#include "nodeendstable.h"
#include "kmernode.h"
#include "library.h"
#include "component.h"
#include <sys/stat.h>
using namespace std;

DSNode* SSNode::nodes = NULL;
const DBGraph* DBGraph::graph = NULL;


DBGraph::DBGraph(const Settings& settings) : table(NULL), settings(settings),
        nodes(NULL), arcs(NULL), numNodes(0), numArcs(0), mapType(SHORT_MAP) {
    DBGraph::graph = this;
    //mahdi comment my
    initialize();
}

void DBGraph::initialize()
{
    //correct it later it should read this information from setting file, but setting dosn't have such information right now

    readLength=100;//settings.getReadLength();
    coverage=100;//settings.getCoverage();
    minCertainVlueCov=2;
    minSafeValueCov=2.5;
    minRedLineValueCov=3;
    estimatedKmerCoverage=(coverage/readLength)*(readLength-Kmer::getK()+1)/12;
    estimatedMKmerCoverageSTD=sqrt( estimatedKmerCoverage);
    certainVlueCov=minCertainVlueCov;
    safeValueCov=minSafeValueCov;
    redLineValueCov=minRedLineValueCov;
    updateCutOffValueRound=1;
    maxNodeSizeToDel=readLength*4;
    cutOffvalue=redLineValueCov;
    stepSize=.2;
    redLineUpsize=1.5;
    redLineDownSize=1;
    safeValueDownSize=.7;
    safeValueUpSize=1;
    certainValueUpSize=.7;
    certainValueDownSize=.3;
}



/**
 * this routine manipulate graph to detect erroneous nodes and delete them
 * the firs routine is filter coverage it removes nodes with coverage 1
 * then Tips will be deleted
 * then bubbles will be deleted
 * some Statistic will be calculated in evey round to find the multiplicity of nodes and our confidence about our guess
 * then Unreliable nodes will be detected based on our Statistic parameter in pre step
 * while loop continue until no changes would be possible
 * cut off value will be inreased by one to delete nodes with low coverage
 * if the size fo graph is lower than genome size modification stops
 *
 *
 */
void DBGraph::graphPurification(string trueMultFilename,
                                const LibraryContainer& libraries)
{
        int round=1;
        size_t maxBubbleDepth=maxNodeSizeToDel;
        size_t increamentDepth = readLength;
        bool simplified = true;
        while (simplified ) {// &&
                //
                //*******************************************************
                updateCutOffValue(round);
                bool tips=clipTips(round);
                if(tips) {
                        mergeSingleNodes(true);
                        #ifdef DEBUG
                        compareToSolution(trueMultFilename, false);
                        #endif
                }
                //*******************************************************
                bool bubble=false;
                size_t depth=settings.getK();
                #ifdef DEBUG
                compareToSolution(trueMultFilename, false);
                #endif
                cout << endl << " ================= Bubble Detection ==================" << endl;
                bubble= bubbleDetection(depth);
                mergeSingleNodes(true);
                bool continuEdit=true;
                size_t maxDepth=(round)*increamentDepth>maxBubbleDepth?maxBubbleDepth:(round)*increamentDepth;
                while(depth<maxDepth&& continuEdit){
                        depth=depth+increamentDepth;
                        cout << "Bubble depth: " << depth << endl;
                        continuEdit= bubbleDetection(depth);
                        if (continuEdit)
                                mergeSingleNodes(true);
                        bubble=false?continuEdit:bubble;
                }
                #ifdef DEBUG
                compareToSolution(trueMultFilename,false);
                #endif
                extractStatistic(round);
                cout << endl << " ============= Delete Unreliable Nodes  ==============" << endl;
                bool deleted=deleteUnreliableNodes();
                mergeSingleNodes(true);
                continuEdit=deleted;
                while(continuEdit){
                        continuEdit=deleteUnreliableNodes();
                        mergeSingleNodes(true);
                        extractStatistic(round);
                }
                while(mergeSingleNodes(true));
                #ifdef DEBUG
                compareToSolution(trueMultFilename,false);
                cout<<"estimated Kmer Coverage Mean: "<<estimatedKmerCoverage<<endl;
                cout<<"estimated Kmer Coverage STD: "<<estimatedMKmerCoverageSTD<<endl;
                #endif
                simplified = tips    || deleted || bubble;
                updateGraphSize();
                round++;
        }
}

/*This routine calculate the value of CutOff-cov in the first call,
 * based on cutOff value three other values will be updated in each round of calling. they increase in each round
 * cutOffvalue*.3< certainValue< cutOffvalue*.8 : the most conservative value , used to remove nodes only based on coverage
 * cutOffvalue*.7< safeValue <cutOffvalue  always lower than cutOff , used for remoing joined tips or isolated nodes.
 * cutOffvalue <redLineValue< cutOffvalue*1.5, increasing incrementally in while loop
 *
 *
 *
 */
void DBGraph::updateCutOffValue(int round)
{
        /*............read line ...........*/
        double ratio=redLineDownSize;
        size_t roundLimit=(redLineUpsize-redLineDownSize)/stepSize;
        if (round<=roundLimit)
                ratio=ratio+(double)round*stepSize;
        else
                ratio=redLineUpsize;
        redLineValueCov= cutOffvalue*ratio;

        /*............safe value ...........*/
        ratio=safeValueDownSize;
        roundLimit=(safeValueUpSize-safeValueDownSize)/stepSize;
        if (round<=roundLimit)
                ratio=ratio+(double)round*stepSize;
        else
                ratio=safeValueUpSize;
        safeValueCov=cutOffvalue*ratio;

        /*............certainValue ...........*/
        ratio=certainValueDownSize;
        roundLimit=(certainValueUpSize-certainValueDownSize)/stepSize ;
        if (round<=roundLimit)
                ratio=ratio+(double)round*stepSize;
        else
                ratio=certainValueUpSize;
        certainVlueCov=cutOffvalue*ratio;

        updateCutOffValueRound++;
        cout << "Certain value: " << this->certainVlueCov << endl;
        cout << "Safe value: " << this->safeValueCov << endl;
        cout << "Red line value: " << this->redLineValueCov << endl;
}



void DBGraph::plotCovDiagram(){
        string dir=settings.getTempDirectory()+"cov";
        const int dir_err = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

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
                #ifdef DEBUG
                nodeMultiplicityR=trueMult[abs(n.getNodeID())];
                if (nodeMultiplicityR>0)
                        correctNode=true;
                #endif
                allNodes.push_back( make_pair(kmerCoverage, make_pair( nodeMultiplicityR, make_pair(correctNode,marginalLength)) ));
        }
        std::sort (allNodes.begin(), allNodes.end());
        size_t  i=0;
        double Interval=.5;
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
                if (St>estimatedKmerCoverage+estimatedMKmerCoverageSTD*10)
                        break;
        }
        //wrriting in file

        ofstream sexpcovFile;
        std::string roundstr="";
        if (updateCutOffValueRound<10)
                roundstr="00"+std::to_string( updateCutOffValueRound);
        else if (updateCutOffValueRound<100)
                roundstr="0"+std::to_string( updateCutOffValueRound);
        else
                roundstr=std::to_string( updateCutOffValueRound);

        string sFileName=dir+"/scov_"+roundstr+".dat";
        sexpcovFile.open(sFileName.c_str());
        sexpcovFile<<"representative_1\trepresentative_2\tsumOfCorrectLength\tsumOfIncorrectLength"<<endl;
        i=0;
        double added=(frequencyArray[1].second.first-frequencyArray[0].second.first)/2;
        while(i<frequencyArray.size()) {
                pair<double,int> pre=frequencyArray[i].second;
                sexpcovFile<<pre.first<< "\t"<< pre.first+added<<"\t"<<frequencyArray[i].first.second<<"\t"<< frequencyArray[i].first.first-frequencyArray[i].first.second <<endl;
                i++;
        }
        sexpcovFile.close();
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
        read.read(ifs);

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
        ofstream nodeFile(nodeFilename.c_str());
        ofstream arcFile(arcFilename.c_str());

        size_t numExtractedNodes = 0, numExtractedArcs = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);

                if (!node.isValid())
                        continue;

                numExtractedNodes++;

                nodeFile << ">NODE" << "\t" << id << "\t"
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
        }

        nodeFile.close();
        arcFile.close();

        ofstream metadataFile(metaDataFilename.c_str());
        metadataFile << numExtractedNodes << "\t" << numExtractedArcs << endl;
        metadataFile.close();

        cout << "Extracted " << numExtractedNodes << " nodes and "
             << numExtractedArcs << " arcs." << endl;
}

void DBGraph::writeGraphBin(const std::string& nodeFilename,
                            const std::string& arcFilename,
                            const std::string& metaDataFilename)
{
        // A) Write node file
        ofstream nodeFile(nodeFilename.c_str(), ios::binary);
        for (NodeID id = 1; id <= numNodes; id++) {
                DSNode& node = getDSNode(id);
                // write the node contents
                node.write(nodeFile);
        }
        nodeFile.close();

        // B) Write arc file
        ofstream arcFile(arcFilename.c_str(), ios::binary);
        for (ArcID i = 0; i < numArcs+2; i++)
                arcs[i].write(arcFile);
        arcFile.close();

        // C) Write metadata file
        ofstream metadataFile(metaDataFilename.c_str());
        metadataFile << numNodes << "\t" << numArcs << endl;
        metadataFile.close();

        cout << "Wrote " << numNodes << " nodes and "
             << numArcs << " arcs" << endl;
}

void DBGraph::loadGraphBin(const std::string& nodeFilename,
                           const std::string& arcFilename,
                           const std::string& metaDataFilename)
{
        // auxiliary variables
        string dS, descriptor;

        // read the metadata
        ifstream metaDataFile(metaDataFilename.c_str());
        if (!metaDataFile)
                throw ios_base::failure("Can't open " + metaDataFilename);
        metaDataFile >> numNodes >> numArcs;
        metaDataFile.close();

        // A) create the nodes
        ifstream nodeFile(nodeFilename.c_str(), ios::binary);
        if (!nodeFile)
                throw ios_base::failure("Can't open " + nodeFilename);

        nodes = new DSNode[numNodes+1];
        SSNode::setNodePointer(nodes);
        for (NodeID id = 1; id <= numNodes; id++) {
                // read the node info

                DSNode& node = getDSNode(id);
                node.read(nodeFile);
        }
        nodeFile.close();

        // B) create the arcs
        ifstream arcFile(arcFilename.c_str());
        if (!arcFile)
                throw ios_base::failure("Can't open " + arcFilename);

        // +2 because index 0 isn't used, final index denotes 'end'.
        arcs = new Arc[numArcs+2];
        DSNode::setArcsPointer(arcs);
        for (ArcID i = 0; i < numArcs+2; i++)
                arcs[i].read(arcFile);
        arcFile.close();
}

size_t DBGraph::updateGraphSize()
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
        sizeOfGraph = sizeOfGraph + node.getMarginalLength() + Kmer::getK() - 1;
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


#ifdef DEBUG
    cout<<"size of graph: "<<sizeOfGraph<<endl;
    cout<<"number of valid Node: "<<numExtractedNodes<<endl;
    cout << "Extracted " << numExtractedNodes << " nodes and "
         << numExtractedArcs << " arcs." << endl;
#endif
    sort(nodeLengths.begin(), nodeLengths.end(), std::greater<int>());

    size_t totalLength = 0;
    for (size_t i = 0; i < nodeLengths.size(); i++)
        totalLength += nodeLengths[i];
    size_t currLength = 0;
    size_t n50;
    for (size_t i = 0; i < nodeLengths.size(); i++) {
        currLength += nodeLengths[i];
        if (currLength >= 0.5*totalLength) {
#ifdef DEBUG
            cout << endl << "N50 is " << nodeLengths[i] << " (total length: " << totalLength << ")" << endl;
            cout << "This was found at node " << i << "/" << numNodes << endl;
#endif
            n50= nodeLengths[i];
            this->n50=n50;
            break;
        }
    }
#ifdef DEBUG
    cout << "The largest node contains " << nodeLengths.back() << " basepairs." << endl;
    cout <<"N50:"<<n50<<endl;
#endif
    return sizeOfGraph;

}
struct greater_than_Component
{
    inline bool operator() (const Component& object1, const Component& object2)
    {
        return (object1.Size > object2.Size);
    }
};
struct less_than_Component
{
    inline bool operator() (const Component& object1, const Component& object2)
    {
        return (object1.Size < object2.Size);
    }
};
void DBGraph::reportSta()
{
        int srcID=1;
        int i=0;
        set<NodeID> nodesHandled;
        set<NodeID> currentSetNodes;
        size_t numberOfcomponents=0;
        ofstream sexpcovFile;
        string dir=settings.getTempDirectory()+"Statistic";
        const int dir_err = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        string sFileName=dir+"/componentStatistics_" + std::to_string( updateCutOffValueRound)+".txt";
        sexpcovFile.open(sFileName.c_str());
        double nkcovAVG=0;
        vector<Component> components;
        cout<<"#Extracting information about the components in the graph"<<endl;

        for (i =srcID ; i <= numNodes; i++){
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
                        nkcovAVG=(nkcovAVG*(nodesHandled.size()-1)+thisNode.getNodeKmerCov())/nodesHandled.size();
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
                if (currentSetNodes.size()>0)
                {
                        Component component;
                        numberOfcomponents++;
                        getComponentSta(currentSetNodes, component);
                        components.push_back(component);
                }
                currentSetNodes.clear();
        }
        sexpcovFile<<endl<<"#number of nodes:\t"<<nodesHandled.size()<<endl;
        sexpcovFile<<"#mean of node coverage:\t"<<nkcovAVG<<endl;
        sexpcovFile<<"#number of disjoint components in this graph\t" <<numberOfcomponents<<endl;
        sexpcovFile<<endl<<"N10"<<'\t'<<"N30"<<'\t'<<"N50"<<'\t'<<"N70"<<'\t'<<"N90"<<'\t'<<"Nodes"<<'\t'<<"Arcs"<<'\t' <<"Size"<<'\t'<<"largest"<<'\t'<<"cov"<<endl;
        sort( components.begin(), components.end(), greater_than_Component());
        size_t num=1;
        for (auto component: components){
                sexpcovFile<<component.get_N(10)<<'\t'<<component.get_N(30)<<'\t'<<component.get_N(50)<<'\t'<<component.get_N(70)<<'\t'<<component.get_N(90)<<'\t'<<component.numOfNodes<<'\t'<<component.numOfArcs<<'\t'<<component.Size<<'\t' <<component.largestNodeSize<<'\t'<<component.nodeKmerCov <<endl;
                if (component.Size>1000){
                        makeN50Files(num, component);
                        num++;
                }
        }
        sexpcovFile.close();
        makeComponentPlotFile(components);
}
void DBGraph::makeN50Files(const size_t num,const Component component)
{
        string dir=settings.getTempDirectory()+"n50";
        const int dir_err = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::string roundstr="";
        if (num<10)
                roundstr="00"+std::to_string( num)+"_"+std::to_string( updateCutOffValueRound);
        else if (num<100)
                roundstr="0"+std::to_string( num)+"_"+std::to_string( updateCutOffValueRound);
        else
                roundstr=std::to_string( num);
        string sFileName=dir+"/N50_"+ roundstr+"_.dat";
        ofstream n50stream;

        n50stream.open(sFileName.c_str());
        n50stream<<"# size of component\t"<< component.Size<<endl;
        n50stream<<"# number of nodes\t"<< component.numOfNodes<<endl;
        n50stream<<"# number of arcs\t"<< component.numOfArcs<<endl;

        n50stream<<"N10\t" <<component.get_N(10)<<"\nN20\t" <<component.get_N(20)<<"\nN30\t" <<component.get_N(30)<<"\nN40\t" <<component.get_N(40)<<"\nN50\t" <<component.get_N(50)
        <<"\nN60\t" <<component.get_N(60)<<"\nN70\t" <<component.get_N(70)<<"\nN80\t" <<component.get_N(80)<<"\nN90\t" <<component.get_N(90);
        n50stream.close();
}

void DBGraph::makeComponentPlotFile(vector<Component> components)
{
        string dir=settings.getTempDirectory()+"Statistic";
        const int dir_err = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        ofstream sexpcovFile;
        string sFileName=dir+"/components_"+std::to_string( updateCutOffValueRound)+".dat";
        sexpcovFile.open(sFileName.c_str());
        size_t  i=0;
        size_t Interval=100;
        sort( components.begin(), components.end(), less_than_Component());
        size_t  St=components[0].Size;
        vector<pair< size_t , vector<Component> > > frequencyArray;
        while (i<components.size()) {
                double intervalCount=0;
                int sumOfMarginalLenght=0;
                vector<Component> componentsInBin;
                while (i<components.size()&&components[i].Size>=St && components[i].Size<St+Interval) {
                        componentsInBin.push_back(components[i]);
                        sumOfMarginalLenght=sumOfMarginalLenght+components[i].Size ;
                        intervalCount++;
                        i++;
                }
                size_t representative=St+Interval/2;
                if (intervalCount!=0)
                        frequencyArray.push_back(make_pair(representative, componentsInBin));
                St=St+Interval;
        }
        i=0;
        sexpcovFile<< "Index" << "\t"<<"Num"<<"\t"<<"avgSize"<<endl ;
        while(i<frequencyArray.size()) {
                vector<Component> bin=frequencyArray[i].second;
                size_t sumOfSize=0;
                for (auto it:bin){
                        sumOfSize=sumOfSize+it.Size;
                }
                sexpcovFile<<fixed<<frequencyArray[i].first << "\t"<< bin.size()<< "\t"<< sumOfSize/bin.size() <<endl;;
                i++;
        }
        sexpcovFile.close();
}


void DBGraph::getComponentSta(set<NodeID> &currentSetNodes, Component &component)
{
        component.numOfNodes=currentSetNodes.size();
        vector<size_t> nodeLengths;
        size_t componentSize =0;
        size_t numExtractedNodes = 0, numExtractedArcs = 0;
        double summedCoverage = 0;
        for (auto id : currentSetNodes) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;
                numExtractedNodes++;
                //add full size of the node, including overlap
                componentSize += node.getMarginalLength() + Kmer::getK() - 1;
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
                summedCoverage += node.getNodeKmerCov();
        }
        component.Size=componentSize;
        component.numOfArcs=numExtractedArcs;
        sort(nodeLengths.begin(), nodeLengths.end(), std::greater<int>());
        size_t totalLength = 0;
        for (size_t i = 0; i < nodeLengths.size(); i++)
                totalLength += nodeLengths[i];
        size_t currLength = 0;
        size_t currentNum=1;
        for (size_t i = 0; i < nodeLengths.size(); i++) {
                currLength += nodeLengths[i];
                while (currLength >= (currentNum * 0.1) * totalLength ) {
                        component.set_N(10 * currentNum, nodeLengths[i]);
                        currentNum++;
                }
        }
        component.nodeKmerCov=round( summedCoverage / numExtractedNodes);
        component.largestNodeSize=nodeLengths[0];
}

void DBGraph::populateTable() {
        table = new KmerNodeTable(settings, numNodes);
        table->populateTable(nodes);
}

void DBGraph::depopulateTable() {
        delete table; table = NULL;
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

void DBGraph::writeGraphFasta() const
{
        vector<size_t> nodeLengths;
        string nodeFileName = settings.getTempDirectory() + "/DBGraph.fasta";
        ofstream nodeFile(nodeFileName);

        size_t numExtractedNodes = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);

                if (!node.isValid())
                        continue;

                numExtractedNodes++;

                nodeFile << ">NODE" << "\t" << numExtractedNodes << "\t"
                         << node.getLength();


                KmerOverlap ol;
                nodeFile << "\t" << (int)node.getNumLeftArcs();
                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
                        nodeFile << "\t" << it->getNodeID();
                }

                nodeFile << "\t" << (int)node.getNumRightArcs();
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        nodeFile << "\t" << it->getNodeID();
                }
                nodeFile << "\n" << node.getSequence() << "\n";;
        }

        nodeFile.close();
}
