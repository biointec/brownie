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

using namespace std;

DSNode* SSNode::nodes = NULL;
const DBGraph* DBGraph::graph = NULL;

// ============================================================================
// GRAPH STATISTICS
// ============================================================================

std::ostream &operator<<(std::ostream &out, const GraphStats &stats)
{
        out << "Number of nodes: " << stats.numNodes << "\n";
        out << "Number of arcs: " << stats.numArcs << "\n";
        out << "N50: " << stats.N50 << "\n";
        out << "Total size (kmers): " << stats.totMargLength;

        return out;
}

// ============================================================================
// GRAPH CLASS
// ============================================================================

DBGraph::DBGraph(const Settings& settings) : table(NULL), settings(settings),
        nodes(NULL), arcs(NULL), numNodes(0), numArcs(0), mapType(SHORT_MAP) {
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
    estimatedKmerCoverage=(coverage/readLength)*(readLength-Kmer::getK()+1)/12;
    estimatedMKmerCoverageSTD=sqrt( estimatedKmerCoverage);
    certainVlueCov=minCertainVlueCov;
    safeValueCov=minSafeValueCov;
    redLineValueCov=minRedLineValueCov;
    updateCutOffValueRound=1;
    maxNodeSizeToDel=readLength*4;
    cutOffvalue=redLineValueCov;
}

void DBGraph::removeNode(NodeID nodeID)
{
        SSNode node = getSSNode(nodeID);
        for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
                SSNode leftNode = getSSNode(it->getNodeID());
                if (leftNode.getNodeID()==-node.getNodeID())
                        continue;
                bool result = leftNode.deleteRightArc (node.getNodeID());
                assert(result);
        }

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                SSNode rightNode = getSSNode ( it->getNodeID() );
                if (rightNode.getNodeID()==-node.getNodeID())
                        continue;
                bool result = rightNode.deleteLeftArc(node.getNodeID());
                assert(result);
        }

        node.deleteAllRightArcs();
        node.deleteAllLeftArcs();
        node.invalidate();
}

NodeID DBGraph::getFirstValidNode(NodeID seed)
{
        for (NodeID id = seed; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                if (getSSNode(id).isValid())
                        return id;
        }
        return 0;
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
        size_t i=0;
        double added=(frequencyArray[1].second.first-frequencyArray[0].second.first)/2;
        while(i<frequencyArray.size()) {
                pair<double,int> pre=frequencyArray[i].second;
                expcovFile<<pre.first<< "\t"<<pre.first+added<<"\t"<<frequencyArray[i].first.first<<"\t"<< frequencyArray[i].first.second<<"\t"<< frequencyArray[i].first.first-frequencyArray[i].first.second <<endl;
                sexpcovFile<<pre.first<< ","<< pre.first+added<<","<<frequencyArray[i].first.second<<","<< frequencyArray[i].first.first-frequencyArray[i].first.second <<endl;
                i++;
        }
        sexpcovFile.close();
        expcovFile.close();
        //#ifdef DEBUG
        //string address =" plot.dem";
        //string command="gnuplot -e  \"filename='"+sFileName+"'\" -e \"outputfile='"+outputFileName+"'\""+address;
        //system(command.c_str());
        //cout<<command<<endl;
        //#endif

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

GraphStats DBGraph::getGraphStats()
{
        vector<size_t> nodeLength;

        size_t numNodes = 0;            // number of (valid) nodes in the graph
        size_t numArcs = 0;             // number of (valid) arcs in the graph
        size_t N50 = 0;                 // N50 of the nodes
        size_t totMargLength = 0;       // total marginal length of all nodes
        size_t totLength = 0;           // total length of all nodes

        for (NodeID id = 1; id <= this->numNodes; id++) {
                SSNode node = getSSNode(id);

                if (!node.isValid())
                        continue;

                numNodes++;

                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++)
                        numArcs++;

                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++)
                        numArcs++;

                totMargLength += node.getMarginalLength();
                totLength += node.getLength();
                nodeLength.push_back(node.getLength());
        }

        sort(nodeLength.begin(), nodeLength.end());

        size_t currLength = 0;
        for (size_t i = 0; i < nodeLength.size(); i++) {
                currLength += nodeLength[i];
                if (currLength >= 0.5*totLength) {
                        N50 = nodeLength[i];
                        break;
                }
        }

        GraphStats stats;
        stats.setMetrics(numNodes, numArcs, N50, totMargLength);
        return stats;
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
