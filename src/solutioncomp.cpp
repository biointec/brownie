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
#include "readfile/fastafile.h"
#include "settings.h"
extern "C" {
#include "suffix_tree.h"
}

void DBGraph::writeCytoscapeGraph(int ID)
{
    /*stringstream ss;//create a stringstream
    ss << ID;//add number to the stream

    cout << "Writing cytoscape graph files..." << endl;
    ofstream ofs(("cytArcs_" + ss.str() + ".txt").c_str());

    for (NodeID id = -numNodes; id <= numNodes; id++) {
        if (id == 0)
            continue;

        SSNode node = getSSNode(id);
        if (!node.isValid())
            continue;

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++)
            ofs << id << "\t" << it->getNodeID() << endl;
    }

    ofs.close();

    ofs.open(("cytColor" + ss.str() + ".txt").c_str());
    ofs << "Color_" << ID << " (class=Integer)" << endl;

    for (NodeID id = 1; id <= numNodes; id++) {

        DSNode& node = getDSNode(id);
        if (!node.isValid())
            continue;

        int value = trueMult[id] >= 1 ? 1 : 0;
        ofs << id << "=" << value << endl;
        ofs << -id << "=" << value << endl;
    }

    ofs.close();

    ofs.open(("cytNodes_seq" + ss.str() + ".txt").c_str());
    ofs << "Sequence_" << ID << " (class=String)" << endl;

    for (NodeID id = 1; id <= numNodes; id++) {
        DSNode& node = getDSNode(id);
        if (!node.isValid())
            continue;

        ofs << id << "=" << getSSNode(id).getSequence() << endl;
        ofs << -id << "=" << getSSNode(-id).getSequence() << endl;
    }

    ofs.close();

    ofs.open(("cytNodes_seqSize" + ss.str() + ".txt").c_str());
    ofs << "NodeLength_" << ID << " (class=Integer)" << endl;

    for (NodeID id = 1; id <= numNodes; id++) {
        DSNode& node = getDSNode(id);
        if (!node.isValid())
            continue;

        ofs << id << "=" << getSSNode(id).getMarginalLength() << endl;
        ofs << -id << "=" << getSSNode(-id).getMarginalLength() << endl;
    }

    ofs.close();*/
    //comment by Mahdi
    //for (srcID = 1; srcID <= numNodes; srcID++)
    //        if (getSSNode(srcID).isValid())
    //                break;
    int srcID=1;
    int i=0;
    for (i =srcID ; i <= numNodes; i++)
        if (getSSNode(i).isValid())
            break;
    //const size_t maxDepth =150;    // maximium sequence depth (exclusive srcID length)
    srcID=i;
    stringstream ss, ss2;   // create a stringstream
    ss << ID;               // add number to the stream
    ss2 << srcID;
    cout << "Writing cytoscape graph files around node " << srcID << endl;
    ofstream ofs(("cytArcs_" + ss2.str() + "_" +ss.str() + ".txt").c_str());
    ofs << "Source node\tTarget node\tArc coverage"<< endl;

    multimap<size_t, NodeID> nodeDepth;     // map of all nodes in the local graph and its depth
    nodeDepth.insert(pair<size_t, NodeID>(0, srcID));
    set<NodeID> nodesHandled;               // nodes that were already handled

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

        // if we're too far in the graph, stop
        //if (thisDepth > maxDepth)
        //    continue;
        if (thisID== 10633) {
            int stop=0;
            stop++;
        }
        SSNode thisNode = getSSNode(thisID);
        for (ArcIt it = thisNode.rightBegin(); it != thisNode.rightEnd(); it++) {
            SSNode rNode = getSSNode(it->getNodeID());
            if (!rNode.isValid())
                continue;
            if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                continue;
            ofs << thisID << "\t" << it->getNodeID() << "\t" << it->getCoverage() << endl;
            nodeDepth.insert(pair<size_t, NodeID>(thisDepth + thisNode.getMarginalLength(), it->getNodeID()));
        }

        for (ArcIt it = thisNode.leftBegin(); it != thisNode.leftEnd(); it++) {
            SSNode lNode = getSSNode(it->getNodeID());
            if (!lNode.isValid())
                continue;
            if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                continue;
            ofs << it->getNodeID() << "\t" << thisID << "\t" << it->getCoverage() << endl;
            nodeDepth.insert(pair<size_t, NodeID>(thisDepth + thisNode.getMarginalLength(), it->getNodeID()));
        }
    }

    ofs.close();
//comment by mahdi , adding new parameter in output file
    ofs.open(("cytNodes_" + ss2.str() + "_" +ss.str() + ".txt").c_str());
    ofs << "Node ID\tTrue multiplicity\tKmer coverage\tMarginal length\tReadStCov\t Multiplicity\t MulCertaintyRatio\t correctnessRatio\t tSequence" << endl;
    for (set<NodeID>::iterator it = nodesHandled.begin(); it != nodesHandled.end(); it++) {
        SSNode node = getSSNode(*it);
        double confidenceRatio=0;
        double nodeMultiplicity=0;
        double correctnessRatio=0;
        if(nodesExpMult.size()>1) {
            pair<int, pair<double,double> >  result=nodesExpMult[abs( node.getNodeID())];
            nodeMultiplicity=result.first;
            confidenceRatio=result.second.first;
            correctnessRatio=result.second.second;
        }
        ofs << *it << "\t" << trueMult[abs(*it)] << "\t" << node.getNodeKmerCov() << "\t"
            << node.getMarginalLength() << "\t"<<double(node.getReadStartCov()/node.getMarginalLength()) <<"\t"<<nodeMultiplicity<<"\t"<< confidenceRatio<<"\t"<< correctnessRatio <<"\t"<<node.getSequence() << endl;
    }
    ofs.close();
}

void DBGraph::writeLocalCytoscapeGraph(int ID, NodeID srcID, size_t maxDepth)
{

    //comment by Mahdi
    //for (srcID = 1; srcID <= numNodes; srcID++)
    //        if (getSSNode(srcID).isValid())
    //                break;

    int i=0;
    for (i =srcID ; i <= numNodes; i++)
        if (getSSNode(i).isValid())
            break;
    srcID=i;
    //const size_t maxDepth =150;    // maximium sequence depth (exclusive srcID length)

    stringstream ss, ss2;   // create a stringstream
    ss << ID;               // add number to the stream
    ss2 << srcID;
    cout << "Writing cytoscape graph files around node " << srcID << endl;
    ofstream ofs(("cytArcs_" + ss2.str() + "_" +ss.str() + ".txt").c_str());
    ofs << "Source node\tTarget node\tArc coverage"<< endl;

    multimap<size_t, NodeID> nodeDepth;     // map of all nodes in the local graph and its depth
    nodeDepth.insert(pair<size_t, NodeID>(0, srcID));
    set<NodeID> nodesHandled;               // nodes that were already handled

    while (!nodeDepth.empty()) {
        // get and erase the current node
        multimap<size_t, NodeID>::iterator e = nodeDepth.begin();
        size_t thisDepth = e->first;
        NodeID thisID = e->second;
        nodeDepth.erase(e);

        // if the node was already handled, skip
        if (nodesHandled.find(thisID) != nodesHandled.end())
            continue;

        // mark this node as handled
        nodesHandled.insert(thisID);

        // if we're too far in the graph, stop
        if (thisDepth > maxDepth)
            continue;

        SSNode thisNode = getSSNode(thisID);
        for (ArcIt it = thisNode.rightBegin(); it != thisNode.rightEnd(); it++) {
            SSNode rNode = getSSNode(it->getNodeID());
            if (!rNode.isValid())
                continue;
            if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                continue;
            ofs << thisID << "\t" << it->getNodeID() << "\t" << it->getCoverage() << endl;
            nodeDepth.insert(pair<size_t, NodeID>(thisDepth + thisNode.getMarginalLength(), it->getNodeID()));
        }

        for (ArcIt it = thisNode.leftBegin(); it != thisNode.leftEnd(); it++) {
            SSNode lNode = getSSNode(it->getNodeID());
            if (!lNode.isValid())
                continue;
            if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                continue;
            ofs << it->getNodeID() << "\t" << thisID << "\t" << it->getCoverage() << endl;
            nodeDepth.insert(pair<size_t, NodeID>(thisDepth + thisNode.getMarginalLength(), it->getNodeID()));
        }
    }

    ofs.close();
//comment by mahdi , adding new parameter in output file
    ofs.open(("cytNodes_" + ss2.str() + "_" +ss.str() + ".txt").c_str());
    ofs << "Node ID\tTrue multiplicity\tKmer coverage\tMarginal length\tReadStCov\t Multiplicity\t MulCertaintyRatio\t correctnessRatio\t tSequence" << endl;
    for (set<NodeID>::iterator it = nodesHandled.begin(); it != nodesHandled.end(); it++) {
        SSNode node = getSSNode(*it);
        double confidenceRatio=0;
        double nodeMultiplicity=0;
        double correctnessRatio=0;
        if(nodesExpMult.size()>1) {
            pair<int, pair<double,double> >  result=nodesExpMult[abs( node.getNodeID())];
            nodeMultiplicity=result.first;
            confidenceRatio=result.second.first;
            correctnessRatio=result.second.second;
        }
        ofs << *it << "\t" << trueMult[abs(*it)] << "\t" << node.getNodeKmerCov() << "\t"
            << node.getMarginalLength() << "\t"<<double(node.getReadStartCov()/node.getMarginalLength()) <<"\t"<<nodeMultiplicity<<"\t"<< confidenceRatio<<"\t"<< correctnessRatio <<"\t"<<node.getSequence() << endl;
    }
    ofs.close();
}

void DBGraph::readReferenceGenome()
{
    if (!reference.empty())
        return;

    cout << "Reading reference genome..." << endl;

    FastAFile ass(false);
//
    ass.open( "genome.fasta");
    string read, desc;
    while (ass.getNextRead(read, desc)) {
        read.append(read);
        reference.push_back(read);
        refST.push_back(ST_CreateTree(read.c_str(), read.size()));
        cout << "Adding a reference of length: " << read.size() / 2 << endl;
    }

    cout << "Done reading " << reference.size() << " reference contigs." << endl;
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

void DBGraph::compareToSolution()
{
    readReferenceGenome();

    ofstream ofs("Staphmult.txt");
    ofs << "Node ID\tMarginal length\tMultiplicity" << endl;
    size_t totalSizeOfIncNodes=0;
    size_t totalSizeOfCNodes=0;
    size_t totalSizeOfNodes=0;
    size_t numOK = 0, numNotOK = 0, numTotal = 0, totalCov = 0;
    for (NodeID id = 1; id <= numNodes; id++) {
        string P = getSSNode(id).getSequence();

        if (!getDSNode(id).isValid())
            continue;

        numTotal++;

        vector<vector<size_t> > pos, posRC;
        int ntOcc = findAllTrueOccurences(P, pos, posRC);

        trueMult[id] = ntOcc;

        if (ntOcc == 0) {
            getDSNode(id).setFlag(0);
            numNotOK++;
            totalSizeOfIncNodes=totalSizeOfIncNodes+getSSNode(id).getMarginalLength();

        } else {
            getDSNode(id).setFlag(1);
            numOK++;
            totalSizeOfCNodes=totalSizeOfCNodes+getSSNode(id).getMarginalLength();
        }
        totalSizeOfNodes=totalSizeOfNodes+getSSNode(id).getMarginalLength();

        ofs << id << "\t" << getSSNode(id).getMarginalLength() << "\t" << trueMult[id] << "\t" << endl;

        if (getDSNode(id).getNumRightArcs() == 0)
            totalCov += trueMult[id] * getDSNode(id).getLength();
        else
            totalCov += trueMult[id] * getDSNode(id).getMarginalLength();
    }

    ofs.close();

    size_t totalGenomeSize = 0;
    for (size_t i = 0; i < reference.size(); i++)
        totalGenomeSize += reference[i].size()/2;

    cout << "\tNumber of nodes that exist: " << numOK << " (" << 100.00 * numOK / numTotal << "%)" <<     "        ===>> "<<" ("<<100.00 * totalSizeOfCNodes / totalSizeOfNodes<<"%)"<<" size of graph" <<endl;
    cout << "\tNumber of nodes that do not exist: " << numNotOK << " (" << 100.00 * numNotOK / numTotal << "%)" << "===>> "<<" ("<<100.00 * totalSizeOfIncNodes / totalSizeOfNodes<<"%)"<<" size of graph" <<endl;
    cout << "\tThe fraction of the genome that is covered CHECK : " << 100.00 * totalCov / totalGenomeSize << "%" << endl;

    size_t TP = 0, FP = 0, TN = 0, FN = 0;
    for (NodeID id = 1; id <= numNodes; id++) {
        DSNode &node = getDSNode(id);

        if (!node.isValid())
            continue;

        /*cout << node.getMarginalLength() << "\t"
             << (int)node.getMultiplicity() << "\t"
             << (int)node.getFlag() << endl;*/

        if (node.getFlag() == 1) {
            if (node.getRoundMult() >= 1)
                TP++;
            else
                FN++;
        } else {
            if (node.getRoundMult() >= 1)
                FP++;
            else
                TN++;
        }
    }

    //mahdi comment
   /* table = new KmerNodeTable ( settings, numNodes );
    table->populateTable ( nodes );
    FastAFile ass(false);
    ass.open(settings.inputDirectory+"genome.fasta");
    size_t numKmers = 0, numFound = 0, numOL = 0;
    string read, desc;
    while (ass.getNextRead(read, desc)) {
        TString Tread=read;
        for ( TStringIt it = Tread.begin(); it != Tread.end(); it++ ) {
            Kmer kmer = *it;
            NodePosPair result = table->find ( kmer );
            if ( result.isValid() ) {
                numFound++;
            }
            numKmers++;
        }
    }*/



    /////////////////////////////////////////////////////////////////


   /* string readFileName="sampleRead.fastq";
    ifstream readsFile;
    string reads100=settings.inputDirectory+ readFileName;
    readsFile.open(reads100.c_str(),ios::in);

    readStructStr readInfo;
    int numOfReads=0;
    int numOfFoundRead=0;
    double avgOfMaxConse=0;
    while(getline(readsFile,readInfo.strID )&&getline(readsFile,readInfo.erroneousReadContent )) {
        getline(readsFile,readInfo.orientation );
        getline(readsFile,readInfo.qualityProfile );
        numOfReads++;

        bool again=true;
        int kmerStart=0;
        bool supported=false;

        int startOfRead=1;
        int startOfNode=0;
        TString read =readInfo.erroneousReadContent;

        if (readInfo.erroneousReadContent.length()<kmerSize)
            continue;
        bool found=false;
        int numOfSuppKmerInRead=0;
        int prevKmerID=0;
        int numOfConsecutive=0;
        int maxConsecutive=0;

        for ( TStringIt it = read.begin(); it != read.end(); it++ ) {
            Kmer kmer = *it;
            NodePosPair result = table->find ( kmer );
            if ( !result.isValid() ) {
                startOfRead++;
                continue;
            }
            else {
                if (startOfRead-prevKmerID==1) {
                    numOfConsecutive++;
                    if (numOfConsecutive>maxConsecutive) {
                        maxConsecutive=numOfConsecutive;
                    }
                } else {
                    numOfConsecutive=0;
                }
                prevKmerID=startOfRead;
                startOfRead++;
                found=true;
            }

        }
        if ( found) {
	    avgOfMaxConse=avgOfMaxConse+maxConsecutive;
            numOfFoundRead++;
        }


    }
    avgOfMaxConse=avgOfMaxConse/(float)numOfFoundRead;*/
    ////////////////////////////////////////////////////////////
    //delete table;
    cout << "Validation report: " << endl;
    //cout<<"avarage of Consecutive kmerHit in sample read: "<< avgOfMaxConse<<endl;
   // cout << "\tfraction of reads which are supported by at least one kmer in Graph: " << numOfFoundRead << "/" << numOfReads << "(" << 100.00*(double)numOfFoundRead/numOfReads << "\%)" << endl;
   // cout << "\tfraction of kmers in genome which is still exist in Graph: " << numFound << "/" << numKmers << "(" << 100.00*(double)numFound/numKmers << "\%)" << endl;

    cout << "\tTrue positives: " <<  TP << " (" << 100.00 * TP / numTotal << "%)" << endl
         << "\tFalse positives: " << FP << " (" << 100.00 * FP / numTotal << "%)" << endl
         << "\tTrue negatives: " << TN << " (" << 100.00 * TN / numTotal << "%)" << endl
         << "\tFalse negatives: " << FN << " (" << 100.00 * FN / numTotal << "%)" << endl;

    cout << "========================================================" << endl;
}

size_t DBGraph::findAllTrueOccurences(const string& P,
                                      vector<vector<size_t> >& positions,
                                      vector<vector<size_t> >& positionsRC) const
{
    positions.clear();
    positionsRC.clear();
    size_t ntOcc = 0;

    // normal strand
    for (size_t i = 0; i < refST.size(); i++) {
        positions.push_back(vector<size_t>());
        SUFFIX_TREE *st = (SUFFIX_TREE*)refST[i];

        int *pos;
        char *c_str = const_cast<char*>(P.c_str());
        int numOcc = ST_FindAllSubstrings(st, c_str, P.size(), &pos);

        for (int j = 0; j < numOcc; j++) {
            if ((pos[j]-1) >= (int)reference[i].size()/2)
                continue;
            positions.back().push_back(pos[j]-1);
            ntOcc++;
        }

        if (numOcc > 0)
            free(pos);
    }

    // calculate the reverse complement, if it is the same, get out now
    string RC = Nucleotide::getRevCompl(P);

    if (RC == P)
        return ntOcc;

    // opposite strand
    for (size_t i = 0; i < refST.size(); i++) {
        positionsRC.push_back(vector<size_t>());
        SUFFIX_TREE *st = (SUFFIX_TREE*)refST[i];

        int *pos;
        char *c_str = const_cast<char*>(RC.c_str());
        int numOcc = ST_FindAllSubstrings(st, c_str, RC.size(), &pos);

        for (int j = 0; j < numOcc; j++) {
            if ((pos[j]-1) >= (int)reference[i].size()/2)
                continue;
            positionsRC.back().push_back(pos[j]-1);
            ntOcc++;
        }

        if (numOcc > 0)
            free(pos);
    }

    return ntOcc;
}

void DBGraph::getActualPath(NodeID srcID, NodeID dstID, size_t distance,
                            vector<NodeID>& result) const
{
    result.clear();
    result.push_back(srcID);
    SSNode curr = getSSNode(srcID);

    vector<vector<size_t> > pos, posRC;
    findAllTrueOccurences(curr.getSequence(), pos, posRC);

    for (size_t i = 0; i < pos.size(); i++) {
        if (pos[i].empty())
            continue;

        size_t currPos = pos[i].front() + curr.getMarginalLength();
        size_t currLen = curr.getMarginalLength();

        while (currLen <= distance) {

            NodeID left = curr.getRightArc(reference[i][currPos+Kmer::getK()-1]);
            if (left == 0)
                return;

            SSNode next = getSSNode(left);
            result.push_back(next.getNodeID());
            currPos += next.getMarginalLength();
            currLen += next.getMarginalLength();
            curr = next;

            if (curr.getNodeID() == dstID)
                break;
        }

        return;
    }

    for (size_t i = 0; i < posRC.size(); i++) {
        if (posRC[i].empty())
            continue;

        size_t currPos = posRC[i].front()+Kmer::getK()-1;
        size_t currLen = curr.getMarginalLength();

        while (currLen <= distance) {

            NodeID right = curr.getRightArc(Nucleotide::getComplement(reference[i][currPos-Kmer::getK()]));
            if (right == 0)
                return;

            SSNode next = getSSNode(right);
            result.push_back(next.getNodeID());
            currPos -= next.getMarginalLength();
            currLen += next.getMarginalLength();
            curr = next;

            if (curr.getNodeID() == dstID)
                break;
        }

        return;
    }
}

bool DBGraph::getDistance(NodeID srcID, NodeID dstID, double maxDistance,
                          double &result) const
{
    vector<NodeID> realpath;
    getActualPath(srcID, dstID, maxDistance, realpath);

    double length = 0;
    for (size_t i = 0; i < realpath.size(); i++) {
        if (realpath[i] == dstID) {
            result = length;
            return true;
        }
        length += getSSNode(realpath[i]).getMarginalLength();
    }

    realpath.clear();
    maxDistance = maxDistance + getSSNode(srcID).getMarginalLength() - getSSNode(dstID).getMarginalLength();
    if (maxDistance < 0)
        return false;
    getActualPath(-srcID, -dstID, maxDistance, realpath);

    length = 0;
    for (size_t i = 0; i < realpath.size(); i++) {
        if (realpath[i] == -dstID) {
            result = -(length - getSSNode(srcID).getMarginalLength() + getSSNode(dstID).getMarginalLength());
            return true;
        }
        length += getSSNode(realpath[i]).getMarginalLength();
    }

    return false;
}

void DBGraph::getTrueOccurencesForReadPair(string &start,
        string &stop,
        vector<string>& output,
        size_t maxSize)
{
    vector<vector<size_t> > pos1, posRC1, pos2, posRC2;
    findAllTrueOccurences(start, pos1, posRC1);
    findAllTrueOccurences(stop, pos2, posRC2);

    // forward strand
    for (size_t i = 0; i < pos1.size(); i++)
        for (size_t j = 0; j < pos1[i].size(); j++)
            for (size_t k = 0; k < pos2[i].size(); k++) {
                if (pos2[i][k] <= pos1[i][j])
                    continue;
                if (pos2[i][k] - pos1[i][j]  + stop.size() <= maxSize)
                    output.push_back(reference[i].substr(pos1[i][j], pos2[i][k] - pos1[i][j] + stop.size()));
            }

    // RC strand
    for (size_t i = 0; i < posRC1.size(); i++)
        for (size_t j = 0; j < posRC2[i].size(); j++)
            for (size_t k = 0; k < posRC1[i].size(); k++) {
                if (posRC1[i][k] <= posRC2[i][j])
                    continue;
                if (posRC1[i][k] - posRC2[i][j] + start.size() <= maxSize)
                    output.push_back(Nucleotide::getRevCompl(reference[i].substr(posRC2[i][j], posRC1[i][k] - posRC2[i][j] + start.size())));
            }
}


