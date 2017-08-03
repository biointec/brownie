#include "findGap.h"
#include "alignment.h"
#include <string>
#include <map>
#include <queue>
#include "library.h"

using namespace std;


class DFSNode
{
public:
        NodeID nodeID;
        size_t readPos;
        int score;
        float relScore;


        DFSNode() : nodeID(0), readPos(0), score(0), relScore(0.0f) {}

        DFSNode(NodeID nodeID_, size_t readPos_, int score_, float relScore_) :
                nodeID(nodeID_), readPos(readPos_), score(score_),
                relScore(relScore_) {}

        bool operator< (const DFSNode& rhs) const {
                if (rhs.relScore != relScore)
                        return rhs.relScore < relScore;
                return rhs.score < score;
        }
};


void FindGap::parameterInitialization(){
        maxSearchSize = 50;
        minComponentSize = 4;
        maxComponentSize = 990000;
        minNumbOfPairs = 10;
        minOverlapSize = 15;
        minSim =50;
        minExactMatchSize = 7;
}

FindGap::FindGap(string nodeFileName, string arcFileName, string metaDataFileName, string alignmentFile,unsigned int kmerVlue  ,string tempDir):settings(kmerVlue,tempDir), dbg(settings),overlapSize(kmerVlue),alignment(1000, 2, 1, -1, -3)
{
        Kmer::setWordSize(settings.getK());
        cout <<"graph is loading ...\n" ;
        dbg.loadGraph(nodeFileName,arcFileName,metaDataFileName);
        cout << dbg.getGraphStats() << endl;
        kmerSize = kmerVlue;
        correctedFile = alignmentFile;
        parameterInitialization();
}


FindGap::FindGap(LibraryContainer& libraries, const Settings& s, DBGraph &graph) :overlapSize(s.getK()),alignment(1000, 2, 1, -1, -3),dbg(graph)
{
        //call it like this :
        //FindGap  findGap (libraries, settings, graph);
        //findGap.closeGaps();
        
        Kmer::setWordSize(settings.getK());
        settings =s;
        kmerSize = settings.getK();
        parameterInitialization();
        ReadLibrary &input =libraries.getInput(0);
        ReadFile *readFile = input.allocateReadFile();
        correctedFile = input.getInputFilename();
}


void FindGap::connectComponents(string nodeFilename, string arcFilename,string metaDataFilename)
{       size_t numberOfJoins = 0 ;
        //dbg.writeCytoscapeGraph("cytoscape");
        Util::startChrono();
        cout << "Creating kmer lookup table... "; cout.flush();
        dbg.buildKmerNPPTable();
        cout << "done (" << Util::stopChronoStr() << ")" << endl;
        Util::startChrono();
        std::map<int, int > nodeComponentMap;
        std::map< pair<int, int>, int> pairedEndJoins;
        ComponentHandler componentHdl(dbg,settings);
        findComponentsInGraph( componentHdl);
        extractPairedComponents(correctedFile , componentHdl, pairedEndJoins);
        cout <<"There are " <<pairedEndJoins.size() <<"number of possible joins between components" <<endl;
        for(auto it = pairedEndJoins.cbegin(); it != pairedEndJoins.cend(); ++it)
        {
                set<int> tipNodes;
                Component c1 = componentHdl.componentsMap[it->first.first];
                Component c2 = componentHdl.componentsMap[it->first.second];
                set <int> nodeSet = c1.getNodeIDSet();

                for (auto nodeID: c2.getNodeIDSet())
                        nodeSet.insert(nodeID);
                findTips(tipNodes , nodeSet);
                //extrac the kmers in the tips and store them in  map which tells you that kmer exist in which nodes (here tips)
                overlapSize = kmerSize /2;
                if (overlapSize%2 ==0)
                        overlapSize = overlapSize +1;
                cout <<"============================================================" <<endl;
                cout << "Finding overlap between component " <<c1.componentID <<" and component " <<c2.componentID <<endl;

                bool overlapDetected = false;
                while (!overlapDetected && overlapSize >= 7){
                        cout << "The overlapSize is set to "<<overlapSize <<endl;
                        std::map<string, set<int> > kmerNodeMap;
                        std::vector<pair<int, int> > potentialJoin;
                        loadKmerMap(tipNodes,kmerNodeMap,overlapSize);
                        getPotentialJoins(componentHdl,kmerNodeMap,potentialJoin);
                        overlapDetected = eliminateFakeJoins(potentialJoin, pairedEndJoins, componentHdl);
                        if (potentialJoin.size() >0)
                                break;
                        overlapSize = overlapSize / 2;
                        if (overlapSize%2 ==0)
                                overlapSize = overlapSize +1;
                        kmerNodeMap.clear();
                }
                if (overlapDetected) numberOfJoins++;
        }
        dbg.concatenateNodes();
        dbg.writeGraph(nodeFilename,arcFilename,metaDataFilename);
        cout << "Number of new links in the graph " <<numberOfJoins << " from " <<pairedEndJoins.size()<<" initial suggestions" <<endl;
}

void FindGap::closeGaps(string nodeFilename, string arcFilename,string metaDataFilename)
{       //dbg.writeCytoscapeGraph("cytoscape");
        ofstream ofs;
        ofs.open("pairedArtificialReads.fasta",std::fstream::out);
        ofs.close();

        size_t numberOfJoins = 0 ;
        set<int> tipNodes;
        vector< pair< pair<int , int> , int > > potentialPairs;
        findPotentialTips(tipNodes);
        cout << "Creating kmer lookup table... "; cout.flush();
        dbg.buildKmerNPPTable();
        streamReads(correctedFile,tipNodes,potentialPairs);
        checkForTipConnection(potentialPairs);
        dbg.concatenateNodes();
        cout << dbg.getGraphStats() << endl;
        if (nodeFilename != "")
                dbg.writeGraph(nodeFilename,arcFilename,metaDataFilename);
}


void FindGap::reorderTips(SSNode &first, SSNode &second)
{
        if (first.getNumRightArcs() !=0 && first.getNumLeftArcs() ==0)
                first = dbg.getSSNode( -first.getNodeID());
        if (second.getNumLeftArcs()!=0 && second.getNumRightArcs() ==0)
                second = dbg.getSSNode( -second.getNodeID());

        // the first tip is a single node
        if (first.getNumRightArcs() ==0 && first.getNumLeftArcs() ==0 && second.getNumRightArcs() !=0){
                string firstRead1="", secondRead1="";
                double sim1 = -100,sim2= -100;
                if (alignTips(first.getNodeID(), second.getNodeID(), firstRead1, secondRead1))
                        sim1 = alignment.align(firstRead1, secondRead1)*100/(int)firstRead1.length();
                string firstRead2="", secondRead2="";
                if (alignTips(-first.getNodeID(), second.getNodeID(), firstRead2, secondRead2))
                        sim2 =alignment.align(firstRead2, secondRead2)*100/(int)firstRead2.length();
                if (sim2 > sim1 )
                        first = dbg.getSSNode( -first.getNodeID());


        }
        //the second tip is a single node
        if (second.getNumLeftArcs()==0 && second.getNumRightArcs() ==0 && first.getNumLeftArcs()!=0){
                string firstRead1="", secondRead1="";
                double sim1 = -100,sim2= -100;
                if (alignTips(first.getNodeID(), second.getNodeID(), firstRead1, secondRead1))
                        sim1 = alignment.align(firstRead1, secondRead1)*100/(int)firstRead1.length();
                string firstRead2="", secondRead2="";
                if (alignTips(first.getNodeID(), -second.getNodeID(), firstRead2, secondRead2))
                        sim2 = alignment.align(firstRead2, secondRead2)*100/(int)firstRead2.length();
                if (sim2 > sim1  )
                        second = dbg.getSSNode( -second.getNodeID());


        }
        // both tips are single nodes
        if (second.getNumLeftArcs()==0 && second.getNumRightArcs() ==0 && first.getNumLeftArcs()==0 && first.getNumRightArcs()==0){
                string firstRead1="", secondRead1="";
                double sim1 = -100,sim2 = -100,sim3 = -100,sim4 = -100;
                if (alignTips(first.getNodeID(), second.getNodeID(), firstRead1, secondRead1))
                        sim1 = alignment.align(firstRead1, secondRead1)*100/(int)firstRead1.length();
                string firstRead2="", secondRead2="";
                if (alignTips(-first.getNodeID(), second.getNodeID(), firstRead2, secondRead2))
                        sim2 = alignment.align(firstRead2, secondRead2)*100/(int)firstRead2.length();
                string firstRead3="", secondRead3="";
                if (alignTips(first.getNodeID(), -second.getNodeID(), firstRead3, secondRead3))
                        sim3 = alignment.align(firstRead3, secondRead3)*100/(int)firstRead3.length();
                string firstRead4="", secondRead4="";
                if (alignTips(-first.getNodeID(), -second.getNodeID(), firstRead4, secondRead4))
                        sim4 = alignment.align(firstRead4, secondRead4)*100/(int)firstRead4.length();
                if (sim2 > sim1 && sim2 > sim3 && sim2 >sim4){
                        first = dbg.getSSNode( -first.getNodeID());
                }else
                        if (sim3 > sim1 && sim3 >sim2 && sim3 > sim4){
                                second = dbg.getSSNode( -second.getNodeID());
                        }else
                                if (sim4 >sim1 && sim4>sim2 && sim4>sim3  )
                                        first = dbg.getSSNode( -first.getNodeID());
                                second = dbg.getSSNode( -second.getNodeID());

        }
}


void FindGap::checkForTipConnection(vector< pair< pair<int , int> , int > >& potentialPairs)
{
        ofstream ofs;
        ofs.open("pairedArtificialReads.fasta",std::fstream::out);
        size_t violateGraphStructure =0;
        cout << "Make connection between tips which alreay have been found ..." <<endl;
        size_t numberOfNewConnection = 0;
        set<int> nodesHandled;
        for(auto it = potentialPairs.cbegin(); it != potentialPairs.cend(); ++it)
        {
                bool found = false;
                SSNode first = dbg.getSSNode( it->first.first);
                SSNode second = dbg.getSSNode( it->first.second);
                //these nodes have been already connected in other combinations
                if (nodesHandled.find(abs(first.getNodeID())) !=nodesHandled.end())
                        continue;
                if (nodesHandled.find(abs(second.getNodeID())) != nodesHandled.end())
                        continue;
                if ( !first.isValid() || !second.isValid() )
                        continue;
                reorderTips(first, second);

                //#ifdef DEBUG
                string firstNodeContent = first.getSequence();
                string secondNodeContent = second.getSequence();
                size_t maxLen = firstNodeContent.length()  < 150 ? firstNodeContent.length():150;
                firstNodeContent = firstNodeContent.substr(firstNodeContent.length()-maxLen,maxLen);
                maxLen = secondNodeContent.length()  < 150 ? secondNodeContent.length():150;
                secondNodeContent = secondNodeContent.substr(0,maxLen);
                //#endif
                string firstRead="", secondRead="";
                double sim = -100;
                if (alignTips(first.getNodeID(), second.getNodeID(), firstRead, secondRead)){
                        sim = alignment.align(firstRead, secondRead)*100/(int)firstRead.length();
                        if (firstRead.length() >minOverlapSize and sim>minSim ){
                                cout <<first.getNodeID() << " : " <<second.getNodeID() <<endl;
                                alignment.printAlignment(firstRead,secondRead);
                                cout <<"similarity:"<<sim <<endl;
                                if (first.getAvgKmerCov() >= second.getAvgKmerCov())
                                        found = connectNodes(first.getNodeID(),second.getNodeID());
                                else
                                        found = connectNodes(-second.getNodeID(), -first.getNodeID());
                                if (found)
                                {
                                        nodesHandled.insert(abs(first.getNodeID()));
                                        nodesHandled.insert(abs(second.getNodeID()));
                                        numberOfNewConnection ++;
                                        cout <<"The connection successfully stablished. This connection is supported by " <<it->second <<endl;
                                }
                                else
                                {
                                        violateGraphStructure ++;
                                        cout <<"We couldn't make this connection. It probably violates the graph structure" <<endl;
                                }

                                cout <<".............................."<<endl;
                        }

                }
                //#ifdef DEBUG
                if (firstRead.length() > firstNodeContent.length())
                        firstRead = firstNodeContent ;
                if (secondRead.length() >secondNodeContent.length())
                        secondNodeContent = secondRead;
                ofs << ">" <<first.getNodeID() <<"_"<<found <<"_" <<sim <<".1\n" << firstNodeContent  <<endl;
                ofs << ">" <<second.getNodeID()<<"_"<<found <<"_" <<sim <<".2\n" << secondNodeContent <<endl;
                //#endif


        }
        cout <<violateGraphStructure << " number of connection were skipped because they make inconsistency in the graph structure." <<endl;
        cout << numberOfNewConnection <<" new connection were stablished\n" <<endl;
        ofs.close();
}


void FindGap::findNPPFast(const string& read, vector<NodePosPair>& nppv)
{
        for (KmerIt it(read); it.isValid(); it++) {
                Kmer kmer = it.getKmer();
                NodePosPair npp = dbg.findNPP(kmer);
                nppv[it.getOffset()] = npp;

                if (!npp.isValid())
                        continue;

                NodeID nodeID = npp.getNodeID();
                const SSNode node = dbg.getSSNode(nodeID);

                size_t readPos = it.getOffset() + Kmer::getK();
                size_t nodePos = npp.getPosition() + Kmer::getK();

                while ((readPos < read.size()) && (nodePos < node.getLength())) {
                        if (read[readPos] == node.getNucleotide(nodePos)) {
                                it++;
                                nppv[it.getOffset()] = NodePosPair(nodeID, nodePos - Kmer::getK() + 1);
                                nodePos++;
                                readPos++;
                        } else
                                break;

                }
        }
}

void FindGap::extractPairedComponents(string readFileName , ComponentHandler& componentHdl, std::map< pair<int, int>, int>& pairedEndJoins)
{

        // it is assumed that the file in fastq format
        // insted of a read, chain of nodes which that read mapped to is written.
        std::map< pair<int, int>, int> pairedEndJoinsTemp;
        cout <<"Finding the possibility of joining components by paired-end reads"<<endl;
        std::ifstream input(readFileName);
        vector< pair<string, string> > breakpoints;
        if (!input.good()) {
                std::cerr << "Error opening: " << readFileName << std::endl;
                return;
        }
        std::string line, Nodes ="" ;
        int i = 0;
        vector<int> firstPairNodes;
        vector<int> secondPairNodes;
        while (std::getline(input, line).good()) {
                if (i % OUTPUT_FREQUENCY == 0) {
                        cout << "\tProcessing read " << (i/4) << "\r";
                        cout.flush();
                }

                if ( i%8 == 1 ){
                        splitLineToNodes(line , firstPairNodes);
                }
                if ( i%8 == 5 ){

                        splitLineToNodes(line , secondPairNodes);
                }
                if ( i%8 == 6){
                        bool found = false;
                        for (auto it1:firstPairNodes){
                                for (auto it2:secondPairNodes){
                                        if (componentHdl.getComponentID(it1) == 0 || componentHdl.getComponentID(it2) == 0)
                                                        continue;
                                        if (componentHdl.getComponentID(it1) != componentHdl.getComponentID(it2)){

                                                Component c1 = componentHdl.componentsMap[componentHdl.getComponentID(it1)];
                                                Component c2 = componentHdl.componentsMap[componentHdl.getComponentID(it2)];
                                                /*if (c1.numOfNodes ==1 || c2.numOfNodes==1)
                                                        continue;*/
                                                pair<int, int> join;
                                                if (c1.componentID < c2.componentID)
                                                        join = make_pair(c1.componentID,c2.componentID);
                                                else
                                                        join = make_pair(c2.componentID,c1.componentID);

                                                if ( pairedEndJoinsTemp.find(join) == pairedEndJoinsTemp.end() ){
                                                        pairedEndJoinsTemp[join] = 1;
                                                        found = true;
                                                }

                                                else{
                                                        pairedEndJoinsTemp[join] = pairedEndJoinsTemp[join] +1;
                                                        found =true;
                                                }
                                                if (found)
                                                        break;

                                        }
                                        if (found)
                                                break;
                                }
                        }

                        firstPairNodes.clear();
                        secondPairNodes.clear();
                }
                i ++ ;
        }


        for(auto it = pairedEndJoinsTemp.cbegin(); it != pairedEndJoinsTemp.cend(); ++it)
        {
               if ( it->second >minNumbOfPairs )
                       pairedEndJoins[it->first] =it->second ;
        }
        cout <<"\n Number of found suggestions: "<<pairedEndJoins.size() <<endl;
        for(auto it = pairedEndJoins.cbegin(); it != pairedEndJoins.cend(); ++it)
        {
                Component c1 = componentHdl.componentsMap[it->first.first];
                Component c2 = componentHdl.componentsMap[it->first.second];
                std::cout <<"component "  <<c1.componentID <<" (size = "<<c1.componentSize<<") and component " << c2.componentID << " (size = "<< c2.componentSize<<") supported by " << it->second << " pairs" <<"\n";
        }

}
void FindGap::splitLineToNodes(string str, vector<NodeID> &nodeIDs)
{
        string next;
        //const char delim = "\t";
        char delim = '\t';

        for (int i=0; i<str.length();i++) {
                char it = str[i];
                // If we've hit the terminal character
                if (it == delim) {
                        // If we have some characters accumulated
                        if (!next.empty()) {
                                // Add them to the result vector
                                nodeIDs.push_back(atoi( next.c_str()));
                                next.clear();
                        }
                } else {
                        // Accumulate the next character into the sequence
                        next += it;
                }
        }
}

bool FindGap::eliminateFakeJoins(std::vector<pair<int, int> > &potentialJoin, std::map< pair<int, int>, int> &pairedEndJoins ,  ComponentHandler& componentHdl)
{
        bool found = false;
        size_t bestSim=0;
        SSNode bestRight, bestLeft;
        string bestleftStr = "", bestrightStr = "";
        cout <<"number  of potential joins : " <<potentialJoin.size() <<endl;
        for (auto it:potentialJoin){
                string firstRead,secondRead = "";
                alignTips(it.first, it.second, firstRead, secondRead);
                double sim =alignment.align(firstRead, secondRead)*100/(int)firstRead.length();
                if (sim > bestSim && firstRead.length() >minOverlapSize ){
                        bestSim = sim;
                        bestLeft = dbg.getSSNode( it.first);
                        bestRight = dbg.getSSNode(it.second);
                        bestleftStr = firstRead;
                        bestrightStr = secondRead;
                }
        }

        if (bestSim > 50){
                cout <<"                next pair               \n=========================" <<endl;
                cout<<"first :"<<bestLeft.getNodeID()<<"\t" <<bestleftStr <<endl;
                cout<<"second:"<<bestRight.getNodeID()<<"\t" <<bestrightStr <<endl;
                cout <<"normalized similarity: " <<bestSim <<endl;
                alignment.align(bestleftStr, bestrightStr);
                alignment.printAlignment(bestleftStr, bestrightStr);
                cout<<"this suggests a connection is between components "<<componentHdl.getComponentID(bestLeft.getNodeID()) <<" and " <<componentHdl.getComponentID(bestRight.getNodeID())<<endl;
                if (componentHdl.getComponentID(bestLeft.getNodeID()) < componentHdl.getComponentID(bestRight.getNodeID())){
                        cout<<"this suggestion is supported by " << pairedEndJoins[make_pair(componentHdl.getComponentID(bestLeft.getNodeID()),componentHdl.getComponentID(bestRight.getNodeID()))] <<" pair end read mapp" <<endl;;
                }else{
                        cout<<"this suggestion is supported by " << pairedEndJoins[make_pair(componentHdl.getComponentID(bestRight.getNodeID()),componentHdl.getComponentID(bestLeft.getNodeID()))] <<" pair end read mapp" <<endl;;
                }
                found = true;
                if (bestLeft.getAvgKmerCov() >= bestRight.getAvgKmerCov())
                        found = connectNodes(bestLeft.getNodeID(),bestRight.getNodeID());
                else
                        found = connectNodes(-bestRight.getNodeID(), -bestLeft.getNodeID());
        }
        if (found)
                cout <<"The connection successfully stablished." <<endl;
        else
                cout <<"We couldn't make this connection." <<endl;

        return found;
}
bool FindGap::alignTips(int firstNodeId, int secondNodeId, string& firstRead, string &secondRead)
{
        SSNode first=dbg.getSSNode(firstNodeId);
        SSNode second=dbg.getSSNode(secondNodeId);
        size_t firsStartIndex =0 , secondStartIndex = 0, firstEndIndex = 0, secondEndIndex = 0;
        if (!extendRead(firsStartIndex,  secondStartIndex,  firstEndIndex ,secondEndIndex, first,second ,firstRead , secondRead))
                return false;
        if (second.getNumRightArcs() !=0 && first.getSequence().length()>firstEndIndex){
                //cout<<"Expand second node to the right for "<<  first.getSequence().length()-firstEndInex <<" bp "<<endl;
                expandReadByGraphToRight(second,first,secondEndIndex,firstEndIndex,secondRead,firstRead);
        }
        if (first.getNumLeftArcs()!=0 && secondStartIndex >0){
                //cout <<"Expand first node to the left for "<< secondStartIndex<<" bp "<<endl;
                expandReadByGraphToLeft(first,second,firsStartIndex,secondStartIndex,firstRead,secondRead);
        }
        return true;
}

bool FindGap::connectNodes(NodeID firstNodeID, NodeID secondNodeID)
{
        SSNode first = dbg.getSSNode(firstNodeID);
        SSNode second = dbg.getSSNode(secondNodeID);
        size_t firsStartIndex =0 , secondStartIndex = 0, firstEndIndex = 0, secondEndIndex = 0;
        string firstRead,secondRead = "";
        string consensusStr="";
        extendRead(firsStartIndex,  secondStartIndex,  firstEndIndex ,
                   secondEndIndex, first,second ,firstRead , secondRead);
        if (firsStartIndex >0)
                consensusStr = first.getSequence().substr(0,firsStartIndex);
        consensusStr = consensusStr +firstRead;
        if (second.getSequence().length()>secondEndIndex )
                consensusStr = consensusStr +second.getSequence().substr(secondEndIndex,second.getSequence().length()-secondEndIndex);
        if (first.getSequence().substr(0,kmerSize-1).compare( consensusStr.substr(0,kmerSize-1)) != 0 )
                return false;
        if (second.getSequence().substr(second.getSequence().length()-kmerSize+1, kmerSize-1).compare(consensusStr.substr(consensusStr.length()-kmerSize+1, kmerSize-1)) !=0 )
                return false;

        for (KmerIt it(consensusStr);it.isValid(); it++)
        {
                Kmer kmer = it.getKmer();
                NodePosPair result= dbg.findNPP(kmer);
                if (result.isValid()&& abs(result.first) != abs(firstNodeID) && abs(result.first) !=abs(secondNodeID))
                        return false;
        }
        first.deleteAllRightArcs();
        first.inheritRightArcs(second);
        first.setSequence(consensusStr);
        for (ArcIt it = second.leftBegin(); it != second.leftEnd(); it++) {
                SSNode leftNode = dbg.getSSNode(it->getNodeID());
                if (leftNode.getNodeID()==-second.getNodeID())
                        continue;
                bool result = leftNode.deleteRightArc (second.getNodeID());
                assert(result);
        }
        for (ArcIt it = second.rightBegin(); it != second.rightEnd(); it++) {
                SSNode rightNode = dbg.getSSNode ( it->getNodeID() );
                if (rightNode.getNodeID()==-second.getNodeID())
                        continue;
                bool result = rightNode.deleteLeftArc(second.getNodeID());
                assert(result);
        }
        second.deleteAllLeftArcs();
        second.deleteAllRightArcs();
        second.invalidate();
        return true;


}



void FindGap::streamReads(string readFileName , set<int> &tipNodes,  vector< pair< pair<int , int> , int > >& potentialPairs )
{
        cout <<"\nStreaming reads to find potential connections between tips" <<endl;
        std::map< pair<int, int>, int> pairedEndJoinsTemp;
        std::ifstream input(readFileName);
        vector< pair<string, string> > breakpoints;
        if (!input.good()) {
                std::cerr << "Error opening: " << readFileName << std::endl;
                return;
        }
        std::string line, Nodes ="" ;
        int i = 0;
        string firstPair, secondPair;
        vector<int> secondPairNodes;
        while (std::getline(input, line).good()) {
                if (i % OUTPUT_FREQUENCY == 0) {
                        cout << "\t Processing read " << (i/4) << "\r";
                        cout.flush();
                }
                if ( i%8 == 1 )
                        firstPair = line;
                if ( i%8 == 5 )
                        secondPair = line;
                if ( i%8 == 6){
                        vector<NodePosPair> nppvFirst(firstPair.length() + 1 - Kmer::getK());
                        findNPPFast(firstPair, nppvFirst);
                        vector<NodePosPair> nppvSecond(secondPair.length() + 1 - Kmer::getK());
                        findNPPFast(secondPair, nppvSecond);
                        int firstTipId=0 , secondTipId =0;
                        bool firstFoundOntip =false, secondFoundOnTip =false;
                        for (size_t i = 0; i < nppvFirst.size(); i++) {
                                if (nppvFirst[i].isValid())
                                        if (tipNodes.find(abs(nppvFirst[i].getNodeID()))!=tipNodes.end()){
                                                firstTipId = nppvFirst[i].getNodeID();
                                                firstFoundOntip = true;
                                                break;
                                        }
                        }
                        for (size_t i = 0; i < nppvSecond.size(); i++) {
                                if (nppvSecond[i].isValid())
                                        if (tipNodes.find(abs(nppvSecond[i].getNodeID()))!=tipNodes.end()){
                                                secondTipId = - nppvSecond[i].getNodeID();
                                                secondFoundOnTip = true;
                                                break;
                                        }
                        }
                        if (firstFoundOntip && secondFoundOnTip )
                                if (abs(firstTipId)!= abs(secondTipId)){
                                        //pair <int, int > join =makePairOfTips( firstTipId, secondTipId);
                                        pair <int, int >join = make_pair(firstTipId,secondTipId);
                                        if (pairedEndJoinsTemp.find(join)!=pairedEndJoinsTemp.end())
                                                pairedEndJoinsTemp[join] = pairedEndJoinsTemp[join] +1 ;
                                        else
                                                pairedEndJoinsTemp[join] = 1 ;

                                }

                }
                i ++ ;
        }

        std::map< pair<int, int>, int> pairedEndJoins;
        //sum the frequency with  the reverse direction
        for(auto it = pairedEndJoinsTemp.cbegin(); it != pairedEndJoinsTemp.cend(); ++it)
        {
                if (abs(it->first.first)<abs(it->first.second)){
                        pair <int, int >reverseJoin = make_pair(-it->first.second, -it->first.first);
                        size_t frequency = it->second;
                        if (pairedEndJoinsTemp.find(reverseJoin)!= pairedEndJoinsTemp.end())
                                frequency = frequency  + pairedEndJoinsTemp[reverseJoin];
                        if (frequency > minNumbOfPairs )
                                pairedEndJoins[it->first] = frequency;
                }
        }
        // remove tip pairs which comes with different combinations
        //sort them first based on the frequency
        struct sort_pairs {
                bool operator()(const std::pair< pair<int,int>, int> &left, const std::pair< pair<int,int>, int> &right) {
                        return left.second > right.second;
                }
        };
        i  =0;
        for(auto it = pairedEndJoins.cbegin(); it != pairedEndJoins.cend(); ++it)
        {
                potentialPairs.push_back( make_pair( make_pair(it->first.first, it->first.second), it->second));
                i++;
        }
        std::sort (potentialPairs.begin(),potentialPairs.end(),sort_pairs());
        set<int> nodesHandled;
        pairedEndJoins.clear();

        cout << "( FirstTipID , SecondTipID ) : Frequency"<<endl;
        for(auto it = potentialPairs.cbegin(); it != potentialPairs.cend(); ++it)
        {
                cout << "( " << it->first.first <<" , " <<it->first.second <<" ) :" <<it->second <<endl;
        }

        cout <<"\nNumber of found suggestions: "<<potentialPairs.size() <<endl;
}


pair <int, int > FindGap::makePairOfTips(int firstTipId, int secondTipId)
{

        if (abs(firstTipId) >abs(secondTipId)){
                int temp = firstTipId;
                firstTipId =secondTipId;
                secondTipId =temp;
        }

        pair <int, int >join = make_pair(firstTipId,secondTipId);
        return join;
}

bool FindGap::extendRead(size_t &firstStartIndex, size_t& secondStartIndex, size_t& firstEndInex ,
                         size_t& secondEndIndex, SSNode first,SSNode second ,string &firstRead , string &secondRead)
{

                string firstSeq = first.getSequence();
                string secondSeq = second.getSequence();
                //firstSeq = firstSeq.substr();
                //secondSeq =secondSeq.substr();

                int maxLen1 = firstSeq.length();
                if (maxLen1 >maxSearchSize)
                        maxLen1 = maxSearchSize;
                firstSeq = firstSeq.substr(firstSeq.length()-maxLen1,maxLen1);

                int maxLen2 = secondSeq.length();
                if (maxLen2 >maxSearchSize)
                        maxLen2 = maxSearchSize;
                secondSeq = secondSeq.substr(0,maxLen2);

                string commonSubstr = getLongestCommonSubStr(firstSeq, secondSeq, firstStartIndex, secondStartIndex ) ;
                if (commonSubstr.length()<minExactMatchSize)
                        return false;
                firstStartIndex = firstStartIndex +first.getSequence().length()-maxLen1;

                //cout<<"common  "<<commonSubstr<<endl;
                firstEndInex = firstStartIndex + commonSubstr.length();
                secondEndIndex = secondStartIndex + commonSubstr.length();

                if (firstStartIndex <secondStartIndex){
                        secondStartIndex = secondStartIndex-firstStartIndex;
                        firstStartIndex =0;
                }
                else{
                        firstStartIndex = firstStartIndex -secondStartIndex;
                        secondStartIndex =0;
                }
                if (first.getSequence().length() -firstEndInex < second.getSequence().length()-secondEndIndex){
                        secondEndIndex = secondEndIndex + first.getSequence().length() -firstEndInex;
                        firstEndInex =   first.getSequence().length();

                }else{
                        firstEndInex = firstEndInex +second.getSequence().length()-secondEndIndex ;
                        secondEndIndex = second.getSequence().length();
                }
                firstRead = first.getSequence().substr(firstStartIndex,firstEndInex-firstStartIndex);
                secondRead = second.getSequence().substr(secondStartIndex, secondEndIndex-secondStartIndex);
                if (firstEndInex <first.getSequence().length() || secondStartIndex >0 )
                        return false;
                return true;
}

void FindGap::expandReadByGraphToRight(SSNode first,SSNode second,size_t& firstEndInex , size_t& secondEndIndex, string &firstRead, string &secondRead)
{
        //cout<<"Expand first node to the right for "<<  second.getSequence().length() - secondEndIndex <<" bp "<<endl;
        vector< pair <string ,vector< NodeID>> > bfs;
        vector<NodeID> nodes;
        nodes.push_back(first.getNodeID());
        string p = first.getSequence().substr(first.getSequence().length()- Kmer::getK()+1);
        bfs.push_back(make_pair(p,nodes));
        vector< pair <string ,vector< NodeID>> > result;
        expandNode( second.getSequence().length() -secondEndIndex +Kmer::getK()-1, bfs ,result);
        string secondRight = second.getSequence().substr(secondEndIndex, second.getSequence().length() - secondEndIndex);
        int bestScore = -(second.getSequence().length() - secondEndIndex)*2;
        string bestPath ="";
        for (auto it:result){
                string currentPathstr= it.first;
                string pathRightPart =currentPathstr.substr( Kmer::getK()-1,second.getSequence().length() - secondEndIndex);
                int score =alignment.align(secondRight,pathRightPart);
                if (score >bestScore){
                        bestScore = score;
                        bestPath = pathRightPart;
                }
        }
        firstRead = firstRead + bestPath;
        firstEndInex = firstEndInex + bestPath.length();
        secondEndIndex = secondEndIndex + bestPath.length();
        secondRead = secondRead + secondRight.substr(0,bestPath.length());
}
void FindGap::expandReadByGraphToLeft(SSNode first,SSNode second,size_t& firsStartIndex , size_t& secondStartIndex, string &firstRead, string &secondRead)
{
        //cout<<"Expand first node to the right for "<<  second.getSequence().length() - secondEndIndex <<" bp "<<endl;
        vector< pair <string ,vector< NodeID>> > bfs;
        vector<NodeID> nodes;
        nodes.push_back(first.getNodeID());
        string p = first.getSequence().substr(first.getSequence().length()- Kmer::getK()+1);
        bfs.push_back(make_pair(p,nodes));
        vector< pair <string ,vector< NodeID>> > result;
        expandNode( secondStartIndex  +Kmer::getK()-1, bfs , result);
        string secondleft=second.getSequence().substr(0, secondStartIndex);
        int bestScore = -(secondStartIndex)*2;
        string bestPath ="";
        for (auto it:result){
                string currentPathstr= it.first;
                string pathRightPart =currentPathstr.substr( Kmer::getK()-1,secondStartIndex);
                Nucleotide::revCompl(pathRightPart);
                int score =alignment.align(secondleft,pathRightPart);
                if (score >bestScore){
                        bestScore = score;
                        bestPath = pathRightPart;
                }
        }
        firstRead =bestPath +firstRead;
        firsStartIndex = 0;
        secondRead = secondleft.substr(secondStartIndex-bestPath.length(),bestPath.length()) + secondRead ;
        secondStartIndex = secondStartIndex  - bestPath.length();

}



void FindGap::expandNode( int length, vector< pair <string ,vector< NodeID>> >& bfs,
                          vector< pair <string ,vector< NodeID>> >& result )
{

        int currReadPos = 0;
        vector< pair <string ,vector< NodeID>> > mustVisit;

        for (auto path:bfs){
                SSNode startNode =dbg.getSSNode( path.second.back());
                string currentPathstr= path.first;
                if (currentPathstr.length()<length){
                        if (startNode.getNumRightArcs() !=0 ){
                                for (ArcIt it = startNode.rightBegin(); it != startNode.rightEnd(); it++) {
                                        vector<NodeID> currentPath = path.second;
                                        NodeID nextID = it->getNodeID();
                                        currentPath.push_back(nextID);
                                        const SSNode nextNode = dbg.getSSNode(nextID);
                                        size_t  OLSize = nextNode.getMarginalLength() ;
                                        string nodeOL = nextNode.substr(Kmer::getK()-1, OLSize);
                                        string nextpath = currentPathstr+ nodeOL;
                                        mustVisit.push_back(make_pair(nextpath,currentPath));
                                }
                        }

                }
                else{
                        result.push_back(path);
                }

        }
        if (mustVisit.size() >0){
                bfs.clear();
                bfs = mustVisit;
                expandNode(length,bfs,result);
        }
}
string FindGap::getLongestCommonSubStr(string str1, string str2, size_t & firsStartIndex, size_t & SecondStartIndex)
{

        int m = str1.length();
        int n = str2.length();
        int LCSuff[m+1][n+1];
        int result = 0;
        int maxIndex1 = 0, maxIndex2 = 0;
        for (int i=0; i<=m; i++)
        {
                for (int j=0; j<=n; j++)
                {
                        if (i == 0 || j == 0)
                                LCSuff[i][j] = 0;

                        else if (str1[i-1] == str2[j-1])
                        {

                                LCSuff[i][j] = LCSuff[i-1][j-1] + 1;
                                if (result < LCSuff[i][j]){
                                        result =LCSuff[i][j];
                                        maxIndex1 = i;
                                        maxIndex2 = j;
                                }
                        }
                        else LCSuff[i][j] = 0;
                }
        }
        firsStartIndex = maxIndex1 -result ;
        SecondStartIndex =maxIndex2 -result;
        string commonSubstr = str1.substr(maxIndex1 -result,result);
        return commonSubstr;
}


void FindGap::getPotentialJoins(ComponentHandler& componentHdl, std::map<string,
                                 set<int> > &kmerNodeMap, std::vector<pair<int, int> > &potentialJoin)
{
        std ::set<int> ::iterator setIt1,setIt2;
        std::map<int,int>::iterator it;
        std::map<string, set<int> >::iterator mapIt;

        for(mapIt = kmerNodeMap.begin(); mapIt != kmerNodeMap.end(); mapIt++ ) {
                set <int> nodeSet =mapIt->second;
                int firstCompID = 0 , secondCopID = 0;
                setIt1 =nodeSet.begin();
                while (setIt1 !=nodeSet.end()){
                        setIt2 = setIt1 ;
                        setIt2 ++;
                        while (setIt2 != nodeSet.end()){
                                if ( componentHdl.getComponentID(abs(*setIt1)) !=componentHdl.getComponentID(abs(*setIt2))){
                                        if (std::find(potentialJoin.begin(), potentialJoin.end(),make_pair(*setIt1,*setIt2))==potentialJoin.end() &&
                                                std::find(potentialJoin.begin(), potentialJoin.end(),make_pair(-*setIt1,-*setIt2))==potentialJoin.end() &&
                                                std::find(potentialJoin.begin(), potentialJoin.end(),make_pair(*setIt2,*setIt1))==potentialJoin.end() &&
                                                std::find(potentialJoin.begin(), potentialJoin.end(),make_pair(-*setIt2,-*setIt1))==potentialJoin.end()
                                        ){
                                                int first =*setIt1 , second=*setIt2;
                                                if (dbg.getSSNode(*setIt1).getNumRightArcs()!= 0 )
                                                        first = -first;
                                                if ( dbg.getSSNode(*setIt2).getNumLeftArcs()!= 0)
                                                        second = -second;

                                                potentialJoin.push_back(make_pair(first,second));
                                        }

                                }
                                setIt2 ++;

                        }
                        setIt1 ++;
                }

        }

        /*for (auto it:potentialJoin){
                std::pair<int, int> v_temp = it;
                int firstNodeID=v_temp.first;
                int secondNodeID = v_temp.second;
                cout <<"match : " << firstNodeID<<"," << secondNodeID << endl;
                SSNode first=dbg.getSSNode(firstNodeID);
                SSNode second=dbg.getSSNode(secondNodeID);
                cout <<first.getSequence() <<endl;
                cout<<second.getSequence()<<endl;
        }*/


}

void FindGap::loadKmerMap(set<int> &tipNodes,std::map<string, set<int> > & kmerNodeMap,size_t overlapSize){
        Kmer::setWordSize(overlapSize);
        cout << "Loading the kmers in the component into the table ..." <<endl;
        std::map<string, set<int> >::iterator mapIt,mapTtR;
        for (NodeID nodeid :tipNodes){
                SSNode node = dbg.getSSNode(nodeid);
                string DNA_sequence = node.getSequence();
                vector<string> sequences;
                int maxLen = DNA_sequence.length();
                if (maxLen >maxSearchSize)
                        maxLen = maxSearchSize;
                if (node.getNumLeftArcs() ==0 && node.getNumRightArcs() >0)
                        sequences.push_back(DNA_sequence.substr(0,maxLen));

                if (node.getNumRightArcs() ==0 && node.getNumLeftArcs() >0)
                        sequences.push_back(DNA_sequence.substr(DNA_sequence.length()-maxLen,maxLen));

                if (node.getNumLeftArcs() ==0 && node.getNumRightArcs() ==0 ){
                        if (maxSearchSize < DNA_sequence.length()*2)
                                sequences.push_back(DNA_sequence);
                        else{
                                sequences.push_back(DNA_sequence.substr(0,maxLen));
                                sequences.push_back(DNA_sequence.substr(DNA_sequence.length()-maxLen,maxLen));
                        }
                }
                for (auto seq:sequences){

                        for (KmerIt it(seq);it.isValid(); it++)
                        {
                                Kmer kmer = it.getKmer();
                                //cout << kmer <<endl;
                                mapIt = kmerNodeMap.find(kmer.str());
                                mapTtR= kmerNodeMap.find(kmer.getReverseComplement().str());
                                if (mapIt != kmerNodeMap.end() ){
                                        mapIt->second.insert(nodeid);

                                }
                                else{
                                        if (mapTtR !=kmerNodeMap.end()){
                                                mapTtR->second.insert(-nodeid);

                                        }else {
                                                set <int> nodeSet;
                                                nodeSet.insert(nodeid);
                                                kmerNodeMap.insert(make_pair(kmer.str(), nodeSet));
                                                set <int> nodeSetR;
                                                nodeSetR.insert(-nodeid);
                                                kmerNodeMap.insert(make_pair(kmer.getReverseComplement().str(), nodeSetR));

                                        }
                                }
                             }
                }
        }

        for(mapIt = kmerNodeMap.begin(); mapIt != kmerNodeMap.end(); mapIt++ ) {
                set <int> nodeSet =mapIt->second;
                if (nodeSet.size()<2){
                        kmerNodeMap.erase(mapIt);

                }/*else{
                std::cout << mapIt->first <<":\t" <<endl;
                for (set<int>::iterator i = nodeSet.begin(); i != nodeSet.end(); i++) {
                        cout << *i <<", ";
        }
        cout << std::endl ;
        }*/
        }

        Kmer::setWordSize(settings.getK());
}
void FindGap::assignTipsToComponent(set<int> &tipNodes,ComponentHandler &componentHdl,std::map<int, int > & nodeComponentMap)
{
        for (int i= 0; i< componentHdl.numberOfcomponents; i++){
                std::set<NodeID> nodes = componentHdl.components[i].getNodeIDSet();
                for (auto it:nodes){
                        if (tipNodes.find(it) != tipNodes.end())
                                nodeComponentMap[it] = i;
                        if (tipNodes.find(-it) != tipNodes.end())
                                nodeComponentMap[-it] = i;
                }
        }
        map<int, int>::iterator it;
        /*for(it = nodeComponentMap.begin(); it != nodeComponentMap.end(); it++ ) {

                std::cout << it->first  // string (key)
                << ':'
                << it->second   // string's value
                << std::endl ;
        }*/
}

void FindGap::findTips(set<int> &tipNodes, set <int>& eligibleNodeSet)
{
         for (NodeID id =1; id <= dbg.getNumNodes(); id++) {
                SSNode node = dbg.getSSNode(id);
                if (!node.isValid())
                        continue;

                // check for dead ends
                bool leftDE = (node.getNumLeftArcs() == 0);
                bool rightDE = (node.getNumRightArcs() == 0);
                if (!leftDE && !rightDE)
                        continue;
                if (eligibleNodeSet.find(node.getNodeID())!=eligibleNodeSet.end() || eligibleNodeSet.find(-node.getNodeID())!=eligibleNodeSet.end())
                        tipNodes.insert(node.getNodeID());
         }
}

void FindGap::findPotentialTips(set<int> &tipNodes)
{
         for (NodeID id =1; id <= dbg.getNumNodes(); id++) {
                SSNode node = dbg.getSSNode(id);
                if (!node.isValid())
                        continue;

                // check for dead ends
                bool leftDE = (node.getNumLeftArcs() == 0);
                bool rightDE = (node.getNumRightArcs() == 0);
                if (!leftDE && !rightDE)
                        continue;
                tipNodes.insert(node.getNodeID());
         }
}

void FindGap:: findComponentsInGraph(ComponentHandler& componentHdl)
{
        cout << "Finding disjoint components in the graph ..." <<endl;
        int srcID=1;
        int i=0;
        set<NodeID> nodesHandled;
        for (i =srcID ; i <= dbg.getNumNodes(); i++){
                set<NodeID> currentSetNodes;
                size_t componentSize =0;
                if (!dbg.getSSNode(i).isValid())
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
                        SSNode thisNode = dbg.getSSNode(thisID);
                        componentSize = componentSize + thisNode.getMarginalLength();
                        for (ArcIt it = thisNode.rightBegin(); it != thisNode.rightEnd(); it++) {
                                SSNode rNode = dbg.getSSNode(it->getNodeID());
                                if (!rNode.isValid())
                                        continue;
                                if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                        continue;
                                nodeDepth.insert(pair<size_t, NodeID>(thisDepth + thisNode.getMarginalLength(), it->getNodeID()));
                        }
                        for (ArcIt it = thisNode.leftBegin(); it != thisNode.leftEnd(); it++) {
                                SSNode lNode =dbg.getSSNode(it->getNodeID());
                                if (!lNode.isValid())
                                        continue;
                                if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                        continue;
                                nodeDepth.insert(pair<size_t, NodeID>(thisDepth + thisNode.getMarginalLength(), it->getNodeID()));
                        }
                }
                if (currentSetNodes.size()>0 && componentSize > minComponentSize && componentSize <maxComponentSize)
                {
                        componentHdl.addComponent(currentSetNodes);
                }
                currentSetNodes.clear();
        }
        cout << "Number of disjoint components in the graph with size higher than "<< minComponentSize <<" bp " <<componentHdl.components.size() <<endl;
}


void FindGap::extractBreakpointSubgraph( std::string breakpointFileName,string  tempDir){
        cout << "Creating kmer lookup table... ";
        dbg.buildKmerNPPTable();

        std::ifstream input(breakpointFileName);
        vector< pair<string, string> > breakpoints;
        if (!input.good()) {
                std::cerr << "Error opening: " << breakpointFileName << std::endl;
                return;
        }
        std::string line, id ="", DNA_sequence ="" ;
        while (std::getline(input, line).good()) {
                if (line[0] == '>') {
                        id = line.substr(1);
                        if (DNA_sequence!=""){
                                breakpoints.push_back( make_pair(id, DNA_sequence));
                        }
                        DNA_sequence.clear();
                }
                else if (line[0] != '>'){
                        DNA_sequence += line;
                }
        }
        if (DNA_sequence!=""){
                breakpoints.push_back( make_pair(id, DNA_sequence));
        }

        for (unsigned int i = 0 ;i <breakpoints.size(); i++){
                pair<string , string> breakPoint= breakpoints [i];
                set<int> nodeSet;
                //extractNodeSetbySequence( breakPoint.second,nodeSet);
                RefComp refComp(breakpointFileName);
                vector<NodeChain> trueNodeChain;
                refComp.getTrueNodeChain(dbg, trueNodeChain);
                writeCytoscapeGraph( tempDir+breakPoint.first,trueNodeChain,1);

        }
}
void FindGap::writeCytoscapeGraph(const std::string& filename,
                                  vector<NodeChain>& nodeChain, size_t maxDepth) const
{
        // a map containing nodeIDs to handle + their depth
        priority_queue<PathDFS, vector<PathDFS>, PathDFSComp> nodeDepth;
        // set of nodes that were handled
        set<NodeID> nodesHandled;

        // if default input values are provided, write the entire graph

        for (unsigned int i =0; i<nodeChain.size();i++){
                std::vector<int>::iterator it;
                NodeChain chain = nodeChain[i];
                for (it =chain.begin() ;it!=chain.end();it++)
                {        // else check if seedNode is valid
                        NodeID seedNodeID = *it;
                        if ((abs(seedNodeID) > dbg.getNumNodes()) || !dbg.getSSNode(seedNodeID).isValid()) {
                                cerr << "WARNING: trying to use an invalid node as a "
                                "seed in writeCytoscapeGraph!" << endl;
                                return;
                        }
                        nodeDepth.push(PathDFS(seedNodeID, 0));
                }
        }

        // map of all nodes in the local graph and its depth
        ofstream ofs((filename + ".arcs").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".arcs" << " for writing" << endl;
        ofs << "Source node\tTarget node\tArc coverage" << endl;

        // A) write all arcs
        while (!nodeDepth.empty()) {
                // get and erase the current node
                PathDFS currTop = nodeDepth.top();
                nodeDepth.pop();
                NodeID thisID = currTop.nodeID;
                size_t thisDepth = currTop.length;

                // if the node was already handled, skip
                if (nodesHandled.find(thisID) != nodesHandled.end())
                        continue;

                // if we're too far in the graph, stop
                if (thisDepth > maxDepth) {
                        nodesHandled.insert(thisID);
                        continue;
                }

                // write the right arcs
                SSNode thisNode = dbg.getSSNode(thisID);
                for (ArcIt it = thisNode.rightBegin(); it != thisNode.rightEnd(); it++) {
                        SSNode rNode = dbg.getSSNode(it->getNodeID());
                        if (!rNode.isValid())
                                continue;
                        if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                continue;
                        ofs << thisID << "\t" << it->getNodeID() << "\t" << it->getCoverage() << "\n";
                        PathDFS nextTop(it->getNodeID(), thisDepth + 1);
                        nodeDepth.push(nextTop);
                }

                // write the left arcs
                for (ArcIt it = thisNode.leftBegin(); it != thisNode.leftEnd(); it++) {
                        SSNode lNode = dbg.getSSNode(it->getNodeID());
                        if (!lNode.isValid())
                                continue;
                        if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                continue;
                        if (it->getNodeID() != thisID)  // avoid writing this arc twice
                                ofs << it->getNodeID() << "\t" << thisID << "\t" << it->getCoverage() << "\n";
                        PathDFS nextTop(it->getNodeID(), thisDepth + 1);
                        nodeDepth.push(nextTop);
                }

                // mark this node as handled
                nodesHandled.insert(thisID);

        }
        ofs.close();

        // B) write all nodes
        ofs.open((filename + ".nodes").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".nodes" << " for writing" << endl;
        ofs << "Node ID\tMarginal length\tNum left arcs\tNum right arcs\tTrue multiplicity\tEstimated multiplicity"
               "\tKmer coverage\tRead start coverage\tSequence\tinPath\tpos" << "\n";
        for (set<NodeID>::iterator it = nodesHandled.begin(); it != nodesHandled.end(); it++) {
                SSNode node = dbg.getSSNode(*it);
                int thisTrueMult = 0;
                bool is_in = false;
                int index = -1 ;
                for (unsigned int i =0; i<nodeChain.size();i++){
                        NodeChain chain = nodeChain[i];
                        std::vector<int>::iterator vit;
                        vit = find (chain.begin(), chain.end(), *it);
                        if (vit != chain.end()) {
                                is_in = true;
                                index = std::distance(chain.begin(), vit);
                                break ;
                        }
                }
                ofs << *it << "\t" << node.getMarginalLength() << "\t"
                    << (int)node.getNumLeftArcs() << "\t" << (int)node.getNumRightArcs() << "\t"
                    << thisTrueMult << "\t" << "0" << "\t"
                    << node.getAvgKmerCov()
                    << "\t" << double(node.getReadStartCov()/node.getMarginalLength())
                    << "\t" << node.getSequence()
                    <<"\t"<<is_in <<"\t"<<index<<"\n";
        }
        ofs.close();
}

int main(int argc, char** args)
{
        enum command {help,cytoscape,breakpoint,closeGap};
        if (argc <7){
                cout << "Correct usage:"<<endl;
                cout << "./findGap -n nodeFileName -a arcFileName -m metaDataFileName " << endl;
                return EXIT_SUCCESS;
        }
        string nodeFileName = "";
        string arcFileName = "";
        string metaDataFileName ="";
        unsigned int kmerSize =21 ;
        string alignmentFile ="";
        string tmpDir = "." ;
        string breakpointFileName ="";
        command c;
        cout <<args[1]<<endl;
        string arg(args[1]);
        if ((arg == "-h") || (arg == "--help")) {
                c = help;
        }else if (arg == "cytoscape"){
                c = cytoscape;
        }else if ( arg =="breakpoint"){
                c = breakpoint;
        }
        else if (arg =="closeGap"){
                c = closeGap;
        }
        if (c==help){
                cout << "Correct usage:"<<endl;
                cout << "./findGap cytoscape -n nodeFileName -a arcFileName -m metaDataFileName " << endl;
                cout << "./findGap breakpoint -n nodeFileName -a arcFileName -m metaDataFileName " << endl;
                return EXIT_SUCCESS;
        }

        for (int i = 2; i < argc; i++) {
                string arg(args[i]);

                if ((arg == "-n") || (arg == "--nodeFileName")) {
                        i = i+1;
                        nodeFileName =  args[i];
                } else if ((arg == "-a") || (arg == "--arcFileName")) {
                        i = i+1;
                        arcFileName = args[i];
                } else if ((arg == "-m") || (arg == "--metaDataFileName")) {
                        i = i+1;
                        metaDataFileName = args[i];
                }
                else if ((arg == "-k") || (arg == "--kmerSize")) {
                        i = i+1;
                        kmerSize = atoi(args[i]);
                }
                else if( (arg == "-c") || (arg == "--correctedFile")){
                        i = i+1;
                        alignmentFile = args[i];
                }
                 else if ((arg == "-p") || (arg == "--tempDir")) {
                        i = i+1;
                        tmpDir = args[i];
                        if (!tmpDir.empty()) {
                                if ((tmpDir.back() != '/') && (tmpDir.back() != '\\'))
                                        tmpDir.push_back('/');
                        }
                }
                else if ((arg == "-f") || (arg == "--breakpointFileName")) {
                        i = i+1;
                        breakpointFileName = args[i];
                }
        }

        if (!nodeFileName.empty() && !arcFileName.empty() && !metaDataFileName.empty()){
                FindGap findGap(nodeFileName,arcFileName,metaDataFileName, alignmentFile,kmerSize,tmpDir);
                switch (c)
                {
                        case cytoscape :
                                findGap.dbg.writeCytoscapeGraph("cytoscape");
                                break;
                        case breakpoint :
                                findGap.extractBreakpointSubgraph(breakpointFileName,tmpDir);
                                break;
                        case closeGap :
                                findGap.closeGaps(nodeFileName, arcFileName,metaDataFileName);
                                break;
                }


        }
        cout <<"Finished successfully"<<endl;
        return EXIT_SUCCESS;
}

