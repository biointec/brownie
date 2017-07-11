

#include <string>
#include "graph.h"
#include <queue>
#include "settings.h"
#include "refcomp.h"
#include "component.h"
class FindGap {
private:
        AlignmentJan alignment;
public :
        DBGraph dbg;
        Settings settings;
        string correctedFile;
        size_t overlapSize;
        size_t maxSearchSize;
        size_t kmerSize;


        /**
         * Default constructor
         */

        FindGap (string nodeFileName, string arcFileName, string metaDataFileName,string alignmentFile,unsigned int kmerSize =21 ,string tempDir=".");


        /**
         * Write a cytoscape graph of the current graph over a set of nodes
         * @param filename Filename of the cytoscape graph
         * @param nodeSet set of seeds node identifier
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void writeCytoscapeGraph(const std::string& filename,vector<NodeChain>& nodeChain, size_t maxDepth=0) const;


        /**
         * extract the the subgraph nearbye given sequence in a fasta file
         * @param dbg A const-ref to the de Bruijn graph
         * @param breakpointFileName the fasta file contains sequences which are supposed to be in breakpoint area
         * @param tempDir the temporary directory to save Cytoscape files
         */
        void extractBreakpointSubgraph(std::string breakpointFileName, std::string tempDir);
        /**
         * Calculate the true node multiplicity
         * @param dbg A const-ref to the de Bruijn graph
         * @param multiplicity Multiplicity vector
         */
        void getNodeMultiplicity(std::vector<size_t>& multiplicity);

        /**
         * find the gap in the graph, recommend connection between tips
         *
         */
        void closeGaps();
        /**
         * extract kmers in the graph, which are in tip nodes and save them in a dictionary
         * @param tipNodes A set of tips nodes
         * @param kmerNodeMap a map which shows each kmer exists in which nodes
         */
        void loadKmerMap(set<int>& tipNodes,std::map<string, set<int> >&  kmerNodeMap, size_t overlapSize);

        /**
         * find the tips in the graph
         * @param tipNodes A set of tips nodes
         * @param eligibleNodeSet the tips should be among these set
         */
        void findTips(set<int> &tipNodes, set <int>& eligibleNodeSet);
        /**
         * find disjoined component in the graph
         * @param componentHdl keeps the specification of the disjoined component in the graph
         */

        void findComponentsInGraph(ComponentHandler& componentHdl);
        /**
         * assign each tip to a component ID which that tip exists there and keep it in a map structure
         * @param tipNodes A set of tips nodes
         * @param componentHdl keeps the specification of the disjoined component in the graph
         * @param nodeComponentMap stors tips are in which component
         */
        void assignTipsToComponent(set<int> &tipNodes,ComponentHandler& componentHdl,std::map<int, int >&  nodeComponentMap);
        /**
         * retun all pairs of tips which in different component which has a common kmer
         * @param potentialJoin A vector contains all the  pairs of tips which has a common kmer
         * @param kmerNodeMap A map which stores all the kmers beside to the nodes which those kmer exist
         * @param componentHdl keeps the specification of the disjoined component in the graph, here we want to check if two tips are in the different component or not
         */
        void getPotentialJoins( ComponentHandler& componentHdl,
                                         std::map<string,set<int> >& kmerNodeMap, std::vector<pair<int, int> >& potentialJoin);
         /**
         * extract the path for each tip, align the sequence in the path and eliminate those tips which cant be join
         * @param potentialJoin A vector contains all the  pairs of tips which has a common kmer
         */
        bool eliminateFakeJoins(std::vector<pair<int, int> > &potentialJoin, std::map< pair<int, int>, int> &pairedEndJoins,  ComponentHandler& componentHdl);
          /**
         * retun the Longest Common Substring
         * @param str1 first input string
         * @param str2 second input string
         * @param firsStartIndex start index of common substring in the first string
         * @param SecondStartIndex start index of common substring in the second string
         */

        string getLongestCommonSubStr(string str1, string str2, size_t& firsStartIndex, size_t & SecondStartIndex);
         /**
         * expands the node in the graph and parse the graph in breadth first search in a given distance from root
         * @param length distance from root that the graph should be parsed
         * @param bfs root node
         * @param result all different path start from root in that distance
         */
        void expandNode( int length, vector< pair <string ,vector< NodeID>> >& bfs , vector< pair <string ,vector< NodeID>> > &result );
         /**
         * expand the common string between two nodes to the right if the first node can be expanded to the right.
         * @param first first node
         * @param second second node
         * @param firsStartIndex start index of common substring in the first read
         * @param secondStartIndex start index of common substring in the second read
         * @param firstRead the first string after expansion
         * @param secondRead the second string after expansion
         */
        void expandReadByGraphToRight(SSNode first,SSNode second, size_t& firsStartIndex,size_t& secondStartIndex, string &firstRead, string &secondRead);
         /**
         * expand the longest common string between two nodes to the left if the first node can be expanded to the left.
         * @param first first node
         * @param second second node
         * @param firsStartIndex start index of common substring in the first read
         * @param secondStartIndex start index of common substring in the second read
         * @param firstRead the first string after expansion
         * @param secondRead the second string after expansion
         */
        void expandReadByGraphToLeft (SSNode first,SSNode second,size_t& firsStartIndex , size_t& secondStartIndex, string &firstRead, string &secondRead);
         /**
         * expand the common substring of two given strings (firstRead nad secondRead) to the left and right by looking at two nodes, without expanding the graph.
         *
         * @param firsStartIndex start index of common substring in the first read
         * @param secondStartIndex start index of common substring in the second read
         * @param firstEndInex end index of common substring in the first read
         * @param secondEndIndex end index of common substring in the second read
         * @param first first node
         * @param second second node
         * @param firstRead the first string after expansion
         * @param secondRead the second string after expansion
         */
        void extendRead(size_t &firsStartIndex, size_t& secondStartIndex, size_t& firstEndInex , size_t& secondEndIndex, SSNode first,SSNode second ,string &firstRead , string &secondRead);
         /**
         * if a paired end read align to two different component, it keeps it as an indiction of a connection between two component
         * @param readFileName read corrected file in fastq format which instead of corrected reads show the path which those reads align
         * @param componentHdl keeps the specification of the disjoined component in the graph
         */
         void extractPairedComponents(string readFileName , ComponentHandler& componentHdl , std::map< pair<int, int>, int> &pairedEndJoins);
        /**
         * it splits the line of nodes to the vector of nodeID
         * @param line the input string which has the list of nodes

         */
        void splitLineToNodes(string line, vector<NodeID> &nodeIDs);


        //int GetBesAlternativePath(  NodeID root, string rightPart ,string& bestAlternative, vector<NodeID> & bestPathNodes);

};
