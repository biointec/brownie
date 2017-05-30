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

#ifndef DBGRAPH_H
#define DBGRAPH_H
#include "global.h"
#include "ssnode.h"
#include "dsnode.h"
#include <deque>
#include "essaMEM-master/sparseSA.hpp"


// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

enum MapType { SHORT_MAP, LONG_MAP };
enum SearchResult { HAVE_SOLUTION, SEARCH_EXHAUSTED, NO_SOLUTION };

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class Settings;
class ReadLibrary;
class DSNode;
class Arc;
class Settings;
class KmerNodeTable;
class NodePosPair;
class NodeEndTable;
class NodeEndRef;
class LibraryContainer;

// ============================================================================
// GRAPH CLASS
// ============================================================================

class DBGraph {

private:
    // table for the kmers in the genome
    typedef google::sparse_hash_set<string> HashSet;
    typedef HashSet::const_iterator HashSetIt;

    KmerNodeTable* table;                   // kmer node table

    /**
     * Parse a buffer of reads and store kmers in temporary buffers per thread
     * @param readBuffer Input read buffer
     */
    void parseReads(size_t thisThread,
                    std::vector<std::string>& readBuffer);

    /**
     * Entry routine for worker thread
     * @param myID Unique threadID
     */
    void workerThread(size_t myID, LibraryContainer* inputs);

    // ====================================================================
    // COVERAGE.CPP PRIVATE
    // ====================================================================

    /**
     * Count the number of kmer occurences per node for a certain readfile
     * @param filename File name of the input read file
     * @param table Table containing the kmer -> nodeID mapping
     */
    void countReadFrequency(const ReadLibrary& input,
                            const KmerNodeTable &table);

    /**
     * Get an initial estimate for the node coverage, based on the
     * 15% largest nodes
     * @return An estimate for the coverage
     * @param kmerFreq The number of kmers per node (output)
     */
    double getInitialEstimateForCoverage(const ReadLibrary& input,
                                         vector<size_t> &kmerFreq) const;

    /**
     * Based on preset node multiplicities and node kmer occurences,
     * estimate the expected kmer occurence (mu)
     * @return The expected average occurence
     */
    double estimateReadStartCoverage(const ReadLibrary& input,
                                     const vector<size_t> & kmerFreq) const;

    // ====================================================================
    // OTHER STUFF
    // ====================================================================


    /**
     * Get the unique node extending a given node to the left
     * @param node Node to be extended (input)
     * @param leftNode Node that left-overlaps with the given node (output)
     * @return True if a unique left-overlapping node is found
     */
    bool getLeftUniqueSSNode(const SSNode &node, SSNode &leftNode) const;

    /**
     * Get the unique node extending a given node to the right
     * @param node Node to be extended (input)
     * @param rightNode Node that right-overlaps with the given node (output)
     * @return True if a unique right-overlapping node is found
     */
    bool getRightUniqueSSNode(const SSNode &node, SSNode &rightNode) const;

    /**
     * Increment the coverage of the arc that span two given kmers
     * @param left Left kmer reference
     * @param right Right kmer reference
     */
    void increaseCoverage(const NodeEndRef &left, const NodeEndRef &right);

    /**
     * Convert a vector of overlapping nodes to a string
     * @param nodeSeq A deque of overlapping nodes
     * @param output An stl string (output)
     */
    void convertNodesToString(const std::vector<NodeID> &nodeSeq,
                              std::string &output);

    /**
     * Convert a deque of overlapping nodes to a string
     * @param nodeSeq A deque of overlapping nodes
     * @param output An stl string (output)
     */
    void convertNodesToString(const std::deque<SSNode> &nodeSeq,
                              std::string &output);

    /**
     *
     *
     */
    void convertNodesToString(const std::vector<NodeID> &nodeSeq,
                              int startPos, int stopPos,
                              string &output);

    /**
     * Thread reads through nodes to obtain arc coverage
     * @param filename File name of the input file
     * @param numReads Number of reads in the input file
     * @param table Node ends table
     */
    void threadThroughReads(const std::string& filename, size_t numReads,
                            const NodeEndTable &table);

    void markPairedArcs(const std::vector<NodeID>& seq);

    // ====================================================================
    // VARIABLES
    // ====================================================================

    const Settings &settings;     // settings object

    DSNode *nodes;          // graph nodes
    Arc *arcs;              // graph arcs

    NodeID numNodes;        // number of nodes
    NodeID numArcs;         // number of arcs

    MapType mapType;

    double coverage;//=100;

    double totalRAM;

#ifdef DEBUG
        std::vector<string> reference;
        std::vector<void* > refST;
        std::vector<int> trueMult;
#endif

public:

    static const DBGraph* graph;
    double estimatedKmerCoverage;
    double estimatedMKmerCoverageSTD;
    double minCertainVlueCov;
    double minSafeValueCov;
    double minRedLineValueCov;
    double certainVlueCov;
    double safeValueCov;
    double redLineValueCov;
    double cutOffvalue;

    double readLength;
    //this variable shows the Maximum node size which can be deleted, which is calculated based on the read readLength
    size_t maxNodeSizeToDel;
    int updateCutOffValueRound;

    size_t n50;
    size_t sizeOfGraph;


    /**
     * Careful concatenation, taking into account the estimated multiplicity
     */
    void initialize();
/**
 *fucntion associated to bubble bubbleDetection
 *
 *
 */
    bool bubbleDetection(int round);
    vector<pair<SSNode, SSNode> >  ExtractBubbles(SSNode rootNode,std::set<NodeID>& visitedNodes , std::set<Arc *>&visitedArc);
    bool removeBubble(SSNode &prevFirstNode ,SSNode& extendFirstNode,size_t &TP,size_t &TN,size_t &FP,size_t &FN,size_t & numOfDel);
    void extractPath(NodeID currID, const vector<NodeID>& prevNode) const;
    vector<NodeID> getPath(NodeID currID, const vector<NodeID>& prevNode) const;
    bool removeNotSingleBubbles(  SSNode &prevFirstNode ,SSNode& extendFirstNode, size_t &TP,size_t &TN,size_t &FP,size_t &FN,size_t & numOfDel);
    bool whichOneIsbubble(SSNode rootNode,bool &first, SSNode &prevFirstNode ,SSNode& extendFirstNode, bool onlySingle, double threshold);
    bool whichOneIsbubble(SSNode rootNode,bool &first, SSNode &prevFirstNode ,SSNode& extendFirstNode, bool onlySingle);
    /**
 * this routine is used to check if the currNode can be a parallel path to some other nodes
 * the root node is also given.
 * this routine is used as a subroutin of deleteUnreliableNodes
 * @return true if the given nodes is a bubble initiate from rootNode
 *
 */
    bool nodeIsBubble(SSNode node, SSNode currNode);
    vector<pair<vector<NodeID>, vector<NodeID>> >  searchForParallelNodes(SSNode node,vector<NodeID> &visited, vector<NodeID> &prevNode,vector<NodeID> &nodeColor, int depth);
    vector<pair<vector<NodeID>, vector<NodeID>> > searchForParallelNodes(SSNode node, int depth);
    bool hasLowCovNode(SSNode root);


 /**
 *fucntion associated to graph graph Purification
 *
 *
 */
    bool removeNode(SSNode &rootNode);
    /**
     * remove links between two nodes
     * @param sourceNode the first node
     * @param rootNode the second node
     */
    bool removeLinks(SSNode & sourceNode, SSNode & rootNode);
    bool mergeSingleNodes();
    /**
     * this routine do the following things
     * 1. calculation of avg and std of node kmer coverage
     * 2. calculation of read start coverage mean
     * 3. for every node it calculates the MULTIPLICITY
     * 4. for every node it calculates the certainity of our guess about MULTIPLICITY
     * 5. for every node it calculates the inCorrctnessRatio of our guess about MULTIPLICITY
     *
     */

    void extractStatistic(int round);
    /**
     * check nodes to see wether is reliable or not
     * a node is reliable if it has a high confidenceRatio and low inCorrctnessRatio
     * these terms are defined in extractStatistic routine
     * @return true if the node is reliable
     */


    bool checkNodeIsReliable(SSNode node);
    /**
     *
     * this routine loops over all nodes, for a reliable node N with MULTIPLICITY(M)
     * it should have at most M outgoing arcs. Therefore we find extra nodes with low coverage
     * or extra nodes wich are appeared as tips or bubbles.
     * @return true if any changes happens
     * the second part of this routine looks for the adjacent nodes with same MULTIPLICITY
     * if they have some outgoing or ingoing arcs with low coverage those nodes should be
     * deleted and these two adjacent reliable nodes should be connected later.
     */
    bool deleteUnreliableNodes();
    bool deleteExtraAttachedNodes();
    bool connectSameMulNodes();
    bool deleteSuspiciousNodes();
    bool clipTips(int round);




    // void filterCoverage();

    /**
     * Default constructor
     */
    DBGraph(const Settings& settings);

    /**
     * Destructor
     */
    ~DBGraph();

    /**
     * Get a const-reference to settings
     * @return A const-reference to settings
     */
    const Settings& getSettings() const {
        return settings;
    }

    // ====================================================================
    // SOLUTIONCOMP.CPP PUBLIC
    // ====================================================================

    void writeCytoscapeGraph(int ID);

    void writeLocalCytoscapeGraph(int ID, NodeID nodeID,size_t maxDepth);

    void readReferenceGenome();

    /**
     * Compare the current graph to the solution
     * @param filename Filename of file containing true multiplicities
     */
    void compareToSolution(const string& filename,bool load);
    /**
     * checking the continuity of kmers in the graph by
     * looping over the kmers in the genome. and reporting breakpoints.
     *
     * */
     void findBreakpoint();
     /**
      * Check if two nodes are connected in the graph
      * @param curNode the current node
      * @param preNode the previous node
      * @return true if two nodes are connected, otherwise false
      *
      */
     bool checkConnectivity(SSNode curNode, SSNode preNode);
     /**
      * looks inside genome set table for a given kmer
      * @param kmer the keyword to search in the table
      * @return true if the word exists in the table otherwise false
      */
     bool lookup( Kmer kmer, HashSet genomeKmerSet);
     /**
      * fill an object of  google sparse hash table with kmers whcih are exist in = the genome
      * @return void
      */
     void fillGenomeSet();

     size_t findAllTrueOccurences(const string& str,
                                 std::vector<std::vector<size_t> >& pos,
                                 std::vector<std::vector<size_t> >& posRC) const;

    /**
     * Find the real path from srcID to dstID
     * @param srcID Source node identifier
     * @param dstID Destination node identifier
     * @param distance Maximum distance between srcID and dstID
     * @param result Vector in which to store the path
     */
    void getActualPath(NodeID srcID, NodeID dstID, size_t distance,
                       std::vector<NodeID>& result) const;

    bool getDistance(NodeID srcID, NodeID dstID, double maxDistance,
                     double &result) const;

    void getTrueOccurencesForReadPair(std::string &start,
                                      std::string &stop,
                                      std::vector<string>& output,
                                      size_t maxSize);

    // ====================================================================
    // COVERAGE.CPP PUBLIC
    // ====================================================================

    /**
     * Count the number of kmer occurences for each node
     */
    void countReadFrequency(const ReadLibrary& input);

    /**
     * Set the node multiplicity based on an estimate for the coverage
     * based on expectation-maximization (EM)
     */
    void countNodeandArcFrequency(LibraryContainer &inputs);

    /**
     * Filter the graph, based on coverage
     */
    bool filterCoverage( float round);
    /*This routine calculate the value of CutOff-cov in the first call,
     * based on cutOff value three other values will be updated in each round of calling. they increase in each round
     * cutOffvalue*.3< certainValue< cutOffvalue*.8 : the most conservative value , used to remove nodes only based on coverage
     * cutOffvalue*.7< safeValue <cutOffvalue  always lower than cutOff , used for remoing joined tips or isolated nodes.
     * cutOffvalue <redLineValue< cutOffvalue*1.5, increasing incrementally in while loop
     *
     *
     *
     */
    void updateCutOffValue(int round);
    void plotCovDiagram(vector<pair< pair< int , int> , pair<double,int> > >& frequencyArray);
    void makeSampleReadFile(float num);
    size_t getLowestArcMultiplicity(NodeID left, NodeID right);

    // ====================================================================
    // COVERAGE.CPP PUBLIC
    // ====================================================================

    /**
     * Map reads to the graph
     */
    void mapReads();

    // ====================================================================
    // GRAPH CONSTRUCTION
    // ====================================================================

    /**
     * Get the number of nodes
     * @return The number of nodes
     */
    NodeID getNumNodes() const {
        return numNodes;
    }

    /**
     * Get the number of arcs
     * @return The number of arcs
     */
    NodeID getNumArcs() const {
        return numArcs;
    }

    /**
     * Create a graph from file
     * @param nodeFilename Filename for the nodes
     * @param arcFilename Filename for the arcs
     * @param metaDataFilename Filename for the metadata
     */
    void createFromFile(const std::string& nodeFilename,
                        const std::string& arcFilename,
                        const std::string& metaDataFilename);

    /**
     * Thread through the reads
     */
    void threadThroughReads();

    /**
     * Clear all nodes and arcs in this graph
     */
    void clear() {
        delete [] nodes;
        delete [] arcs;
        nodes = NULL;
        arcs = NULL;
        numNodes = numArcs = 0;
    }

    /**
     * Perform graph simplification
     */
    bool simplyfyGraph();

    /**
     * Get a reference to a double stranded node, given the nodeID
     * @param nodeID Identifier for the node
     * @return Reference to the node
     */
    DSNode& getDSNode(NodeID nodeID) const {
        assert(nodeID > 0 && nodeID <= numNodes);
        return nodes[nodeID];
    }

    /**
     * Get a single stranded node, given the nodeID
     * @param nodeID Identifier for the node
     * @return A Single Stranded node
     */
    const SSNode getSSNode(NodeID nodeID) const {
        NodeID uNodeID = abs(nodeID);
        assert(uNodeID != 0 && uNodeID <= numNodes);
        return SSNode(nodes + uNodeID, nodeID);
    }

    /**
     * Get a single stranded node, given the nodeID
     * @param nodeID Identifier for the node
     * @return A Single Stranded node
     */
    SSNode getSSNode(NodeID nodeID) {
        NodeID uNodeID = abs(nodeID);
        assert(uNodeID != 0 && uNodeID <= numNodes);
        return SSNode(nodes + uNodeID, nodeID);
    }

    /**
     * Simply write graph to file
     */
    void writeGraph(const std::string& nodeFilename,
                    const std::string& arcFilename,
                    const std::string& metaDataFilename);

    /**
     * Write graph to file (binary version)
     * @param nodeFilename node filename
     * @param arcFilename arc filename
     * @param metaDataFilename Metadata filename
     */
    void writeGraphBin(const std::string& nodeFilename,
                       const std::string& arcFilename,
                       const std::string& metaDataFilename);

    /**
     * Load graph to file (binary version)
     * @param nodeFilename node filename
     * @param arcFilename arc filename
     * @param metaDataFilename Metadata filename
     */
    void loadGraphBin(const std::string& nodeFilename,
                      const std::string& arcFilename,
                      const std::string& metaDataFilename);

    size_t updateGraphSize();

    /**
     * Extract nodes from graph and write to file
     */



    /**
     * Check a graph for consistency
     */
    void sanityCheck();
    bool removeIncorrectLink();

    //comment by mahdi
    //map<NodeID,set<int> > kmerToReadSearch;
    map<size_t, set<size_t> > kmerToReadMap;
    //map<int,set<NodeID> >readToKmerSearch;
    //the first argument is node multiplicity,in the second pair the first argumument is confidense ratio for this multiplicity and the second argument is for correctntss ratio of gueess
    typedef pair<int, pair<double, double> > pair_k;
    map<NodeID, pair_k> nodesExpMult;
    typedef multimap<NodeID, pair_k>::iterator mapIterator;

    /**
     * Populates the table as required by the ReadCorrection procedure
     */
    void populateTable();
    /**
     * deletes the table again
     */
    void depopulateTable();
    /**
     * Find Kmer in the Kmernodetable
     */
    NodePosPair getNodePosPair(Kmer const &kmer) const;
    /**
     * Checks if the Kmer exists in the KmerNodeTable
     */
    bool kmerExistsInGraph(Kmer const &kmer) const;
    /**
     * get readLength
     */
    double getReadLength() const;
    /*
     *
     */
    void writeGraphFasta() const;
    /**
     *
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
    void graphPurification(string trueMultFilename,
                           const LibraryContainer& libraries);


};

#endif
