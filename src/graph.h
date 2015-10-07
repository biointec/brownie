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
#include "gaussval.h"
#include "observations.h"
//#include "anchorpath.h"
#include <cstdlib>
#include <deque>
#include <map>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include "sparseSA.hpp"
#include "fasta.hpp"
#include <queue>
typedef std::pair<NodePair, GaussVal> Observation;

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
    KmerNodeTable* table;                   // kmer node table

    /**
     * Parse a buffer of reads and store kmers in temporary buffers per thread
     * @param readBuffer Input read buffer
     */
    void parseReads(size_t thisThread,
                    std::vector<std::string>& readBuffer) const;

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
     * Estimate the node multiplicity based on the kmer occurrence and the
     * expected kmer occurce (mu)
     * @return True of at least one node has a new multiplicity
     */
    bool setNodeMultiplicity(const vector<vector<size_t> > &kmerFreq);

    /**
     * Based on preset node multiplicities and node kmer occurences,
     * estimate the expected kmer occurence (mu)
     * @return The expected average occurence
     */
    double estimateReadStartCoverage(const ReadLibrary& input,
                                     const vector<size_t> & kmerFreq) const;

    /**
     * Coverage-based arc filtering
     * @return True if at least one arc was removed
     */
    bool flowConservationBasedArcFiltering();

    /**
     * Filter all arcs to/from nodes with zero expected multiplicity
     * @return True if at least one arc was removed
     */
    bool arcFilteringBasedOnExpMultiplicity();

    /**
     * Filter all arcs to/from nodes with zero expected multiplicity
     * @return True if at least one arc was removed
     */
    bool arcFilteringBasedOnSignExpMultiplicity();

    /**
     * Starting from a seed node with coverage zero, find a cluster with
     * coverage zero
     * @param curr Current seed node
     * @param cluster Vector of nodeIDs belonging to that cluster
     * @return True of such cluster is found, false otherwise
     */
    bool findIsolatedCluster(SSNode& curr, std::vector<NodeID>& cluster);

    /**
     * Mark the anchor nodes
     */
    void setAnchorNodes();

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
     * Remove a stretch of nodes from the graph
     * @param lNode Leftmost node to remove
     * @param rNode Rightmost node to remove
     */
    void removeTip(SSNode &lNode, SSNode &rNode);

    /**
     * Increment the coverage of the arc that span two given kmers
     * @param left Left kmer reference
     * @param right Right kmer reference
     */
    void increaseCoverage(const NodeEndRef &left, const NodeEndRef &right);

    /**
     * Clip dead ends from the graph
     * @param hard Aggressive clipping or not
     * @return True if the graph was modified
     */
    bool clipTips(bool hard = false);

    /**
     * Check if a tip is eligible for clipping and clip if so
     * @param startNode SSNode that has no left arcs
     * @return True of the tip was clipped from the graph, false otherwise
     */
    bool clipTipFromNode(SSNode &startNode);

    /**
     * Check if a tip is eligible for clipping and clip if so (hard!)
     * @param startNode SSNode that has no left arcs
     * @return True of the tip was clipped from the graph, false otherwise
     */
    bool clipTipFromNodeHard(SSNode &startNode);

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

    bool checkReduction(NodeID lID, vector<NodeID> vmID, NodeID rID);


    void markPairedArcs(const std::vector<NodeID>& seq);

    bool reduction(NodeID lID, vector<NodeID> vmID, NodeID rID);
    bool covBasedReduction(NodeID lID, vector<NodeID> vmID, NodeID rID);

    // ====================================================================
    // VARIABLES
    // ====================================================================

    const Settings &settings;     // settings object

    DSNode *nodes;          // graph nodes
    Arc *arcs;              // graph arcs

    NodeID numNodes;        // number of nodes
    NodeID numArcs;         // number of arcs

    MapType mapType;

    map<NodePair, Observations> nodeDist;

    std::multimap<NodePair, GaussVal> signLinks;
    std::multimap<NodePair, GaussVal> nonSignLinks;
    std::multimap<NodePair, GaussVal> dubiousLinks;
    double readLength;//=100;
    double coverage;//=100;
    int kmerSize;//=31;

    double totalRAM;
    // SOLUTION COMP
    std::vector<string> reference;
    std::vector<void* > refST;
    std::map<NodeID, size_t> trueMult;



public:
    enum readCorrectionStatus {

        kmerfound,
        fullHealing,
        parHealing,
        kmerNotfound,
        graphIsMissing,
        anotherKmer

    };
    static const DBGraph* graph;
    double estimatedKmerCoverage;//=(coverage/readLength)*(readLength-kmerSize+1)/12;
    double estimatedMKmerCoverageSTD;//=sqrt( estimatedKmerCoverage);
    double minCertainVlueCov;
    double minSafeValueCov;
    double minRedLineValueCov;
    double certainVlueCov;
    double safeValueCov;
    double redLineValueCov;
    double cutOffvalue;
    size_t maxNodeSizeToDel;
    int updateCutOffValueRound;

    size_t n50;
    int numberOfValidNodes;
    size_t sizeOfGraph;
    double erroneousClusterMean;
    double correctClusterMean;
    /**
     * Careful concatenation, taking into account the estimated multiplicity
     */
    void initialize();
    bool mergeSingleNodes(bool force);

    bool bubbleDetection(int round);
    bool bubbleDetection();
    vector<pair<SSNode, SSNode> >  ExtractBubbles( SSNode rootNode);
    bool removeBubble(SSNode &prevFirstNode ,SSNode& extendFirstNode,size_t &TP,size_t &TN,size_t &FP,size_t &FN,size_t & numOfDel);
    bool removeNotSingleBublles(  SSNode &prevFirstNode ,SSNode& extendFirstNode, size_t &TP,size_t &TN,size_t &FP,size_t &FN,size_t & numOfDel);
    bool removeNode(SSNode &rootNode);
    void extractStatistic(int round);

    bool deleteUnreliableNodes( int round);
    bool deleteExtraRightLink(SSNode leftNode, int round);
    bool deleteExtraLeftLink(SSNode leftNode, int round);

    void canonicalBubble();
    double getMemoryUsedByProc();
    int  parseLine(char* line);
    bool deleteSuspiciousNodes();
    void adjustKmerReadRelation(SSNode targetNode,SSNode trustedNode);
    void correctReads();
    void makeBinaryFile();
    void updateReadFile(const map<long , string> &cacheTable  );
    string getReadByLine(long num, string fileName);
    void WriteReadsToTextFile();

    void test();

    double findDifference2(string a, string b);
    bool alignToKmer(string& read,const string& kmer, int readLine, int& startPoint, bool& kmerReverse);
    bool votingCorrection(vector<pair< pair<double,double>, bool> >  &tempArray,map<long , string > &cacheTable, string kmer);
    char findMax(int a, int c, int g, int t,int maxNum, int &num);


    char getReverseChar(char c);
    void findAvgCov();
    bool charDiff(string a, string b, int k);

    bool filterChimeric(int round);
    bool continueEdit(size_t& bigestN50, string& nodeFileName, string &arcFileName,string &metaDataFilename);
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
     */
    void compareToSolution();

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
     * Filter disconnected nodes with multiplicity zero
     */
    void filterIsolatedChaff();

    /**
     * Count the number of kmer occurences for each node
     */
    void countReadFrequency(const ReadLibrary& input);

    /**
     * Validate existing kmer frequencies
     */
    void validateKmerFrequency();

    /**
     * Set the node multiplicity based on an estimate for the coverage
     * based on expectation-maximization (EM)
     */
    void countNodeandArcFrequency(LibraryContainer &inputs);

    /**
     * Filter the graph, based on coverage
     */
    bool filterCoverage( float round);
    bool updateCutOffValue(int round);
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
    // PERM.CPP PUBLIC
    // ====================================================================

    /**
     * Perform graph reductions using paired-end reads
     */
    void perm();

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
     * Concatenate consecutive nodes
     */
    bool concatenation();

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
    SSNode getSSNode(NodeID nodeID) const {
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

    size_t getN50();

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
    void writeGraphExplicit() const;
};

#endif
