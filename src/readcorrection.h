/***************************************************************************
 *   Copyright (C) 2014 - 2016 Jan Fostier (jan.fostier@intec.ugent.be)    *
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

#ifndef READCORRECTION_H
#define READCORRECTION_H

#include "settings.h"
#include "graph.h"
#include "alignment.h"
#include "essaMEM-master/sparseSA.hpp"

#include <mutex>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class NodePosPair;

// ============================================================================
// ALIGNMENT METRICS CLASS
// ============================================================================

class AlignmentMetrics
{
private:
        size_t numReads;                // number of reads handled
        size_t numCorrReads;            // number of reads corrected
        size_t numCorrByMEM;            // number of times MEM procedure was used
        size_t numSubstitutions;        // number of substitutions made to the reads
        std::mutex metricMutex;         // mutex for merging metrics

public:
        /**
         * Default constructor
         */
        AlignmentMetrics() : numReads(0), numCorrReads(0), numCorrByMEM(0),
                numSubstitutions(0) {}

        /**
         * Update the statistics
         * @param corrected True if (part of) the read was corrected
         * @param corrByMEM True if MEM procedure was used
         * @param numSubstitutions_ Number of substitutions made to this read
         */
        void addObservation(bool corrected, bool corrByMEM,
                            size_t numSubstitutions_) {
                numReads++;
                if (corrected)
                        numCorrReads++;
                if (corrected && corrByMEM)
                        numCorrByMEM++;
                numSubstitutions += numSubstitutions_;
        }

        /**
         * Add other metrics (thread-safe)
         * @param metrics Metrics to add
         */
        void addMetrics(const AlignmentMetrics& rhs);

        /**
         * Get the number of reads
         * @return The number of reads
         */
        size_t getNumReads() const {
                return numReads;
        }

        /**
         * Get the number of corrected reads
         * @return The number of corrected reads
         */
        size_t getNumCorrReads() const {
                return numCorrReads;
        }

        /**
         * Get the number of corrected reads using the MEM procedure
         * @return The number of corrected reads using the MEM procedure
         */
        size_t getNumCorrByMEM() const {
                return numCorrByMEM;
        }

        /**
         * Get the number of substitutions made to the reads
         * @return The number of substitutions made to the reads
         */
        size_t getNumSubstitutions() const {
                return numSubstitutions;
        }

        /**
         * Output statistics to the stdout
         */
        void printStatistics() const;
};

// ============================================================================
// SEED CLASS
// ============================================================================

class Seed
{
public:
        std::vector<NodeID> nodeID;     // chain of node IDs
        size_t nodeFirst;               // offset within first node
        size_t readFirst;               // start position in read
        size_t readEnd;                 // end position in read

        Seed() : nodeID(0), nodeFirst(0), readFirst(0), readEnd(0) {}

        Seed(NodeID nodeID_, size_t nodeFirst_, size_t readFirst_, size_t readEnd_) :
             nodeFirst(nodeFirst_), readFirst(readFirst_), readEnd(readEnd_) {
                nodeID.push_back(nodeID_);
        }

        bool operator< (const Seed& rhs) const {
                if (nodeID.front() != rhs.nodeID.front())
                        return nodeID.front() < rhs.nodeID.front();
                return nodeFirst < rhs.nodeFirst;
        }

        /**
         * @param dbg Const-reference to the De Bruijn graph
         * @param npp Node Position pair vector (output)
         */
        void createNodePosition(const DBGraph& dbg,
                                std::vector<NodePosPair>& npp) const;

        static void mergeSeeds(const std::vector<Seed>& seeds,
                               std::vector<Seed>& mergedSeeds);
};

// ============================================================================
// DFSNode CLASS
// ============================================================================

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

// ============================================================================
// READ CORRECTION CLASS
// ============================================================================

class ReadCorrection
{
private:
        const DBGraph &dbg;
        const Settings &settings;
        AlignmentJan alignment;
        const sparseSA& sa;
        const std::vector<long>& startpos;

        /**
         * Get the marginal length of a string
         * @param str String under consideration
         */
        size_t getMarginalLength(const std::string& str) const {
                return str.length() + 1 - Kmer::getK();
        }

        /**
         * Find the node position pairs for a read
         * @param read Reference to the read
         * @param npp Vector of node position pairs
         */
        void findNPPSlow(const std::string& read, std::vector<NodePosPair>& npp);

        /**
         * Find the node position pairs for a read
         * @param read Reference to the read
         * @param npp Vector of node position pairs
         */
        void findNPPFast(const std::string& read, std::vector<NodePosPair>& npp);

        /**
         * Correct a specific read record
         * @param TODO
         */
        void correctRead(std::string& read, std::vector<NodePosPair> npp,
                         size_t& first, size_t& last);

        /**
         * Correct a specific read record
         * @param record Record to correct (input/output)
         * @param metric Alignment metric to update (input/output)
         */
        void correctRead(ReadRecord& record, AlignmentMetrics& metric);

        /**
         * Correct the records in one chunk
         * @param npp Vector of node position pairs
         * @param seeds Vector containining consistent seeds (output)
         */
        void extractSeeds(const std::vector<NodePosPair>& nppv,
                              std::vector<Seed>& seeds);

        /**
         * Find the largest consecutive run and use this as seed
         * @param read Reference to the read
         * @param npp Vector of node position pairs
         * @param first First position of the seed (output)
         * @param last Last position of the seed (output)
         */
        void extendSeed(std::string& read, std::vector<NodePosPair>& npp,
                        size_t& first, size_t& last);

        /**
         * Find the largest consecutive run and use this as seed
         * @param read Reference to the read
         * @param npp Vector of node position pairs
         * @param first First position of the seed
         * @param last Last position of the seed
         * @return false if no seed was found
         */
        void applyReadCorrection(std::string& read,
                              const std::vector<NodePosPair>& npp,
                              size_t first, size_t last);

        void recSearch(NodeID curr, string& read, vector<NodePosPair>& npp,
                       size_t currPos, size_t& counter, int score,
                       int& bestScore, size_t& seedLast);

        void revCompl(vector<NodePosPair>& npp);

        void findSeedKmer(const std::string& read,
                          std::vector<Seed>& seeds);

        /**
         * Find the node position pairs for a read using EssaMEM
         * @param read Reference to the read
         * @param seeds TODO
         */
        void findSeedMEM(const std::string& read, std::vector<Seed>& seeds);

        int correctRead(const std::string& read,
                        std::string& bestCorrectedRead,
                        const std::vector<Seed>& seeds);

public:
        /**
         * Default constructor
         * @param dbg_ Reference to the De Bruijn graph
         * @param settings_ Reference to the settings class
         */
        ReadCorrection(const DBGraph& dbg_, const Settings& settings_,
                          const sparseSA& sa_, const std::vector<long>& startpos_) :
                          dbg(dbg_), settings(settings_),
                          alignment(100, 2, 1, -1, -3), sa(sa_), startpos(startpos_) {}

        /**
         * Correct the records in one chunk
         * @param readChunk Chunk of records to correct
         * @param metric Alignment metric to update
         */
        void correctChunk(std::vector<ReadRecord>& readChunk,
                          AlignmentMetrics& metric);
};

// ============================================================================
// READ CORRECTION HANDLER CLASS
// ============================================================================

class ReadCorrectionHandler
{
private:
        DBGraph &dbg;
        const Settings &settings;
        sparseSA *sa;
        std::string reference;
        std::vector<long> startpos;

        void initEssaMEM();

        /**
         * Entry routine for worker thread
         * @param myID Unique threadID
         * @param libaries Library container with libraries to be corrected
         * @param metrics Alignment metrics accross threads
         */
        void workerThread(size_t myID, LibraryContainer& libraries,
                          AlignmentMetrics& metrics);

public:
        /**
         * Default constructor
         */
        ReadCorrectionHandler(DBGraph& g, const Settings& s);

        /**
         * Destructor
         */
        ~ReadCorrectionHandler();

        /**
         * Perform error correction in the libaries
         * @param libraries Library container with libraries to be corrected
         */
        void doErrorCorrection(LibraryContainer &libraries);
};

#endif
