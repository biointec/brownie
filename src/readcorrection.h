/***************************************************************************
 *   Copyright (C) 2015 Jan Fostier (jan.fostier@intec.ugent.be)           *
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

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class NodePosPair;


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

class ReadCorrectionJan
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
        void findNPPSlow(std::string& read, std::vector<NodePosPair>& npp);

        /**
         * Find the node position pairs for a read
         * @param read Reference to the read
         * @param npp Vector of node position pairs
         * @return True of at least one kmer hit was found, false otherwise
         */
        bool findNPPFast(std::string& read, std::vector<NodePosPair>& npp);

        /**
         * Find the node position pairs for a read using EssaMEM
         * @param read Reference to the read
         * @param npp Vector of node position pairs
         * @return True of at least one kmer hit was found, false otherwise
         */
        bool findNPPEssaMEM(std::string& read, std::vector<NodePosPair>& npp);

        /**
         * Correct a specific read record
         * @param TODO
         */
        void correctRead(std::string& read, std::vector<NodePosPair> npp,
                         size_t& first, size_t& last);

        /**
         * Correct a specific read record
         * @param record Record to correct (input/output)
         * @return true of the read was corrected, false otherwise
         */
        bool correctRead(ReadRecord& record);

        /**
         * Correct the records in one chunk
         * @param npp Vector of node position pairs
         * @param consistency Vector containining the consistency (output)
         */
        void checkConsistency(std::vector<NodePosPair>& npp,
                              vector<bool>& consistency);

        /**
         * Find the largest consecutive run and use this as seed
         * @param npp Vector of node position pairs
         * @param consistency Vector containining the consistency
         * @param seeds Vector containing start/end positions of seeds
         */
        void findSeed(const std::vector<NodePosPair>& npp,
                      const std::vector<bool>& consistency,
                      std::vector<std::pair<size_t, size_t> >& seeds);

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

public:
        /**
         * Default constructor
         * @param dbg_ Reference to the De Bruijn graph
         * @param settings_ Reference to the settings class
         */
        ReadCorrectionJan(const DBGraph& dbg_, const Settings& settings_,
                          const sparseSA& sa_, const std::vector<long>& startpos_) :
                          dbg(dbg_), settings(settings_),
                          alignment(100, 3, 1, -1, -3), sa(sa_), startpos(startpos_) {}

        /**
         * Correct the records in one chunk
         * @param readChunk Chunk of records to correct
         */
        void correctChunk(std::vector<ReadRecord>& readChunk);
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
         */
        void workerThread(size_t myID, LibraryContainer& libraries);

public:
        /**
         * Default constructor
         */
        ReadCorrectionHandler(DBGraph& g, const Settings& s);

        /**
         * Destructor
         */
        ~ReadCorrectionHandler() { delete sa; }

        /**
         * Perform error correction in the libaries
         * @param libraries Library container with libraries to be corrected
         */
        void doErrorCorrection(LibraryContainer &libraries);
};

#endif
