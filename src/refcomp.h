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

#ifndef REFCOMP_H
#define REFCOMP_H

#include <string>
#include <vector>

#include "nodechain.h"

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class DBGraph;

// ============================================================================
// BREAKPOINT CLASS
// ============================================================================

class BreakPoint {

private:
        size_t refID;
        size_t begin;
        size_t end;

public:
        /**
         * Default constructor
         * @param refID_ Reference identifier
         * @param begin_ Start position of the breakpoint
         */
        BreakPoint(size_t refID_, size_t begin_) :
                refID(refID_), begin(begin_) { end = begin_ + 1; };


        /**
         * Extend the current breakpoint by one
         */
        void extendBreakPoint() {
                end++;
        }

        /**
         * Operator << overloading
         * @param out Output stream
         * @param bp Breakpoint to display
         * @return Output stream
         */
        friend std::ostream &operator<<(std::ostream &out, const BreakPoint &bp);

        /**
         * Get the reference ID of the breakpoint
         * @return The reference ID of the breakpoint
         */
        size_t getRefID() const {
                return refID;
        }

        /**
         * Get the begin point of the breakpoint
         * @return The begin point of the breakpoint
         */
        size_t getBegin() const {
                return begin;
        }

        /**
         * Get the end point of the breakpoint
         * @return The end point of the breakpoint
         */
        size_t getEnd() const {
                return end;
        }
};

// ============================================================================
// REFERENCE SEGMENT CLASS
// ============================================================================

enum class RefSegmentType { CONTIG, PARALLEL, BREAK, DELETION, INSERTION, UNKNOWN };

class RefSegment {

private:
        size_t refID;           // reference ID of this segment
        size_t refBegin;        // first position in the reference
        size_t refEnd;          // last position in the reference
        NodeID startNodeID;     // node ID of the first node
        NodeID endNodeID;       // node ID of the last node
        size_t startNodeBegin;  // first position in the start node
        size_t endNodeEnd;      // last position in the end node
        RefSegmentType type;    // type of segment

public:
        /**
         * @brief Default constructor
         * @param refID_ Reference identifier
         */
        RefSegment(size_t refID_) : refID(refID_) {}


        void setBeginPos(size_t refBegin_, NodeID startNodeID_, size_t startNodeBegin_) {
                refBegin = refBegin_;
                startNodeID = startNodeID_;
                startNodeBegin = startNodeBegin_;
        }

        void setEndPos(size_t refEnd_, NodeID endNodeID_, size_t endNodeEnd_) {
                refEnd = refEnd_;
                endNodeID = endNodeID_;
                endNodeEnd = endNodeEnd_;
        }

        void setType(RefSegmentType type_) {
                type = type_;
        }

        RefSegmentType getType() const {
                return type;
        }

        /**
         * Operator << overloading
         * @param out Output stream
         * @param bp Breakpoint to display
         * @return Output stream
         */
        friend std::ostream &operator<<(std::ostream &out, const RefSegment &bp);

        /**
         * Get the reference ID
         * @return The reference ID
         */
        size_t getRefID() const {
                return refID;
        }

        /**
         * @brief Get the reference begin point
         * @return The reference begin point
         */
        size_t getRefBegin() const {
                return refBegin;
        }

        /**
         * @brief Get the start node ID
         * @return The start node ID
         */
        NodeID getStartNodeID() const {
                return startNodeID;
        }

        /**
         * @brief Get the begin position on the start node
         * @return The begin position on the start node
         */
        size_t getStartNodeBegin() const {
                return startNodeBegin;
        }

        /**
         * @brief Get the reference end point
         * @return The reference end point
         */
        size_t getRefEnd() const {
                return refEnd;
        }

        /**
         * @brief Get the end node ID
         * @return The end node ID
         */
        NodeID getEndNodeID() const {
                return endNodeID;
        }

        /**
         * @brief Get the end position on the end node
         * @return The end position on the end node
         */
        size_t getEndNodeEnd() const {
                return endNodeEnd;
        }
};

// ============================================================================
// REFERENCE COMPARISON CLASS
// ============================================================================

class RefComp {

private:
        std::vector<std::string> refSeq_v;      // reference sequences

        /**
         * @brief Align reference sequences to a de Bruijn graph
         * @param dbg de Bruijn graph
         * @param refSegment_v Reference segments (output)
         */
        void alignReference(const DBGraph& dbg,
                            std::vector<RefSegment>& refSegment_v) const;

        /**
         * @brief Annotate the breakpoints
         * @param dbg de Bruijn graph
         * @param refSegment_v Reference segments (input/output)
         */
        void annotateBreakpoints(const DBGraph& dbg,
                                 std::vector<RefSegment>& refSegment_v) const;

        /**
         * @brief Output segments to the screen
         * @param segment_v Segments to output
         */
        void printSegments(const std::vector<RefSegment>& segment_v) const;

public:
        /**
         * Default constructor
         * @param refFilename File name of the reference genome
         */
        RefComp(const std::string& refFilename);

        /**
         * @brief Get the size of the reference contigs
         * @return The size of the reference contigs
         */
        size_t getSize() const;

        /**
         * @brief Get the number of contigs
         * @return The number of contigs
         */
        size_t getNumContigs() const { return refSeq_v.size(); }

        /**
         * @brief Align reference to a de Bruijn graph
         * @param dbg de Brijn graph
         */
        void findAlignableContigs(const DBGraph& dbg);

        /**
         * @brief Validate a de Bruijn graph
         * @param dbg A const-ref to the de Bruijn graph
         */
        void validateGraph(const DBGraph& dbg, size_t minContigSize = 1);

        /**
         * Get the true node chains from the reference sequence
         * @param dbg A const-ref to the de Bruijn graph
         * @param nodeChain Node chains (output)
         */
        void getTrueNodeChain(const DBGraph& dbg,
                              std::vector<NodeChain>& nodeChain);

        /**
         * extract the set of nodes which contains the kmers of a given sequence
         * @param dbg A const-ref to the de Bruijn graph
         * @param refSeq given sequence
         * @param nodeSet set of ndoes (output)
         * @return the size of the subgraph
         */
        size_t extractNodeSetbySequence(const DBGraph& dbg,
                                      std::string &refSeq, std::set<int> &nodeSet);
        /**
         * extract the the subgraph nearbye given sequence in a fasta file
         * @param dbg A const-ref to the de Bruijn graph
         * @param breakpointFileName the fasta file contains sequences which are supposed to be in breakpoint area
         * @param tempDir the temporary directory to save Cytoscape files
         */
        void extractBreakpointSubgraph(const DBGraph& dbg,
                                                std::string breakpointFileName, std::string tempDir);
        /**
         * Calculate the true node multiplicity
         * @param dbg A const-ref to the de Bruijn graph
         * @param multiplicity Multiplicity vector
         */
        void getNodeMultiplicity(const DBGraph& dbg,
                                 std::vector<size_t>& multiplicity);

        /**
         * Get the true node chains from the given file as a reference
         * @param dbg A const-ref to the de Bruijn graph
         * @param breakpointFileName the fasta file contains sequences which are supposed to be in breakpoint area
         * @return number of breakpoint
         */
        size_t findbreakpoints(const DBGraph& dbg,std::string breakpointFileName);
};

#endif
