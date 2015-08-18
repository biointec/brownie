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

#ifndef SSNODE_H
#define SSNODE_H

#include "global.h"
#include "dsnode.h"
#include "nucleotide.h"
#include "tkmer.h"

// ============================================================================
// NODE CLASS
// ============================================================================

class SSNode {

private:
        static DSNode *nodes;   // pointer to the double stranded nodes

        NodeID nodeID;          // identiffier of the node
        DSNode *dsNode;         // reference to the double stranded node

public:
        /**
         * Default constructor
         */
        SSNode() : nodeID(0), dsNode(NULL) {}

        /**
         * Constructor
         * @param dsNode Double stranded node
         * @param ID Unique identifier of the node
         */
        SSNode(DSNode* dsNode, NodeID nodeID) : nodeID(nodeID), dsNode(dsNode) {
                assert(nodeID != 0);
        }

        void setVisited(bool target) {
                dsNode->setVisited(target);
        }

        bool isVisited() const {
                return dsNode->isVisited();
        }

        /**
         * Constructor
         * @param it Iterator pointing to this node
         */
        SSNode(const ArcIt& it) : nodeID(it->getNodeID()), dsNode(nodes + abs(nodeID)) {
                assert(nodeID != 0);
        }

        /**
         * Check whether this node is an anchor node or not
         * @return True of false
         */
        bool isAnchor() const {
                return dsNode->isAnchor();
        }

        /**
         * Check whether this node is an anchor node or not
         * @return True of false
         */
        void setAnchor(bool value) {
                dsNode->setAnchor(value);
        }

        /**
         * Set the loop flag
         * @param isLoop True of false
         */
        void setLoop(bool isLoop) {
                dsNode->setLoop(isLoop);
        }

        /**
         * Check whether the node is a loop or not
         * @return True of false
         */
        bool isLoop() const {
                return dsNode->isLoop();
        }

        uint8_t getFlag() const {
                return dsNode->getFlag();
        }

        /**
         * Get the expected multiplicity
         * @return The expected multiplicity
         */
        double getExpMult() const {
                return dsNode->getExpMult();
        }

        /**
         * Set the expected multiplicity
         * @param target Target multiplicity
         */
        void setExpMult(double target) {
                dsNode->setExpMult(target);
        }

        /**
         * Set the read start coverage
         * @param target The target read start coverage
         */
        void setReadStartCov(Coverage target) {
                dsNode->setReadStartCov(target);
        }

        /**
         * Get the read start coverage
         * @return The read start coverage
         */
        Coverage getReadStartCov() const {
                return dsNode->getReadStartCov();
        }

        /**
         * Atomically increment the read start coverage
         */
        void incReadStartCov() {
                dsNode->incReadStartCov();
        }

        /**
         * Set the kmer coverage
         * @param target The kmer coverage
         */
        void setKmerCov(Coverage target) {
                dsNode->setKmerCov(target);
        }

        /**
         * Get the kmer coverage
         * @return The kmer coverage
         */
        Coverage getKmerCov() const {
                return dsNode->getKmerCov();
        }

        /**
         * Atomically increment the kmer coverage
         */
        void incKmerCov() {
                dsNode->incKmerCov();
        }

        /**
         * Get the multiplicity, rounded to the closest integer
         * @return The multiplicity
         */
        size_t getRoundMult() const {
                return dsNode->getRoundMult();
        }

        /**
         * Get the low side estimation of the multiplicity
         * @return The low side estimation of the multiplicity
         */
        size_t getLoExpMult() const {
                return dsNode->getLoExpMult();
        }

        /**
         * Get the high side estimation of the multiplicity
         * @return The high side estimation of the multiplicity
         */
        size_t getHiExpMult() const {
                return dsNode->getHiExpMult();
        }

        /**
         * Check whether the multiplicity estimate is dubious
         * @return True of false
         */
        bool multIsDubious() const {
                return dsNode->multIsDubious();
        }

        /**
         * Invalidate this node
         */
        void invalidate() {
                dsNode->invalidate();
        }

        /**
         * Check if a node is invalidated
         * @return True or false
         */
        bool isValid() {
                if (nodeID == 0)
                        return false;
                return dsNode->isValid();
        }

        /**
         * Get the identifier of this node
         */
        NodeID getNodeID() const {
                return nodeID;
        }

        /**
         * Operator '==' overloading
         * @param rhs Right hand side SSNode
         * @return True if they're equal
         */
        bool operator==(const SSNode &rhs) const {
                if (dsNode != rhs.dsNode)
                        return false;
                return (nodeID == rhs.nodeID);
        }

        /**
         * Operator '!=' overloading
         * @param rhs Right hand side kmer
         * @return True if they're different
         */
        bool operator!=(const SSNode &rhs) const {
                return !(*this == rhs);
        }

        /**
         * Get the length of the node (# nucleotides in DNA string)
         * @return The length of the node
         */
        NodeLength getLength() const {
                return dsNode->getLength();
        }

        /**
         * The marginal length == length - k + 1
         * @return The marginal length of the node
         */
        NodeLength getMarginalLength() const {
                return dsNode->getMarginalLength();
        }

        /**
         * Get the number of left arcs
         * @return The number of left arcs
         */
        uint8_t getNumLeftArcs() const {
                return (nodeID > 0) ?
                        dsNode->getNumLeftArcs() : dsNode->getNumRightArcs();
        }

        /**
         * Get the number of right arcs
         * @return The number of right arcs
         */
        uint8_t getNumRightArcs() const {
                return (nodeID > 0) ?
                        dsNode->getNumRightArcs() : dsNode->getNumLeftArcs();
        }

        /**
         * Get the number of left arcs
         * @param The number of left arcs
         */
        void setNumLeftArcs(uint8_t numArcs) const {
                if (nodeID > 0)
                        dsNode->setNumLeftArcs(numArcs);
                else
                        dsNode->setNumRightArcs(numArcs);
        }

        /**
         * Get the number of right arcs
         * @param The number of right arcs
         */
        void setNumRightArcs(uint8_t numArcs) const {
                if (nodeID > 0)
                        dsNode->setNumRightArcs(numArcs);
                else
                        dsNode->setNumLeftArcs(numArcs);
        }

        void swapRightArcsSign() {
                if (nodeID > 0)
                        dsNode->swapRightArcsSign();
                else
                        dsNode->swapLeftArcsSign();
        }

        void copyRightArcs(SSNode &source) {
                int numRightArcs = source.getNumRightArcs();
                source.setNumRightArcs(0);
                setNumRightArcs(numRightArcs);
                setFirstRightArcID(source.getFirstRightArcID());
                // if they have opposite sign
                if ((nodeID < 0) != (source.getNodeID() < 0))
                        swapRightArcsSign();
        }

        /**
         * Get an iterator pointing to the first left arc
         * @return an iterator pointing to the first left arc
         */
        ArcIt leftBegin() const {
                return (nodeID > 0) ?
                        dsNode->leftBegin(false) : dsNode->rightBegin(true);
        }

        /**
         * Get an iterator pointing past the last left arc
         * @retur nan iterator pointing to the last left arc
         */
        ArcIt leftEnd() const {
                return (nodeID > 0) ?
                        dsNode->leftEnd(false) : dsNode->rightEnd(true);
        }

        /**
         * Get an iterator pointing to the first right arc
         * @return an iterator pointing to the first left arc
         */
        ArcIt rightBegin() const {
                return (nodeID > 0) ?
                        dsNode->rightBegin(false) : dsNode->leftBegin(true);
        }

        /**
         * Get an iterator pointing past the last right arc
         * @retur nan iterator pointing to the last right arc
         */
        ArcIt rightEnd() const {
                return (nodeID > 0) ?
                        dsNode->rightEnd(false) : dsNode->leftEnd(true);
        }

        /**
         * Delete the left arcs
         */
        void deleteAllLeftArcs() {
                if (nodeID > 0)
                        dsNode->deleteLeftArcs();
                else
                        dsNode->deleteRightArcs();
        }

        /**
         * Delete the right arcs
         */
        void deleteAllRightArcs() {
                if (nodeID > 0)
                        dsNode->deleteRightArcs();
                else
                        dsNode->deleteLeftArcs();
        }

        /**
         * Get the identifier for the left right arc
         * @return The identifier for the left right arc
         */
        ArcID getFirstLeftArcID() {
                if (nodeID > 0)
                        return dsNode->getFirstRightArcID();
                return dsNode->getFirstLeftArcID();
        }

        /**
         * Get the identifier for the first right arc
         * @return The identifier for the first right arc
         */
        ArcID getFirstRightArcID() {
                if (nodeID > 0)
                        return dsNode->getFirstRightArcID();
                return dsNode->getFirstLeftArcID();
        }

        /**
         * Set the identifier for the left right arc
         * @param The identifier for the left right arc
         */
        void setFirstLeftArcID(ArcID target) {
                if (nodeID > 0)
                        return dsNode->setFirstRightArcID(target);
                return dsNode->setFirstLeftArcID(target);
        }

        /**
         * Set the identifier for the first right arc
         * @param The identifier for the first right arc
         */
        void setFirstRightArcID(ArcID target) {
                if (nodeID > 0)
                        return dsNode->setFirstRightArcID(target);
                return dsNode->setFirstLeftArcID(target);
        }

        /**
         * Delete a specific left arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteLeftArc(NodeID targetID) {
                if (nodeID > 0)
                        return dsNode->deleteLeftArc(targetID);
                return dsNode->deleteRightArc(-targetID);
        }

        /**
         * Delete a specific right arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteRightArc(NodeID targetID) {
                if (nodeID > 0)
                        return dsNode->deleteRightArc(targetID);
                return dsNode->deleteLeftArc(-targetID);
        }

        /**
         * Get a specific left arc
         * @param targetID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* getLeftArc(NodeID targetID) {
                if (nodeID > 0)
                        return dsNode->getLeftArc(targetID);
                return dsNode->getRightArc(-targetID);
        }

        /**
         * Get a specific right arc
         * @param targetID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* getRightArc(NodeID targetID) {
                if (nodeID > 0)
                        return dsNode->getRightArc(targetID);
                return dsNode->getLeftArc(-targetID);
        }

        /**
         * Get a specific left arc
         * @param nucleotide Nucleotide under consideration
         * @return Pointer to the specific arc, NULL if not found
         */
        NodeID getLeftArc(char nucleotide) {
                for (ArcIt it = leftBegin(); it != leftEnd(); it++) {
                       // cout << "This is node " << getNodeID() << " left arc to " << it->getNodeID() << " with right charachter " << SSNode(it).peekNucleotideRight() << endl;
                        if (SSNode(it).peekNucleotideMarginalRight() == nucleotide)
                                return it->getNodeID();
                }

                return 0;
        }

        /**
         * Get a specific right arc
         * @param nucleotide Nucleotide under consideration
         * @return Pointer to the specific arc, NULL if not found
         */
        NodeID getRightArc(char nucleotide) {
                for (ArcIt it = rightBegin(); it != rightEnd(); it++)
                        if (SSNode(it).peekNucleotideMarginalLeft() == nucleotide)
                                return it->getNodeID();

                return 0;
        }

        /**
         * Replace a left arc with a different one
         * @param origID Identifier for the original target node
         * @param newID Identifier for the new target node
         */
        void replaceLeftArc(NodeID origID, NodeID newID) {
                if (nodeID > 0) {
                        if (dsNode->getLeftArc(origID) == NULL)
                                std::cout << "Paniek ! " << std::endl;
                } else
                        if (dsNode->getRightArc(-origID) == NULL)
                                std::cout << "Paniek ! " << std::endl;
                if (nodeID > 0)
                        dsNode->getLeftArc(origID)->setNodeID(newID);
                else
                        dsNode->getRightArc(-origID)->setNodeID(-newID);
        }

        /**
         * Replace a right arc with a different one
         * @param origID Identifier for the original target node
         * @param newID Identifier for the new target node
         */
        void replaceRightArc(NodeID origID, NodeID newID) {
                if (nodeID > 0)
                        dsNode->getRightArc(origID)->setNodeID(newID);
                else
                        dsNode->getLeftArc(-origID)->setNodeID(-newID);
        }

        /**
         * Get the sequence of this node
         * @return stl string containing the sequence
         */
        std::string getSequence() const {
                std::string seq = dsNode->getSequence();
                if (nodeID < 0)
                        Nucleotide::revCompl(seq);

                return seq;
        }

        /**
         * Get a nucleotide at a specified position ('-' for out-of bounds)
         * @param pos Position in the sequence
         * @return Nucleotide at specified position
         */
        char getNucleotide(PositionID pos) const {
                // check for out-of-bounds
                if (pos >= getLength())
                        return '-';
                if (nodeID < 0)
                        return Nucleotide::getComplement(dsNode->getNucleotide(getLength() - pos - 1));
                else
                        return dsNode->getNucleotide(pos);
        }

        /**
         * Set the sequence of this node
         * @param str String containing only 'A', 'C', 'G' and 'T'
         */
        void setSequence(const std::string& str) {
                if (nodeID > 0)
                        dsNode->setSequence(str);
                else
                        dsNode->setSequence(Nucleotide::getRevCompl(str));
        }

        /**
         * Get the left kmer of this node
         * @return The left kmer of this node
         */
        Kmer getLeftKmer() {
                std::string seq = getSequence();
                return Kmer(seq);
        }

        /**
         * Get the right kmer of this node
         * @return The right kmer of this node
         */
        Kmer getRightKmer() {
                std::string seq = getSequence();
                return Kmer(seq, seq.size() - Kmer::getK());
        }

        /**
         * Get the leftmost nucleotide of this node
         * @return The leftmost nucleotide
         */
        char peekNucleotideLeft() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideLeft();
                return Nucleotide::getComplement(dsNode->peekNucleotideRight());
        }

        /**
         * Get the rightmost nucleotide of this node
         * @return The rightmost nucleotide
         */
        char peekNucleotideRight() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideRight();
                return Nucleotide::getComplement(dsNode->peekNucleotideLeft());
        }

        /**
         * Get the nucleotide at position k - 1
         * @return The nucleotide at position k - 1
         */
        char peekNucleotideMarginalLeft() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideMarginalLeft();
                return Nucleotide::getComplement(dsNode->peekNucleotideMarginalRight());
        }

        /**
         * Get the nucleotide at position size - k
         * @return The nucleotide at position size - k
         */
        char peekNucleotideMarginalRight() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideMarginalRight();
                return Nucleotide::getComplement(dsNode->peekNucleotideMarginalLeft());
        }

        void inheritRightArcs(SSNode& target) {
                // make sure the node has currently no right arcs
                assert(getNumRightArcs() == 0);

                // update the arc information for the connected nodes
                for (ArcIt it = target.rightBegin(); it != target.rightEnd(); it++) {
                        SSNode rightNode(it);
                        rightNode.replaceLeftArc(target.getNodeID(), nodeID);
                }

                // copy the arcs
                copyRightArcs(target);
        }

        /**
         * Set the static node pointer
         * @param nodes Pointer to the double stranded nodes
         */
        static void setNodePointer(DSNode *nodes_) {
                nodes = nodes_;
        }
};

#endif
