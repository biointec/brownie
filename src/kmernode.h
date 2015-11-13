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

#ifndef KMERNODETABLE_H
#define KMERNODETABLE_H

#include "global.h"
#include "tkmer.h"
#include "dsnode.h"
#include <google/sparse_hash_map>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class Settings;
class KmerNode;

// ============================================================================
// NODE EVENT CLASS
// ============================================================================

class NodePosPair : public std::pair<NodeID, PositionID>
{
public:
        /**
         * Default constructor
         */
        NodePosPair() :
                std::pair<NodeID, PositionID>(0, 0) {}

        /**
         * Default constructor
         * @param id Node identifier
         * @param pos Position identifier
         */
        NodePosPair(NodeID id, PositionID pos) :
                std::pair<NodeID, PositionID>(id, pos) {}

        /**
         * Get the node identifier
         * @return Node identifier
         */
        NodeID getNodeID() const {
                return first;
        }

        /**
         * Check whether the object points to valid position
         * @return True of false
         */
        bool isValid() const {
                return first != 0;
        }

        /**
         * Get the offset position
         * @return Offset position
         */
        PositionID getOffset() const {
                return second;
        }
};

// ============================================================================
// NODE EVENT CLASS
// ============================================================================

class NodeEvent {

private:
        NodeID target;          // target node
        PositionID lOffset;     // left offset position
        PositionID rOffset;     // right offset position
        bool duplication;       // duplication event or not
        size_t timeStamp;       // timeStamp of the event

public:
        /**
         * Default constructor
         */
        NodeEvent() : target(0), lOffset(0), rOffset(0),
                      duplication(false), timeStamp(0) {}

        /**
         * Default constructor
         * @param target Node identifier for the new node
         * @param lOffset Left offset in the new node
         * @param rOffset Right offset in the new node
         * @param duplication True of the node is to be duplicated
         * @param timeStamp Timestamp of the event
         */
        NodeEvent(NodeID target_, PositionID lOffset_, PositionID rOffset_,
                  bool duplication_, size_t timeStamp_) : target(target_),
                  lOffset(lOffset_), rOffset(rOffset_),
                  duplication(duplication_), timeStamp(timeStamp_) {};

        /**
         * Get the identifier of the target node
         * @return Identifier of the target node
         */
        NodeID getTargetID() const {
                return target;
        }

        /**
         * Get the relative left position in the new node
         * @return Position in the new node
         */
        PositionID getLeftOffset() const {
                return lOffset;
        }

        /**
         * Get the relative right position in the new node
         * @return Position in the new node
         */
        PositionID getRightOffset() const {
                return rOffset;
        }

        /**
         * True if the node is deleted
         * @return True of false
         */
        bool isDuplication() const {
                return duplication;
        }

        /**
         * Get the timestamp of the event
         * @return Timestamp of the event
         */
        size_t getTimeStamp() const {
                return timeStamp;
        }
};

// ============================================================================
// TYPEDEFS
// ============================================================================

typedef google::sparse_hash_map<Kmer, KmerNode, KmerHash> GoogleKmerNodeTable;

// shortcut notation for a const iterator
typedef GoogleKmerNodeTable::const_iterator KmerNodeIt;

// shortcut notation for a <Key, Data> pair
typedef std::pair<Kmer, KmerNode> KmerNodeValue;

// ============================================================================
// KMER NODE CLASS
// ============================================================================

class KmerNode {

private:
        NodeID nodeID;
        PositionID pos;

public:
        /**
         * Default constructor
         */
        KmerNode() : nodeID(0), pos(0) {}

        /**
         * Constructor
         * @param nodeID Node identifier
         * @param pos Position offset of the kmer within the node
         */
        KmerNode(NodeID nodeID_, PositionID pos_) :
                nodeID(nodeID_), pos(pos_) {}

        /**
         * Obtain node identifier
         * @return node identifier
         */
        NodeID getNodeID() const {
                return nodeID;
        }

        /**
         * Obtain position
         * @return node position
         */
        PositionID getPosition() const {
                return pos;
        }
};

// ============================================================================
// KMER NODE REF
// ============================================================================

class KmerNodeRef : public std::pair<KmerNodeIt, bool> {

private:
        static const DSNode *nodes;

public:
        /**
         * Default constructor
         */
        KmerNodeRef() {
                // make sure the static member is set before attempting to
                // instantiate KmerNodeRef
                assert(nodes != NULL);
        };

        /**
         * Constructor
         * @param it Iterator to the table
         * @param reverse True if the iterator points to the reverse complement
         */
        KmerNodeRef(KmerNodeIt it, bool reverse) :
                std::pair<KmerNodeIt, bool>(it, reverse) {}

        /**
         * Get the nodeID this kmer belongs to
         * @return The node identifier
         */
        NodeID getNodeID() {
                if (second)
                        return -first->second.getNodeID();
                else
                        return first->second.getNodeID();
        }

        /**
         * Get the nodeID this kmer belongs to
         * @return The node identifier
         */
        PositionID getPosition() {
                if (second)
                        return nodes[abs(first->second.getNodeID())].getLength() -
                               first->second.getPosition() - Kmer::getK();
                else
                        return first->second.getPosition();
        }

        /**
         * Set the static nodes member
         * @param nodes Pointer to the nodes
         */
        static void setNodes(const DSNode *nodes_) {
                nodes = nodes_;
        }
};

// ============================================================================
// KMER NODE TABLE CLASS
// ============================================================================

class KmerNodeTable {

private:

        /**
         * Insert a kmer in the table
         * @param kmer Kmer to insert
         * @param nodeID Node identifier
         * @param pos Position in the node
         * @param node Double stranded node reference
         * @return KmerRef containing iterator to the kmer and reversed flag
         */
        KmerNodeRef insert(const Kmer &kmer, NodeID id,
                           PositionID pos, const DSNode &node);

        const Settings &settings;               // reference to the settings
        NodeID numNodes;                        // number of nodes
        GoogleKmerNodeTable *table;             // actual table
        std::vector<NodeEvent> *remapInfo;      // remapping of nodes
        size_t timeStamp;                       // current timestamp

public:

        /**
         * Default constructor
         * @param settings Reference to the settings object
         */
        KmerNodeTable(const Settings& settings, NodeID numNodes);

        /**
         * Destructor
         */
        ~KmerNodeTable();

        /**
         * Create a kmer node table
         * @param nodes Pointer to the double stranded nodes
         */
        void populateTable(const DSNode *nodes);

        /**
         * Find a kmer in the graph
         * @param kmer Kmer to look for
         * @return A vector containing all nodes and their position in that node
         */
        void find(const Kmer& kmer, std::vector<NodePosPair>& npp) const;

        /**
         * Recursively find all occurrences of a kmer in the graph
         * @param curr Current node, position pair
         * @param currTimeStamp Current timestamp
         * @param npp Vector containg all solutions
         */
        void recFindInTable(NodePosPair curr, size_t currTimeStamp,
                            std::vector<NodePosPair>& npp) const;

        /**
         * Find a kmer in the table
         * @param kmer Kmer to look for
         * @return The node, position pair of that kmer
         */
        NodePosPair find(const Kmer& kmer) const;

        /**
         * Merge left node to right node
         * @param leftID Identifier for the left node
         * @param righID Identifier for the right node
         * @param leftSize Marginal size of the left node
         * @param rightSize Marginal size of the right node
         * @param deleteLeft Delete left node after merging?
         */
        void mergeLeftToRight(NodeID leftID, NodeID rightID,
                              size_t leftSize, size_t rightSize,
                              bool deleteLeft);

        /**
         * Create a kmer node table
         * @param nodes Pointer to the double stranded nodes
         * @param numNodes Number of nodes
         */
        void sanityCheck(const DSNode *nodes, NodeID numNodes) const;
};

#endif
