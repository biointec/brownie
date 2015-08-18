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

#ifndef NODENBH_H
#define NODENBH_H

#include "nodeposcoll.h"

#include <map>
#include <set>

// ============================================================================
// TYPEDEFS
// ============================================================================

typedef std::map<NodeID, NodePosColl>::iterator NodeNBHIt;
typedef std::map<NodeID, NodePosColl>::const_iterator CNodeNBHIt;
typedef std::map<double, NodeID>::iterator SortedNodesIt;

class ReadLibrary;
class MappedRead;
class Observations;
class DBGraph;
class Scaffold;

// ============================================================================
// NODE NEIGHBORHOOD CLASS
// ============================================================================

class NodeNBH
{
private:
        NodeID srcID;           // source node identifier
        double srcPos;          // source position
        GaussVal RoIB, RoIE;    // region of interest begin/end

        std::map<NodeID, NodePosColl> NPCv;    // node positions collections

        /**
         * Compute the normal distribution
         * @param x Point in which to evaluate the normal distribution
         * @param mu Average of the normal distribution
         * @param sigma Standard deviation of the normal distribution
         */
        double normal(double x, double mu, double sigma) const;

        /**
         * Compute the logarithm of the normal distribution
         * @param x Point in which to evaluate the lognormal distribution
         * @param mu Average of the lognormal distribution
         * @param sigma Standard deviation of the lognormal distribution
         */
        double logNormal(double x, double mu, double sigma) const;

        /**
         *
         *
         */
        void correctObs(const Input& input, bool PER, double lA, double lB,
                        std::vector<NodePos>& corrPosv) const;

        /**
         * Apply a mixture model and EM on the samples
         * @param order Number of components in the mixture model
         * @param std Standard deviation of the components
         * @param mu Average of the components (output)
         * @param N Number of observation per component (output)
         */
        void mixtureModel(const std::vector<double>& obs,
                          size_t order, double std,
                          std::vector<NodePos>& positions) const;

        /**
         * Apply a mixture model and EM on the samples
         * @param std Standard deviation of the components
         * @param mu Average of the components (output)
         * @param N Number of observations per component (output)
         * @param nodeSize Marginal length of the node
         */
        void mixtureModel(const std::vector<double>& obs, const Input& input, bool PER,
                          double lA, double lB,
                          std::vector<NodePos>& nodePosv) const;

        /**
         * Does a given node overlap with the region of interest
         * @param nodePos Position of that node (Gaussian value)
         * @param RoDE Region of data endpoint
         * @return False if they significantly do not overlap
         */
        bool nodeInRoI(const GaussVal& nodePos, double nodeLen) const;

        /**
         * Create a node neighborhood based on individual reads
         * @param input Library input
         * @param mRead Mapped reads
         * @param graph DB graph
         */
        void createNBHfromReads(const Input& input,
                                const std::vector<MappedRead>& mReads,
                                const DBGraph& graph);

        /**
         * Create a node neighborhood based on paired-end reads
         * @param input Library input
         * @param mRead Mapped reads
         * @param graph DB graph
         */
        void createNBHfromPER(const Input& input,
                              const std::vector<MappedRead>& mReads,
                              const DBGraph& graph);

        void setExpMult(const Input& input, bool PER, double lA, double lB,
                        std::vector<NodePos>& positions) const;

        /**
         * Add a node position to the neighborhood
         * @param dstID Destination identifier
         * @param pos A vector of node positions
         */
        void addPositions(NodeID dstID, const std::vector<NodePos>& pos);

public:
        /**
         * Constructor
         * @param srcID Source node identifier
         * @param scoreBegin Score start point, relative to srcID
         * @param scoreEnd Score end point, relative to endID
         * @param insertLength Insert size which was used
         * @param insertStd Insert std which was used
         */
        NodeNBH(NodeID srcID_, const GaussVal& SB_, const GaussVal& SE_) :
                srcID(srcID_), srcPos(0.0), RoIB(SB_), RoIE(SE_) {};

        /**
         * Create a node neighborhood
         * @param input Library input
         * @param mRead Mapped reads
         * @param graph DB graph
         */
        void createNBH(const Input& input,
                       const std::vector<MappedRead>& mReads,
                       const DBGraph& graph);

        /**
         * Operator [] overloading
         * @param ns Node specification
         * @return A reference to the node position
         */
        NodePos& operator[](const NodeSpec& ns) {
                return NPCv[ns.first].getPositions()[ns.second];
        }

        /**
         * Get the source node
         * @return Source identifier
         */
        NodeID getSrcID() const {
                return srcID;
        }

        /**
         * Get the begin score position
         * @return score begin
         */
        GaussVal getSB() const {
                return RoIB;
        }

        /**
         * Get the end score position
         * @return score begin
         */
        GaussVal getSE() const {
                return RoIE;
        }

        /**
         * Add paired node information to a specific destination node
         * @param dstID Identifier for the destination node
         * @param pnData Paired node information to the dstID node
         */
        void addNode(NodeID dstID, const NodePosColl& npColl) {
                assert(NPCv.find(dstID) == NPCv.end());
                NPCv[dstID] = npColl;
        }

        /**
         * Given the overlap of a node with the score window
         * @param position Tentative position for the node
         * @param size Size of the node
         * @return The overlap of the node with the score window
         */
        double getCoveredSize(double position, double size) const;

        /**
         * Get an iterator to the first paired node
         * @return An iterator to the first paired node
         */
        NodeNBHIt begin() {
                return NPCv.begin();
        }

        /**
         * Get an iterator past the first paired node
         * @return An iterator past the first paired node
         */
        NodeNBHIt end() {
                return NPCv.end();
        }

        /**
         * Get an iterator to the first paired node
         * @return An iterator to the first paired node
         */
        CNodeNBHIt begin() const {
                return NPCv.begin();
        }

        /**
         * Get an iterator past the first paired node
         * @return An iterator past the first paired node
         */
        CNodeNBHIt end() const {
                return NPCv.end();
        }

        /**
         * Get an iterator to a specific paired node
         * @param nodeID Node identifier
         * @return An iterator to a specific paired node
         */
        NodeNBHIt find(NodeID nodeID) {
                return NPCv.find(nodeID);
        }

        /**
         * Clear all paired nodes
         */
        void clearNodes() {
                NPCv.clear();
        }

        /**
         * Shift the positions by an offset
         * @param gv Offset specification
         */
        void shiftPositions(const GaussVal& gv);

        /**
         * Find unique nodes
         * @param uniqueNodes Map containing candidate unique nodes (output)
         */
        void findAnchors(std::map<NodeID, NodePos>& uniqueNodes) const;

        /**
         * Find the links to the anchors
         * @param signAnchors Links to anchors that are significant
         * @param nonSignAnchors Links to anchors that are not significant
         * @param dubiousAnchors Anchors that appear twice
         */
        void addLinksBetweenAnchors(std::map<NodePair, GaussVal>& signAnchors,
                         std::map<NodePair, GaussVal>& nonSignAnchors,
                         std::map<NodePair, GaussVal>& dubiousAnchors) const;

        void addObservations(Scaffold& scaffold);

        /**
         * Reverse all links (srcID => -srcID)
         */
        void reverse();

        /**
         * Print all paired node info to the screen
         */
        void printInfo() const;

        /**
         * Get the expected multiplicity, based on node collisions
         */
        double getExpSrcMultiplicity() const;
};

#endif
