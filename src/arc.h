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

#ifndef ARC_H
#define ARC_H

#include "global.h"
#include <atomic>

// ============================================================================
// ARC CLASS
// ============================================================================

class Arc {

private:
        NodeID nodeID;                  // ID of node to which the arc points
        std::atomic<Coverage> cov;      // arc coverage

public:
        /**
         * Default constructor
         */
        Arc() : nodeID(0), cov(0) {};

        /**
         * Copy constructor
         */
        Arc(const Arc& rhs) : nodeID(rhs.nodeID), cov(rhs.cov.load()) {}

        /**
         * Assignment operator
         */
        Arc& operator=(const Arc& rhs) {
                nodeID = rhs.nodeID;
                cov = rhs.cov.load();
                return *this;
        }

        /**
         * Set the target nodeID the arc is pointing to
         * @param targetNodeID The ID of the node the arc is pointing to
         */
        void setNodeID(NodeID targetNodeID) {
                nodeID = targetNodeID;
        }

        /**
         * Get the target nodeID and the side of the target node
         * @return targetNodeID The ID of the node the arc is pointing to
         */
        NodeID getNodeID() const {
                return nodeID;
        }

        /**
         * Set the coverage of the arc
         * @param coverage Coverage of the arc
         */
        void setCoverage(Coverage coverage) {
                cov = coverage;
        }

        /**
         * Get the arc coverage
         * @return The arc coverage
         */
        Coverage getCoverage() const {
                return cov;
        }

        /**
         * Increment of the coverage
         */
        void incReadCov() {
                cov++;
        }

        /**
         * Delete arc (mark as invalid)
         */
        void deleteArc() {
                nodeID = 0;
        }

        /**
         * Check if arc is valid (== not deleted)
         * @return True of false
         */
        bool isValid() const {
                return nodeID != 0;
        }

        /**
         * Write a node to file
         * @param ofs Open output file stream
         */
        void write(std::ofstream& ofs) const {
                ofs.write((char*)&nodeID, sizeof(nodeID));
                ofs.write((char*)&cov, sizeof(cov));
        }

        /**
         * Load a node from file
         * @param ifs Open input file stream
         */
        void read(std::ifstream& ifs) {
                ifs.read((char*)&nodeID, sizeof(nodeID));
                ifs.read((char*)&cov, sizeof(cov));
        }
};

#endif
