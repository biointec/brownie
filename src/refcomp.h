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
         * @param refID Reference identifier
         * @param begin Start position of the breakpoint
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
// REFERENCE COMPARISON CLASS
// ============================================================================

class RefComp {

private:
        std::vector<std::string> reference;     // actual reference genome

public:
        /**
         * Default constructor
         * @param refFilename File name of the reference genome
         */
        RefComp(const std::string& refFilename);

        /**
         * Validate a de Bruijn graph
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
         * Calculate the true node multiplicity
         * @param dbg A const-ref to the de Bruijn graph
         * @param multiplicity Multiplicity vector
         */
        void getNodeMultiplicity(const DBGraph& dbg,
                                 std::vector<size_t>& multiplicity);
};

#endif
