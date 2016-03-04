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

#ifndef CORRECTGRAPH_H
#define CORRECTGRAPH_H

#include "global.h"

// ============================================================================
// DIJKSTRA AUXILIARY CLASSES
// ============================================================================

class PathDFS {
public:
        NodeID nodeID;
        size_t length;

        /**
         * Default constructor
         */
        PathDFS(NodeID nodeID, size_t length) :
                nodeID(nodeID), length(length) {};
};

struct PathDFSComp {

        /**
         * Compare two paths (by length)
         */
        bool operator()(const PathDFS& f, const PathDFS& s) {
                return f.length >= s.length;
        }
};

#endif
