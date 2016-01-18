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

#ifndef COMPONENT_H
#define COMPONENT_H

#include <map>

class Component{
private:
        std::map<int, size_t> N;
public:
        size_t Size;
        size_t numOfNodes;
        size_t numOfArcs;
        double nodeKmerCov;
        size_t largestNodeSize;
        Component() :  N(), Size(0),numOfNodes(0),numOfArcs(0),nodeKmerCov(0), largestNodeSize(0)  {}
        std::map<int, size_t> get_N() const {
                return N;
        }
        void set_N(int x, size_t Nx) {
                N[x] = Nx;
        }
        size_t get_N(int x) const {
                return N.at(x);
        }
};

#endif
