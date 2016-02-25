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

#include "graph.h"
#include "settings.h"
#include "library.h"

#include <queue>

using namespace std;

void DBGraph::extractPath(NodeID currID, const vector<NodeID>& prevNode) const
{
        vector<NodeID> path;
        path.push_back(currID);

        while (true) {
                currID = prevNode[currID + numNodes];
                if (currID == 0)
                        break;
                path.push_back(currID);
        }

        reverse(path.begin(), path.end());
        for (auto it : path)
                cout << it << " ";
}

bool DBGraph::hasLowCovNode(SSNode root){

        ArcIt it = root.rightBegin();
        for (int i=0;i<root.getNumRightArcs();i++){
                SSNode node =getSSNode(it->getNodeID());
                if (node.getAvgKmerCov()<this->cutOffvalue)
                        return true;
                it++;
        }
        return false;

}



