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

#include "graph.h"
#include "kmernode.h"
#include "readfile/fastafile.h"
#include "settings.h"
#include <cmath>
#include <fstream>

using namespace std;



bool DBGraph::consecutiveNPP(const NodePosPair& npp1,
                              const NodePosPair& npp2, size_t offset) const
{
        // if one of the nnp's is not found, the answer is no
        if ((npp1.getNodeID() == 0) || (npp2.getNodeID() == 0))
                return false;

        if (npp1.getNodeID() == npp2.getNodeID()) {
                // if they map to the same node, check if they're consecutive
                return (npp1.getPosition() == (npp2.getPosition()-offset));
        } else {
                // if they don't map to the same node,
                // check if the nodes are connected properly
                SSNode node = getSSNode(npp1.getNodeID());
                if ((node.getMarginalLength() - offset) != npp1.getPosition())
                        return false;

                if (npp2.getPosition() != 0)
                        return false;

                if (node.getRightArc(npp2.getNodeID()) == NULL)
                        return false;
        }

        return true;
}

bool DBGraph::connectedNodes(NodeID leftID, NodeID rightID) const
{
        SSNode left = getSSNode(leftID);
        return (left.getRightArc(rightID) != NULL);
}

struct sortBySize
{
        bool operator()(const int& left, const int& right) const {
                return left > right;
        }
};


// ============================================================================
// PERM.CPP PUBLIC
// ============================================================================

void DBGraph::perm()
{
      /*  table = new KmerNodeTable(settings, numNodes);

#ifndef NDEBUG
        readReferenceGenome();
#endif

        // =========== NEW IDEA ========

        // read the paired-end reads in memory
        nodeDist.clear();
        mReads.resize(settings.getInputCount());
        for (size_t i = 0; i < settings.getInputCount(); i++) {
                const Input &input = settings.getInput(i);

                cout << "Mapping paired-end file " << i+1 << "/"
                     << settings.getInputCount() << ": "
                     << input.getFilename()
                     << " (type: " << input.getFileType()
                     << ", reads: " << input.getReadType()
                     << ", read direction: " << input.getReadDirType() << ")"
                     << endl;

                cout << "Insert lenght: " << input.getFragmentSize() << " +/- "
                     << input.getFragmentStd() << endl;

                readPairedReads(input, mReads[i]);
        }

        dbReductions();*/
}
