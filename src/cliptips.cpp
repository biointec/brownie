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


using namespace std;

bool DBGraph::clipTips(int round)
{
        cout << " ===== Removing tips =====" << endl;

#ifdef DEBUG
        size_t tp=0, tn=0, fp=0,fn=0;
        size_t tps=0, tns=0, fps=0,fns=0;
        size_t tpj=0, tnj=0, fpj=0,fnj=0;
#endif

        size_t numDeleted = 0, numTotal = 0;
        cout << "Cut-off value for removing tips is: " << redLineValueCov << endl;
        double threshold = redLineValueCov;

        for (NodeID id = 1; id <= numNodes; id++) {

                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;
                numTotal++;

                // check for dead ends
                bool leftDE = (node.getNumLeftArcs() == 0);
                bool rightDE = (node.getNumRightArcs() == 0);
                if (!leftDE && !rightDE)
                        continue;

                SSNode startNode = (rightDE) ? getSSNode(-id) : getSSNode(id);
                bool isolated = rightDE && leftDE;
                bool joinedTip = startNode.getNumRightArcs() > 1;

                bool remove = false;
                if ((startNode.getNodeKmerCov() <threshold) &&
                    (startNode.getMarginalLength() < maxNodeSizeToDel))
                        remove = removeNode(startNode);

                if (remove)
                        numDeleted++;

#ifdef DEBUG

                if (remove) {
                        if (trueMult[id] > 0) {
                                if (isolated)
                                        fps++;
                                else if(joinedTip)
                                        fpj++;
                                else
                                        fp++;
                        } else {
                                if(isolated)
                                        tps++;
                                else if(joinedTip)
                                        tpj++;
                                else
                                        tp++;
                        }
                } else {
                        if (trueMult[id] > 0) {
                                if(isolated)
                                        tns++;
                                else if(joinedTip)
                                        tnj++;
                                else
                                        tn++;
                        } else {
                                if(isolated)
                                        fns++;
                                else if(joinedTip)
                                        fnj++;
                                else
                                        fn++;
                        }
                }
#endif
        }

        cout << "Clipped " << numDeleted << "/" << numTotal << " nodes" << endl;
#ifdef DEBUG
        cout << "****************************************" << endl;
        cout << "Isolated TP: " << tps << "\tTN: "<< tns << "\tFP: " << fps << "\tFN: "<< fns << endl;
        cout << "Sensitivity: " << 100.0 * Util::getSensitivity(tps, fns) << "%" << endl;
        cout << "Specificity: " << 100.0 * Util::getSpecificity(tns, fps) << "%" << endl;
        cout << "****************************************" << endl;
        cout << "Joined TP: " << tpj << "\tTN: " << tnj << "\tFP: " << fpj << "\tFN: " << fnj << endl;
        cout << "Sensitivity: " << 100.0 * Util::getSensitivity(tpj, fnj) << "%" << endl;
        cout << "Specificity: " << 100.0 * Util::getSpecificity(tnj, fpj) << "%" << endl;
        cout << "****************************************" << endl;
        cout << "Tip TP: " << tp << "\tTN: " << tn << "\tFP: " << fp << "\tFN: " << fn << endl;
        cout << "Sensitivity: " << 100.0 * Util::getSensitivity(tp, fn) << "%" << endl;
        cout << "Specificity: " << 100.0 * Util::getSpecificity(tn, fp) << "%" << endl;
#endif

        return  numDeleted > 0;
}

