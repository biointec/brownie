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
#include "settings.h"


using namespace std;

bool DBGraph::clipTips(int round)
{
        cout <<endl<< " =================== Removing tips ===================" << endl;

#ifdef DEBUG
        size_t tp=0, tn=0, fp=0,fn=0;
        size_t tps=0, tns=0, fps=0,fns=0;
        size_t tpj=0, tnj=0, fpj=0,fnj=0;
#endif

        size_t numDeleted = 0, numTotal = 0;
        cout << "Cut-off value for removing tips is: " << redLineValueCov << endl;
        double threshold = redLineValueCov;
        size_t numThreads= settings.getNumThreads();
        #pragma omp parallel num_threads( numThreads)
        {
                size_t numDeletedLocal = 0, numTotalLocal = 0;
                #pragma omp for
                for (NodeID id = 1; id <= numNodes; id++) {
                        SSNode node = getSSNode(id);
                        if (!node.isValid())
                                continue;
                        numTotalLocal++;

                        // check for dead ends
                        bool leftDE = (node.getNumLeftArcs() == 0);
                        bool rightDE = (node.getNumRightArcs() == 0);
                        if (!leftDE && !rightDE)
                                continue;

                        SSNode startNode = (rightDE) ? getSSNode(-id) : getSSNode(id);
                        bool isolated = rightDE && leftDE;
                        bool joinedTip = startNode.getNumRightArcs() > 1;
                        if (isolated||joinedTip)
                                threshold=this->safeValueCov;
                        bool remove = false;
                        if ((startNode.getNodeKmerCov() <threshold) &&
                                (startNode.getMarginalLength() < maxNodeSizeToDel))
                                #pragma omp critical
                                remove = removeNode(startNode);

                        if (remove)
                                numDeletedLocal++;

                        #ifdef DEBUG
                        updateStaInClipTip(id , remove, isolated, joinedTip, fps,fpj, fp, tps, tpj,tp, tns, tnj,tn, fns, fnj, fn );
                        #endif
                }
                #pragma omp critical
                {
                        numDeleted =+numDeletedLocal ;
                        numTotal = +numTotalLocal;

                }
        }

        cout << "Clipped " << numDeleted << "/" << numTotal << " nodes" << endl;
        #ifdef DEBUG
        printStatisticInClipTip(fps,fpj, fp, tps, tpj,tp, tns, tnj,tn, fns, fnj, fn );
        #endif
        return  numDeleted > 0;
}
void DBGraph::updateStaInClipTip(NodeID id , bool remove, bool isolated, bool joinedTip, size_t& fps
        , size_t &fpj, size_t &fp, size_t& tps, size_t &tpj,size_t &tp, size_t &tns, size_t &tnj,
        size_t &tn, size_t &fns, size_t& fnj,size_t& fn ){
        #ifdef DEBUG
        if (remove) {
                if (trueMult.size()>0&& trueMult[id] > 0) {
                        if (isolated)
                                #pragma omp critical
                                fps++;
                        else if(joinedTip)
                                #pragma omp critical
                                fpj++;
                        else
                                #pragma omp critical
                                fp++;
                } else {
                        if(isolated)
                                #pragma omp critical
                                tps++;
                        else if(joinedTip)
                                #pragma omp critical
                                tpj++;
                        else
                                #pragma omp critical
                                tp++;
                }
        } else {
                if (trueMult.size()>0&& trueMult[id] > 0) {
                        if(isolated)
                                #pragma omp critical
                                tns++;
                        else if(joinedTip)
                                #pragma omp critical
                                tnj++;
                        else
                                #pragma omp critical
                                tn++;
                } else {
                        if(isolated)
                                #pragma omp critical
                                fns++;
                        else if(joinedTip)
                                #pragma omp critical
                                fnj++;
                        else
                                #pragma omp critical
                                fn++;
                }
        }
        #endif
}

void DBGraph::printStatisticInClipTip(size_t& fps
, size_t &fpj, size_t &fp, size_t& tps, size_t &tpj,size_t &tp, size_t &tns, size_t &tnj,
size_t &tn, size_t &fns, size_t& fnj,size_t& fn){
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
}
