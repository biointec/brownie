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

#include "leastsquares.h"
#include "nodenbh.h"
#include "scaffold.h"
#include "settings.h"

// ============================================================================
// OBSERVATIONS
// ============================================================================

void DBGraph::reverseObs(Observation& obs)
{
        NodeID srcID = obs.first.first;
        NodeID dstID = obs.first.second;

        obs.first = NodePair(dstID, srcID);
        obs.second = GaussVal(-obs.second.getAverage(), obs.second.getStd());
}

void DBGraph::complementObs(Observation& obs)
{
        NodeID srcID = obs.first.first;
        NodeID dstID = obs.first.second;

        double srcLen = DBGraph::graph->getSSNode(srcID).getMarginalLength();
        double dstLen = DBGraph::graph->getSSNode(dstID).getMarginalLength();

        obs.first = NodePair(-srcID, -dstID);
        obs.second = GaussVal(-obs.second.getAverage() + srcLen -
                              dstLen, obs.second.getStd());
}

void DBGraph::revComplObs(Observation& obs)
{
        reverseObs(obs);
        complementObs(obs);
}

void DBGraph::normalizeObs(Observation& obs)
{
        if (obs.first.first < 0)
                complementObs(obs);
}

// ============================================================================
// SCAFFOLD
// ============================================================================

void DBGraph::getLinksFromAnchor(NodeID srcID,
                                 const multimap<NodePair, GaussVal>& iLinks,
                                 multimap<NodePair, GaussVal>& oLinks)
{
        typedef multimap<NodePair, GaussVal>::const_iterator AnchorLinkIt;

        oLinks.clear();

        AnchorLinkIt it = iLinks.lower_bound(NodePair(abs(srcID), -numNodes));
        for ( ; it != iLinks.end(); it++) {
                if (it->first.first != abs(srcID))
                        break;

                Observation obs = *it;
                if (srcID < 0)
                        complementObs(obs);

                oLinks.insert(obs);
        }
}

void DBGraph::findLinksBetweenScaffolds(const Scaffold& scafA,
                                        const Scaffold& scafB,
                                        const multimap<NodePair, GaussVal>& iLinks,
                                        multimap<NodePair, GaussVal>& oLinks)
{
        typedef multimap<NodePair, GaussVal>::const_iterator AnchorLinkIt;

        oLinks.clear();

        set<NodeID> setA;
        for (CZeroMapIt it = scafA.zeroBegin(); it != scafA.zeroEnd(); it++)
                setA.insert(it->first);

        for (set<NodeID>::iterator it = setA.begin(); it != setA.end(); it++) {
                multimap<NodePair, GaussVal> tL;
                getLinksFromAnchor(*it, iLinks, tL);

                for (AnchorLinkIt e = tL.begin(); e != tL.end(); e++) {
                        NodeID dstID = e->first.second;

                        if (scafB.hasNode(dstID))
                                oLinks.insert(*e);
                }
        }
}

void DBGraph::filterAnchors()
{
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);
                if (!node.isAnchor())
                        continue;

                ObsMap links;
                getLinksFromAnchor(id, signLinks, links);

                for (CObsMapIt it1 = links.begin(); it1 != links.end(); it1++) {
                        CObsMapIt it2 = it1;

                        for (it2++; it2 != links.end(); it2++) {
                                NodeID nodeA = it1->first.second;
                                NodeID nodeB = it2->first.second;

                                if (nodeA == nodeB)
                                        continue;

                                const GaussVal &posA = it1->second;
                                const GaussVal &posB = it2->second;

                                if (!compareForOverlap(nodeA, posA, nodeB, posB)) {
                                        cout << "Src node: " << id << ", dst node " << nodeA << " in interval: " << posA << " - " << posA + getSSNode(nodeA).getMarginalLength()  << " overlaps with " << nodeB << ": " << posB << endl;
                                        node.setAnchor(false);
                                }
                        }
                }

                if (node.isAnchor())
                        continue;

                // if the node is no longer an anchor: erase the links
                for (CObsMapIt it1 = links.begin(); it1 != links.end(); it1++) {
                        cout << "Erasing: " << it1->first << endl;
                        signLinks.erase(it1->first);

                        NodePair rev = (it1->first.second > 0) ?
                                NodePair(it1->first.second, it1->first.first) :
                                NodePair(-it1->first.second, -it1->first.first);

                        cout << "Erasing: " << rev << endl;
                        signLinks.erase(rev);
                }
        }
}

void DBGraph::addLinksBetweenAnchors(NodeID srcID, const ReadLibrary& input,
                                     const vector<MappedRead>& mRead,
                                     map<NodePair, GaussVal>& sign,
                                     map<NodePair, GaussVal>& nonSign,
                                     map<NodePair, GaussVal>& dubious)
{
        SSNode src = getSSNode(srcID);
        double srcLen = src.getMarginalLength();

        GaussVal SB(srcLen, input.getFragmentStd());
        GaussVal SE(input.getPERKmerEnd(srcLen), input.getFragmentStd());

        NodeNBH nodeNBH(srcID, SB, SE);

        nodeNBH.createNBH(input, mRead, *this);

       /* if (srcID == 194) {
                nodeNBH.printInfo();
                exit(0);
        }*/
        nodeNBH.addLinksBetweenAnchors(sign, nonSign, dubious);
}

void DBGraph::findLinksBetweenAnchors()
{
    /*    for (size_t i = 0; i < settings.getInputCount(); i++) {
                const Input& input = settings.getInput(i);

                map<NodePair, GaussVal> thisSign;
                map<NodePair, GaussVal> thisNonSign;
                map<NodePair, GaussVal> thisDubious;

                for (NodeID id = -numNodes; id <= numNodes; id++) {
                        if (id == 0)
                                continue;

                        if (!getSSNode(id).isAnchor())
                                continue;

                        if (mReads[i].find(id) == mReads[i].end())
                                continue;

                        const vector<MappedRead>& mRead = mReads[i].find(id)->second;
                        addLinksBetweenAnchors(id, input, mRead, thisSign,
                                               thisNonSign, thisDubious);
                }

                makeLinksConsistent(thisSign);
                signLinks.insert(thisSign.begin(), thisSign.end());
                makeLinksConsistent(thisNonSign);
                nonSignLinks.insert(thisNonSign.begin(), thisNonSign.end());
        }*/
}

void DBGraph::buildPrimaryScaffold(set<Scaffold*>& scaffoldSet)
{
        typedef std::multimap<GaussVal, NodePair, SortGVByStd>::iterator SortedObsMapIt;

        // step 1: find the anchor nodes
        setAnchorNodes();

        // step 2: find the links between the anchors
        findLinksBetweenAnchors();

        // step 3: try to filter duplicated nodes
        filterAnchors();

        // step 3: mark the dubious nodes as non-significant
        /*cout << "Number of dubious links: " << dubiousLinks.size() << endl;
        for (ObsMapIt it = dubiousLinks.begin(); it != dubiousLinks.end(); it++)
                cout << "Dubious: " << it->first << ", " << it->second << endl;*/
       //         getSSNode(it->first.second).setAnchor(false);

        map<NodeID, Scaffold*> scaffoldMap;
        for (NodeID id = 1; id <= numNodes; id++) {
                if (!getSSNode(id).isAnchor())
                        continue;
                scaffoldMap[id] = new Scaffold(id);
        }

        // sort the links by link significance
        multimap<GaussVal, NodePair, SortGVByStd> sSignLinks;
        for (ObsMapIt it = signLinks.begin(); it != signLinks.end(); it++)
                sSignLinks.insert(pair<GaussVal, NodePair>(it->second, it->first));

      //  for (ObsMapIt it = signLinks.begin(); it != signLinks.end(); it++)
      //          cout << it->first << ": " << it->second << endl;

        cout << "Number of links between scaffolds: " << sSignLinks.size() / 2 << endl;
        for (SortedObsMapIt it = sSignLinks.begin(); it != sSignLinks.end(); it++) {
               // cout << it->second << ": " << it->first << endl;

                NodeID srcID = it->second.first;
                NodeID dstID = it->second.second;

                double distance = it->first.getAverage();
                double std = it->first.getStd();

                double actualLen;
                if (getDistance(srcID, dstID, abs(distance) + 20 * std, actualLen)) {
                        double dx = abs(distance - actualLen);
                        if (dx > 5.0 * std)
                               cout << "DISTANCE NOT OK from " << srcID << " to " << dstID << "; actual: " << actualLen << ", predicted: " << distance << " +/- " << std  << endl;
                } else {
                       cout << " NOT FOUND from " << it->second << endl;
                }
        }

        cout << " =================================== " << endl;

        size_t cnt = 0;
        bool verbose = false;

        for (SortedObsMapIt it = sSignLinks.begin(); it != sSignLinks.end(); it++) {

                //cout << "Handling link: " << cnt++ << " / " << sSignLinks.size() << "\r";
                //std::flush(cout);

                NodePair np = it->second;
                GaussVal gv = it->first;

                double lenA = getSSNode(np.first).getMarginalLength();
                double lenB = getSSNode(np.second).getMarginalLength();

                // find the first scaffold
                NodeID nodeA = np.first;
                Scaffold *scaffoldA = scaffoldMap[abs(nodeA)];
                if (scaffoldA->hasNode(-nodeA)) {
                        np = NodePair(-np.first, -np.second);
                        gv = GaussVal(-gv.getAverage() + lenA - lenB, gv.getStd());
                }

                // find the second scaffold
                NodeID nodeB = np.second;
                Scaffold *scaffoldB = scaffoldMap[abs(nodeB)];

                // this link was already added to the scaffold!
                if (scaffoldA == scaffoldB)
                        continue;

                // reverse complement the second scaffold, if necessary
                if (scaffoldB->hasNode(-nodeB))
                        scaffoldB->reverseComplement();

                multimap<NodePair, GaussVal> thisSignLinks;
                findLinksBetweenScaffolds(*scaffoldA, *scaffoldB, signLinks, thisSignLinks);
                multimap<NodePair, GaussVal> thisNonSignLinks;
                findLinksBetweenScaffolds(*scaffoldA, *scaffoldB, nonSignLinks, thisNonSignLinks);

         /*       if (!scaffoldA->canMerge(*scaffoldB, np, gv, thisSignLinks, thisNonSignLinks)) {
                        cout << "LINK: " << np << ", " << gv << endl;
                        cout << "A) Merging: " <<  endl;
                        cout << *scaffoldA << endl;
                        cout << "B) " << nodeB << endl;
                        cout << *scaffoldB << endl << endl;
                        for (multimap<NodePair, GaussVal>::iterator e = thisSignLinks.begin(); e != thisSignLinks.end(); e++) {
                                cout << "Otherlinks: " << e->first << ", " << e->second << endl;
                        }
                }*/

                scaffoldA->merge(*scaffoldB, np, gv, thisSignLinks);

                /*if (!validateScaffold(*scaffoldA)) {
                        cout << "LINK: " << it->first << ", " << it->second << endl;
                        cout << *scaffoldA << endl;
                        cout << " =========== " << endl << *scaffoldB << endl;
                        exit(0);
                }*/

                if (verbose) {
                        cout << "After merging: " << endl;
                        cout << *scaffoldA << endl;
                }

                for (CZeroMapIt n = scaffoldB->zeroBegin(); n != scaffoldB->zeroEnd(); n++)
                        scaffoldMap[abs(n->first)] = scaffoldA;

                delete scaffoldB;
        }

        cout << "Handling link: " << cnt++ << " / " << sSignLinks.size() << endl;

        scaffoldSet.clear();
        for (NodeID id = 1; id <= numNodes; id++) {
                if (!getSSNode(id).isAnchor())
                        continue;

                scaffoldSet.insert(scaffoldMap[id]);
        }
}

