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
#include <cmath>
#include <fstream>

using namespace std;

bool DBGraph::concatenation()
{
        bool simplifiedGraph = false;
        for (NodeID lID = -numNodes; lID <= numNodes; lID++) {
                if (lID == 0)
                        continue;

                SSNode left = getSSNode(lID);
                if (!left.isValid())
                        continue;

                if (left.getNumRightArcs() != 1)
                        continue;

                NodeID rID = left.rightBegin()->getNodeID();
                SSNode right = getSSNode(rID);

                // don't merge palindromic repeats
                if (rID == -lID)
                        continue;

                if (right.getNumLeftArcs() != 1)
                        continue;

                // at this point, we can merge left en right
                table->mergeLeftToRight(-rID, -lID, right.getMarginalLength(),
                                        left.getMarginalLength(), false);

                left.deleteRightArc(rID);
                right.deleteLeftArc(lID);

                left.inheritRightArcs(right);

                right.invalidate();

                deque<SSNode> deq;
                deq.push_back(left);
                deq.push_back(right);

                string str;
                convertNodesToString(deq, str);

                left.setSequence(str);
                simplifiedGraph = true;
                lID--;
        }

        return simplifiedGraph;
}

bool DBGraph::checkReduction(NodeID lID, vector<NodeID> vmID, NodeID rID)
{
        if (getSSNode(lID).getRightArc(vmID.front()) == NULL)
                cerr << "Trying a reduction of a non-existing path" << endl;

        // find all occurences of the complete path
        deque<SSNode> deq;
        deq.push_back(getSSNode(lID));
        for (size_t i = 0; i < vmID.size(); i++)
                deq.push_back(getSSNode(vmID[i]));
        deq.push_back(getSSNode(rID));

        string str;
        convertNodesToString(deq, str);

        vector<vector<size_t> > pos, posRC;
        findAllTrueOccurences(str, pos, posRC);

        // find all occurences of the lID and vmID.front()
        deq.clear();
        deq.push_back(getSSNode(lID));
        deq.push_back(getSSNode(vmID.front()));
        convertNodesToString(deq, str);

        vector<vector<size_t> > fpos, fposRC;
        findAllTrueOccurences(str, fpos, fposRC);

        // find all occurences of the vmID.back() and rID
        deq.clear();
        deq.push_back(getSSNode(vmID.back()));
        deq.push_back(getSSNode(rID));
        convertNodesToString(deq, str);

        vector<vector<size_t> > lpos, lposRC;
        findAllTrueOccurences(str, lpos, lposRC);

        if ((pos.size() + posRC.size()) != (fpos.size() + fposRC.size())) {
                cerr << "Error while deleting arc from " << lID << " to " << vmID.front() << endl;
                return false;
        }
        if ((pos.size() + posRC.size()) != (lpos.size() + lposRC.size())) {
                cerr << "Error while deleting arc from " << vmID.back() << " to " << rID << endl;
                return false;
        }

        return true;
}

bool DBGraph::reduction(NodeID lID, vector<NodeID> vmID, NodeID rID)
{
        // check if the reduction is a palindrome
        bool palindrome = ((vmID.size() == 2) && (vmID[0] == -vmID[1]));

        // make sure all nodes in the path are unique
        vector<NodeID> path;
        path.reserve(vmID.size() + 2);
        path.push_back(lID);
        path.insert(path.begin() + 1, vmID.begin(), vmID.end());
        path.push_back(rID);
        for (size_t i = 0; i < path.size(); i++)
                path[i] = abs(path[i]);
        sort(path.begin(), path.end());
        if (set<NodeID>(path.begin(), path.end()).size() != path.size()) {
                cout << "DUPLICATE FOUND " << palindrome << endl;
                if (!palindrome)
                        return false;
        }

        /*cout << "Reduction " << lID << " ";
        for (size_t i = 0; i < vmID.size(); i++)
                cout << vmID[i] << " ";
        cout << rID << endl;*/

        checkReduction(lID, vmID, rID);

        SSNode lNode = getSSNode(lID);
        SSNode rNode = getSSNode(rID);
        vector<SSNode> mNode;
        for (size_t i = 0; i < vmID.size(); i++)
                mNode.push_back(getSSNode(vmID[i]));

        size_t lSize = lNode.getMarginalLength();
        size_t rSize = rNode.getMarginalLength();

        bool leftMerge = (lNode.getNumRightArcs() == 1);
        bool rightMerge = (rNode.getNumLeftArcs() == 1);

        vector<bool> duplication(vmID.size(), true);
        if ((vmID.size() == 1) &&
            (mNode[0].getNumLeftArcs() <= 1) &&
            (mNode[0].getNumRightArcs() <= 1))
                duplication[0] = false;

        if (palindrome)
                duplication[1] = false;

        if (leftMerge && rightMerge) {

                //   cout << "Attempting to merge: " << lID << " "
                //        << mID << " " << rID << endl;

                // update kmerNode table (tracking of kmers)
                size_t newSize = lSize;
                for (size_t i = 0; i < vmID.size(); i++) {
                        NodeID mID = vmID[i];
                        SSNode mNode = getSSNode(mID);
                        size_t mSize = mNode.getMarginalLength();

                        table->mergeLeftToRight(-mID, -lID, mSize,
                                                newSize, duplication[i]);

                        newSize += mSize;
                }

                table->mergeLeftToRight(-rID, -lID, rSize, newSize, false);

                mNode.front().deleteLeftArc(lID);
                mNode.back().deleteRightArc(rID);
                lNode.deleteRightArc(vmID.front());
                rNode.deleteLeftArc(vmID.back());

                lNode.inheritRightArcs(rNode);
                rNode.invalidate();

                deque<SSNode> deq;
                deq.push_back(lNode);
                for (size_t i = 0; i < vmID.size(); i++)
                        deq.push_back(mNode[i]);
                deq.push_back(rNode);

                string str;
                convertNodesToString(deq, str);

                lNode.setSequence(str);

                if (palindrome) {
                        mNode.front().deleteRightArc(vmID.back());
                        mNode.front().invalidate();
                }
                return true;
        } else if (leftMerge) {

                //   cout << "Attempting to LEFT merge: " << lID
                //        << " " << mID << " " << rID << endl;

                // update kmerNode table (tracking of kmers)
                size_t newSize = lSize;
                for (size_t i = 0; i < vmID.size(); i++) {
                        NodeID mID = vmID[i];
                        SSNode mNode = getSSNode(mID);
                        size_t mSize = mNode.getMarginalLength();

                        table->mergeLeftToRight(-mID, -lID, mSize,
                                                newSize, duplication[i]);

                        newSize += mSize;
                }

                lNode.replaceRightArc(vmID.front(), rID);
                rNode.replaceLeftArc(vmID.back(), lID);
                mNode.front().deleteLeftArc(lID);
                mNode.back().deleteRightArc(rID);

                deque<SSNode> deq;
                deq.push_back(lNode);
                for (size_t i = 0; i < vmID.size(); i++)
                        deq.push_back(mNode[i]);

                string str;
                convertNodesToString(deq, str);

                lNode.setSequence(str);
                if (palindrome) {
                        mNode.front().deleteRightArc(vmID.back());
                        mNode.front().invalidate();
                }
                return true;

        } else if (rightMerge) {

                //   cout << "Attempting to RIGHT merge: " << lID << " "
                //        << mID << " " << rID << endl;

                // update kmerNode table (tracking of kmers)
                size_t newSize = rSize;
                for (ssize_t i = vmID.size() - 1; i >= 0; i--) {
                        NodeID mID = vmID[i];
                        SSNode mNode = getSSNode(mID);
                        size_t mSize = mNode.getMarginalLength();

                        table->mergeLeftToRight(mID, rID, mSize,
                                                newSize, duplication[i]);

                        newSize += mSize;
                }

                lNode.replaceRightArc(vmID.front(), rID);
                rNode.replaceLeftArc(vmID.back(), lID);
                mNode.front().deleteLeftArc(lID);
                mNode.back().deleteRightArc(rID);

                deque<SSNode> deq;
                for (size_t i = 0; i < vmID.size(); i++)
                        deq.push_back(mNode[i]);
                deq.push_back(rNode);

                string str;
                convertNodesToString(deq, str);

                rNode.setSequence(str);
                if (palindrome) {
                        mNode.front().deleteRightArc(vmID.back());
                        mNode.front().invalidate();
                }

                return true;
        }

        return false;
}

bool DBGraph::covBasedReduction(NodeID lID, vector<NodeID> vmID, NodeID rID)
{
        // make sure the initial and last node are unique
        SSNode lNode = getSSNode(lID);
        SSNode rNode = getSSNode(rID);

      //  if (lNode.getMultiplicity() != 1)
      //          cerr << "Doing reduction with left node multiplicity != 1" << endl;
      //  if (rNode.getMultiplicity() != 1)
      //          cerr << "Doing reduction with right node multiplicity != 1" << endl;

        // in principle, unique nodes only have a single arc
        bool leftMerge = (lNode.getNumRightArcs() == 1);
        bool rightMerge = (rNode.getNumLeftArcs() == 1);

        vector<SSNode> mNode;
        for (size_t i = 0; i < vmID.size(); i++)
                mNode.push_back(getSSNode(vmID[i]));

        size_t lSize = lNode.getMarginalLength();
        size_t rSize = rNode.getMarginalLength();

        for (size_t i = 0; i < vmID.size(); i++) {
                if (abs(vmID[i]) == abs(lID))
                        cerr << "Left node also occurs in middle nodes" << endl;
                if (abs(vmID[i]) == abs(rID))
                        cerr << "Right node also occurs in middle nodes" << endl;
        }

        // count the number of times a node is visited
        map<NodeID, size_t> nodeOcc;
        for (size_t i = 0; i < vmID.size(); i++) {
                map<NodeID, size_t>::iterator it = nodeOcc.find(abs(vmID[i]));
                if (it == nodeOcc.end())
                        nodeOcc[abs(vmID[i])] = 1;
                else
                        it->second++;
        }

        vector<bool> duplication(vmID.size(), true);
        for (int i = vmID.size() - 1; i >= 0; i--) {
                if (nodeOcc[abs(vmID[i])] == getSSNode(vmID[i]).getRoundMult()) {
                        duplication[i] = false;
                        nodeOcc[abs(vmID[i])] = getSSNode(vmID[i]).getRoundMult() + 1;
                }
        }

        if (leftMerge && rightMerge) {

                //   cout << "Attempting to merge: " << lID << " "
                //        << mID << " " << rID << endl;

                // update kmerNode table (tracking of kmers)
                size_t newSize = lSize;
                for (size_t i = 0; i < vmID.size(); i++) {
                        NodeID mID = vmID[i];
                        SSNode mNode = getSSNode(mID);
                        size_t mSize = mNode.getMarginalLength();

                        table->mergeLeftToRight(-mID, -lID, mSize,
                                                newSize, duplication[i]);

                        newSize += mSize;
                }

                table->mergeLeftToRight(-rID, -lID, rSize, newSize, false);


                mNode.front().deleteLeftArc(lID);
                mNode.back().deleteRightArc(rID);
                lNode.deleteRightArc(vmID.front());
                rNode.deleteLeftArc(vmID.back());

                lNode.inheritRightArcs(rNode);
                rNode.invalidate();

                // invalidate middle nodes
                for (size_t i = 0; i < vmID.size(); i++) {
                        if (duplication[i])
                                continue;
                        mNode[i].invalidate();
                }

                for (size_t i = 0; i < vmID.size(); i++) {
                        if (mNode[i].isValid())
                                continue;
                        if (i > 0) {
                                mNode[i].deleteLeftArc(vmID[i-1]);
                                mNode[i-1].deleteRightArc(vmID[i]);
                        }

                        if (i < vmID.size()-1) {
                                mNode[i].deleteRightArc(vmID[i+1]);
                                mNode[i+1].deleteLeftArc(vmID[i]);
                        }

                        //mNode[i].setMultiplicity(mNode[i].getMultiplicity() - 1);
                }

                deque<SSNode> deq;
                deq.push_back(lNode);
                for (size_t i = 0; i < vmID.size(); i++)
                        deq.push_back(mNode[i]);
                deq.push_back(rNode);

                string str;
                convertNodesToString(deq, str);

                lNode.setSequence(str);

                return true;
        }

        return false;
}
