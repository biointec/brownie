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

void DBGraph::removeTip ( SSNode &lNode, SSNode &rNode )
{
        // first, disconnect the right node from all its right neighbors
        for ( ArcIt it = rNode.rightBegin(); it != rNode.rightEnd(); it++ ) {
                SSNode rrNode = getSSNode ( it->getNodeID() );
                bool result = rrNode.deleteLeftArc ( rNode.getNodeID() );
                assert ( result );
        }
        rNode.deleteAllRightArcs();

        // the left node should not be left connected (it's a tip)
        assert ( lNode.getNumLeftArcs() == 0 );

        // delete all nodes
        SSNode curr = lNode;
        while ( curr != rNode ) {
                SSNode next = getSSNode ( curr.rightBegin()->getNodeID() );

                //comment by mahdi

                curr.deleteAllLeftArcs();
                curr.deleteAllRightArcs();
                curr.invalidate();

                if ( abs(curr.getNodeID()) ) {
                        // cout << "\tClipped a tip that should exist, namely node " << curr.getNodeID() << endl;
                }

                // if the statement below is true, we have an isolated hairpin
                // we can simply stop here, because the nodes and arcs on the
                // return path have already been destroyed
                if ( next.getNodeID() == -curr.getNodeID() ) {
                        break;
                }

                curr = next;
        }
        if ( trueMult[abs(rNode.getNodeID())] ) {
                // cout << "\tClipped a tip that should exist, namely node " << rNode.getNodeID() << endl;
        }
        rNode.deleteAllLeftArcs();
        rNode.invalidate();
}

bool DBGraph::clipTipFromNode ( SSNode &startNode )
{
        // size_t cutOffLength = Kmer::getWordSize() * 5;
        size_t totalLength = 0;
        //DBGraph graph ( settings );
        SSNode currNode = startNode, prevNode = startNode;
        while ( currNode.getNumLeftArcs() < 2 &&
                currNode.getNumRightArcs() < 2 ) {
                // update the total length, taking the overlap into account
                totalLength += currNode.getLength() - Kmer::getK() + 1;
        //if(currNode.getExpMult()/currNode.getMarginalLength()>10)
        //  return false;

        //          if (totalLength >= estimatedKmerCoverage)
        //return false;

        //if (currNode.getRoundMult() > 0)
        //        return false;

        // isolated snipset

        //comment by mahdi
        //  if (currNode.getExpMult() / currNode.getMarginalLength() >  this->estimatedKmerCoverage*.2 )
        //    return false;
        double  probalility=  gsl_cdf_gaussian_P((currNode.getNodeKmerCov()-estimatedKmerCoverage)/estimatedMKmerCoverageSTD,1);
        if (probalility>.2)
                return false;

        //cout << "\tClipped a tip that should exist, namely node " << curr.getNodeID() << endl;
        if ( currNode.getNumRightArcs() == 0 ) {

                removeTip ( startNode, currNode );
                return true;


        }

        // move to the next node
        prevNode = currNode;
        currNode = getSSNode ( currNode.rightBegin()->getNodeID() );

                }

                // joined tips
                if ( currNode.getNumLeftArcs() < 2 ) {
                        return false;
                }





                //    if (prevNode.getExpMult() / prevNode.getMarginalLength() >  this->estimatedKmerCoverage*.2 )
                //         return false;
                double  probalility=  gsl_cdf_gaussian_P((prevNode.getNodeKmerCov()-estimatedKmerCoverage)/estimatedMKmerCoverageSTD,1);
                if (probalility>.2)
                        return false;

                removeTip ( startNode, prevNode );
                return true;


}

bool DBGraph::clipTipFromNodeHard ( SSNode &startNode )
{
        // assert it is a tip
        assert ( startNode.getNumLeftArcs() == 0 );

        size_t totalLength = Kmer::getK() - 1;
        SSNode currNode = startNode, prevNode = startNode;

        while ( currNode.getNumLeftArcs() < 2 &&
                currNode.getNumRightArcs() < 2 ) {
                totalLength += currNode.getMarginalLength();

        // isolated snipset
        if ( currNode.getNumRightArcs() == 0 ) {
                break;
        }

        prevNode = currNode;
        currNode = getSSNode ( currNode.rightBegin()->getNodeID() );
                }

                // special case: WARNING is this correct ?!
                if ( currNode == startNode ) {
                        totalLength = startNode.getLength();
                }

                if ( totalLength >= 2*Kmer::getK() ) {
                        return false;
                }

                removeTip ( startNode, prevNode );
                return true;
}

bool DBGraph::clipTips ( bool hard )
{
        bool graphModified = false;

        // choose the correct function pointer
        bool ( DBGraph::*clipTipFromNodePtr ) ( SSNode & ) = ( hard ) ?
        &DBGraph::clipTipFromNodeHard :
        &DBGraph::clipTipFromNode;

        for ( NodeID id = 1; id <= numNodes; id++ ) {
                if (id==0)
                        continue;

                DSNode &node = getDSNode ( id );
                if ( !node.isValid() ) {
                        continue;
                }


                // check for dead ends
                bool leftDE = ( node.getNumLeftArcs() == 0 );
                bool rightDE = ( node.getNumRightArcs() == 0 );
                if ( !leftDE && !rightDE ) {
                        continue;
                }
                SSNode startNode = ( rightDE ) ? getSSNode ( -id ) : getSSNode ( id );



                if(startNode.getNodeID()==1) {

                        int stop=1;
                        stop++;
                }
                if ( ( *this.*clipTipFromNodePtr ) ( startNode ) ) {
                        graphModified = true;
                }
        }

        return graphModified;
}

bool DBGraph::clipTips(int round)
{
        cout<<"*********************<<clip Tips starts>>......................................... "<<endl;
        int numInitial = 0;
        double tp=0, tn=0, fp=0,fn=0;
        cout<<"cut off value for removing tips is: "<<this->redLineValueCov<<endl;

        for ( NodeID id = 1; id <= numNodes; id++ ) {


                SSNode node = getSSNode ( id );
                if ( !node.isValid() ) {
                        continue;
                }
                // check for dead ends
                bool leftDE = ( node.getNumLeftArcs() == 0 );
                bool rightDE = ( node.getNumRightArcs() == 0 );

                if ( !leftDE && !rightDE ) {
                        continue;
                }
                SSNode startNode = ( rightDE ) ? getSSNode ( -id ) : getSSNode ( id );
                bool singleNode=rightDE&&leftDE;
                double threshold=singleNode?this->certainVlueCov:this->redLineValueCov;
                if (startNode.getNodeKmerCov()<threshold )
                {
                        if (removeNode(startNode)){
                                #ifdef DEBUG
                                if (trueMult[abs( startNode.getNodeID())]>0) {
                                        fp++;
                                } else {
                                        tp++;
                                }
                                #endif
                        }
                        else{
                                #ifdef DEBUG
                                if (trueMult[abs( startNode.getNodeID())]>0) {
                                        tn++;
                                } else {
                                        fn++;
                                }
                                #endif
                        }
                }else{
                        #ifdef DEBUG
                        if (trueMult[abs( startNode.getNodeID())]>0) {
                                tn++;
                        } else {
                                fn++;
                        }
                        #endif
                }
        }



        cout << "Clipped " << tp+fp << "/" << numInitial << " nodes." << endl;
        #ifdef DEBUG
        cout<< "TP:	"<<tp<<"	TN:	"<<tn<<"	FP:	"<<fp<<"	FN:	"<<fn<<endl;
        cout << "Sensitivity: ("<<100*((double)tp/(double)(tp+fn))<<"%)"<<endl;
        cout<<"Specificity: ("<<100*((double)tn/(double)(tn+fp))<<"%)"<<endl;
        #endif
        return  tp+fp>0;

}

