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

#include "settings.h"
#include "graph.h"



class Component{
private:

        std::set<NodeID> NodesID;


public:
        size_t componentID;
        size_t componentSize;
        size_t numOfNodes;
        double componentKmerCov;

        Component (){

        }

        Component(set<NodeID> &IDs, size_t ID, size_t size ,double cov) : NodesID(IDs),componentID(ID), componentSize(size),numOfNodes(IDs.size()) ,componentKmerCov(cov){
        }
        /**
         * keep the nodes ID in component in a set of nodes Id
         * @param IDS nodesID in a set data structure
         */
        void setNodeIDSet(set<NodeID> IDs){
                NodesID=IDs;
        }
        /**
         * returns the set contains Nodes ID in the comonent
         *
         */
        set<NodeID> getNodeIDSet(){
                return NodesID;
        }
        int getComponentID(){
                return componentID;
        }

};


class ComponentHandler{
private:

        DBGraph &dbg;
        Settings &settings;
        std::map<NodeID, size_t > nodeComponentMap;

public:
        std::map< size_t, Component> componentsMap;
        vector<Component> components;
        size_t numberOfcomponents;

        ComponentHandler( DBGraph& g,Settings& s):dbg(g),settings(s),numberOfcomponents(0){}

        
        /**
         * Add new coponent to the list of comonent
         * @param nodesIdSet set of nodes which are in the component
         * @return make a new component and returns it
         *
         */
        size_t getComponentID(NodeID nodeID){
                auto it = nodeComponentMap.find(nodeID);
                if (it !=nodeComponentMap.end())
                        return it->second;
                else{
                        it = nodeComponentMap.find(-nodeID);
                        if (it !=nodeComponentMap.end())
                                return it->second;
                        else
                                return 0;
                }
        }
        void addComponent(set<NodeID>&  nodesIdSet){
                numberOfcomponents++;
                size_t size =0;
                double cov = 0;
                for (auto it:nodesIdSet){
                        size = size + dbg.getSSNode(it).getMarginalLength();
                        cov = cov + dbg.getSSNode(it).getKmerCov();
                }
                cov = cov / (double) size;
                Component newComonent(nodesIdSet, numberOfcomponents,size, cov);
                components.push_back(newComonent);
                for (auto it:nodesIdSet){
                        nodeComponentMap[it] = newComonent.componentID ;

                }
                componentsMap[newComonent.getComponentID()] = newComonent;
        }
        void reportStatistics ();
        /**
         * find number of disjoint components in the graph
         * @return number of comonents in the graph
         *
         */
        void findComponentsInGraph( size_t minSize = 0 );

        void detectErroneousComponent (double covCutoff, size_t maxMargLength , vector<Component> &tobeRemoved);
        struct less_than_Component
        {
                inline bool operator() (const Component& object1, const Component& object2)
                {
                        return (object1.componentSize < object2.componentSize);
                }
        };
        struct greater_than_Component
        {
                inline bool operator() (const Component& object1, const Component& object2)
                {
                        return (object1.componentSize > object2.componentSize);
                }
        };

};



