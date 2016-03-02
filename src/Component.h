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
        std::map<int, size_t> N;
        std::set<NodeID> NodesID;

public:
        size_t Size;
        size_t numOfNodes;
        size_t numOfArcs;
        double nodeKmerCov;
        size_t largestNodeSize;
        DBGraph  &dbg;
        Settings &settings;


        Component( DBGraph& g,  Settings& s, set<NodeID> &IDs) :  N(), Size(0),numOfNodes(0),numOfArcs(0),nodeKmerCov(0), largestNodeSize(), dbg(g),settings(s) {setNodeIDSet(IDs);}
        std::map<int, size_t> get_N() const {
                return N;
        }
        /**
         * set the value of Nx
         * @param x the size in Nx can be 10, 20, ...100
         * @param Nx the value associated to the size
         **/
        void set_N(int x, size_t Nx) {
                N[x] = Nx;
        }
        /**
         * returns the value of Nx
         * @param x the size of Nx, can be 10, 20, ...100
         * @return the value of Nx
         */
        size_t get_N(int x) const {
                return N.at(x);
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
        /**
         * makes essential files for Cytoscape  program to show the Component
         * @param ID referes to the name of the component, lower number shows bigger size.
         *
         **/
        void writeCytoscapeComponent( size_t ID);
        /**
         * reports statistical information about the component
         **/
        void getComponentSta();

         /**
          * writes back to the file the Nx of the component (x :10, 20, ...100)
          * @param size id or the name for the component, lower number are bigger components
          *
          **/
         void makeNxFiles(const size_t num);
         /**
          * makes a file shows frequency distribution of nodes in a component
          * based on thier node Kmer Coverage
          * @param size id or the name for the component, lower number are bigger components
          */
         void plotCovDiagram(size_t num);


};


class ComponentHandler{
private:
        vector<Component> components;
        DBGraph &dbg;
        Settings &settings;
public:
        size_t numberOfcomponents;
        ComponentHandler( DBGraph& g,Settings& s):dbg(g),settings(s),numberOfcomponents(0){}

        /**
         * Add new coponent to the list of comonent
         * @param nodesIdSet set of nodes which are in the component
         * @return make a new component and returns it
         */
        Component addComponent(set<NodeID> nodesIdSet){
                Component newComonent(dbg,settings, nodesIdSet);
                components.push_back(newComonent);
                numberOfcomponents++;
                return newComonent;
        }
        /**
         * sorts Components list in a descending order by their size
         * @param components vector of unsorted component should be sorted
         *
         */

        void  sortDes(vector<Component>& components){
                sort( components.begin(), components.end(), greater_than_Component());
        }
        /**
         * sorts Components list in a ascending order by their size
         * @param components vector of unsorted component should be sorted
         *
         */
        void sortAsc(vector<Component>&components){
                sort( components.begin(), components.end(), less_than_Component());
        }

        /**
         * returns the sorted list of all components in the graph in descending order
         * @return vector of component
         *
         */
        vector<Component> getAllComponentsSortedDes(){
                sortDes(components);
                return components;
        }
         /**
         * returns the sorted list of all components in the graph in ascending order
         * @return vector of component
         *
         */
        vector<Component> getAllComponentsSortedAsc(){
                sortAsc(components);
                return components;
        }
        /**
         * finds the frequency Of component in different size.
         * normally there are more small components and few big
         * @param components set of all components in graph
         *
         **/

        void makeComponentPlotFile();

        struct less_than_Component
        {
                inline bool operator() (const Component& object1, const Component& object2)
                {
                        return (object1.Size < object2.Size);
                }
        };
        struct greater_than_Component
        {
                inline bool operator() (const Component& object1, const Component& object2)
                {
                        return (object1.Size > object2.Size);
                }
        };

};
