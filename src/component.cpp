
#include "component.h"
#include "set"
using namespace std;

void ComponentHandler:: findComponentsInGraph(size_t minSize)
{
        cout << "Finding disjoint components in the graph ..." <<endl;
        int srcID=1;
        int i=0;
        set<NodeID> nodesHandled;
        for (i =srcID ; i <= dbg.getNumNodes(); i++){
                set<NodeID> currentSetNodes;
                size_t componentSize =0;
                if (!dbg.getSSNode(i).isValid())
                        continue;
                if (nodesHandled.find(i) != nodesHandled.end())
                        continue;
                if (nodesHandled.find(i) != nodesHandled.end())
                        continue;
                srcID=i;
                multimap<size_t, NodeID> nodeDepth;     // map of all nodes in the local graph and its depth
                nodeDepth.insert(pair<size_t, NodeID>(0, srcID));
                // nodes that were already handled

                while (!nodeDepth.empty()) {
                        // get and erase the current node
                        multimap<size_t, NodeID>::iterator
                        e = nodeDepth.begin();
                        size_t thisDepth = e->first;
                        NodeID thisID = e->second;
                        nodeDepth.erase(e);
                        // if the node was already handled, skip
                        if (nodesHandled.find(thisID) != nodesHandled.end())
                                continue;
                        if (nodesHandled.find(-thisID) != nodesHandled.end())
                                continue;
                        // mark this node as handled
                        nodesHandled.insert(thisID);
                        currentSetNodes.insert(thisID);
                        SSNode thisNode = dbg.getSSNode(thisID);
                        componentSize = componentSize + thisNode.getMarginalLength();
                        for (ArcIt it = thisNode.rightBegin(); it != thisNode.rightEnd(); it++) {
                                SSNode rNode = dbg.getSSNode(it->getNodeID());
                                if (!rNode.isValid())
                                        continue;
                                if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                        continue;
                                nodeDepth.insert(pair<size_t, NodeID>(thisDepth + thisNode.getMarginalLength(), it->getNodeID()));
                        }
                        for (ArcIt it = thisNode.leftBegin(); it != thisNode.leftEnd(); it++) {
                                SSNode lNode =dbg.getSSNode(it->getNodeID());
                                if (!lNode.isValid())
                                        continue;
                                if (nodesHandled.find(it->getNodeID()) != nodesHandled.end())
                                        continue;
                                nodeDepth.insert(pair<size_t, NodeID>(thisDepth + thisNode.getMarginalLength(), it->getNodeID()));
                        }
                }
                if ( componentSize >minSize && currentSetNodes.size() >1)
                        this->addComponent(currentSetNodes);

                currentSetNodes.clear();
        }
        reportStatistics();
}
void ComponentHandler::reportStatistics ()
{

        vector<Component> components;
        for (auto it : componentsMap)
                components.push_back(it.second);
        cout << "There are " << components.size() << "number of disjoint components in the graph" <<endl;
        /*sort(components.begin(), components.end(),greater_than_Component());
        for (auto  it : components)
        {
            Component c = it;
            cout <<"Component ID   : " <<c.componentID << endl;
            cout <<"Coverage       : " <<c.componentKmerCov << endl;
            cout <<"Num Of Nodes   : " <<c.numOfNodes <<endl;
            cout <<"Component Size : " <<c.componentSize <<endl;
            cout <<"**************" <<endl;
        }*/
}
void ComponentHandler::detectErroneousComponent (double covCutoff, size_t maxMargLength , vector<Component> &tobeRemoved)
{
        findComponentsInGraph();
        vector<Component> components;
        for (auto it : componentsMap)
                components.push_back(it.second);
        for (auto  it : components)
        {
            Component c = it;
            if (c.componentKmerCov < covCutoff && c.componentSize < maxMargLength)
                    tobeRemoved.push_back(c);
        }
}

