
#include "component.h"
#include "set"
using namespace std;

size_t ComponentHandler:: findComponentsInGraph(size_t minSize)
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
        cout << "Number of disjoint components in the graph with more than 1 nodes and larger than : " <<minSize << " is " <<this->numberOfcomponents <<endl;
        return this->numberOfcomponents ;
}
