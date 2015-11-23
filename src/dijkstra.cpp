#include "dijkstra.h"

Dijkstra::Dijkstra(const DBGraph &g):dbg(g)
{
}

/* if the node isn't exist in set, add it,
 * if exist, update the value, if the path make it shorter change the value
 * otherwise keep the previous value
 *
 *
 */
void Dijkstra::updateChecklist(NodeID nodeID, double newDistance)
{
        bool find=false;
        for(auto it :checkedList){
                if (it.first==nodeID){
                        find=true;
                        it.second=newDistance<it.second? newDistance:it.second;
                }
        }
        if (!find)
                checkedList.push_back(make_pair(newDistance,nodeID ));

}
double Dijkstra::getDistance(NodeID id )
{
        for(auto it :checkedList)
                if (it.second==id)
                        return it.first;
        return 0;
}


double Dijkstra::shortestPath(SSNode start, SSNode end)
{
        element root;
        root.first=0;
        root.second=start.getNodeID();
        visited.clear();
        checkedList.clear();
        checkedList.push_back( root);
        //double min=std::numeric_limits<double>::infinity();
        while(!checkedList.empty()){
                element expand=checkedList[0];
                SSNode expandNode=dbg.getSSNode( expand.second);
                 for ( ArcIt it = expandNode.rightBegin(); it != expandNode.rightEnd(); it++ ) {
                        SSNode rightNode=dbg.getSSNode(it->getNodeID());
                        //new distance is equal to the cost of reaching to this point + the length of new node
                        if (visited.find(it->getNodeID())==visited.end())
                                updateChecklist(rightNode.getNodeID(),rightNode.getMarginalLength()  +getDistance(rightNode.getNodeID()));
                 }
                 if (checkedList[0].second==end.getNodeID())
                         return( checkedList[0].first);
                 else
                         visited.insert(checkedList[0].second);
                 checkedList.erase(checkedList.begin());
                 std::sort (checkedList.begin(), checkedList.end());
        }
        return -1;
}
