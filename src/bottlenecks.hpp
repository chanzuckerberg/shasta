#ifndef SHASTA_BOTTLENECKS_HPP
#define SHASTA_BOTTLENECKS_HPP

#include <boost/graph/iteration_macros.hpp>

#include <queue>
#include <set>
#include "vector.hpp"

namespace shasta {

    // Given a directed graph and a start vertex,
    // find vertices (the "bottlenecks") that are on all paths
    // starting at the start vertex.
    template<class Graph> void findForwardBottlenecks(
        const Graph&,
        typename Graph::vertex_descriptor startVertex,
        vector<typename Graph::vertex_descriptor>& bottlenecks);
}



// Given a directed graph and a start vertex,
// find vertices (the "bottlenecks") that are on all paths
// starting at the start vertex.
// This is done using a BFS starting at the start vertex.
// The bottlenecks are the vertex that, when dequeued,
// leave the queue empty.
template<class Graph> void shasta::findForwardBottlenecks(
    const Graph& graph,
    typename Graph::vertex_descriptor startVertex,
    vector<typename Graph::vertex_descriptor>& bottlenecks)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;

    bottlenecks.clear();

    std::queue<vertex_descriptor> q;
    q.push(startVertex);

    std::set<vertex_descriptor> discoveredVertices;
    discoveredVertices.insert(startVertex);

    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();
        if(q.empty() and ( v0!= startVertex)) {
            bottlenecks.push_back(v0);
        }

        BGL_FORALL_OUTEDGES_T(v0, e, graph, Graph) {
            const vertex_descriptor v1 = target(e, graph);
            bool wasInserted = false;
            tie(ignore, wasInserted) = discoveredVertices.insert(v1);
            if(wasInserted) {
                q.push(v1);
            }
        }
    }


}




#endif

