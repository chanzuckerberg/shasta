#ifndef SHASTA_SHORTEST_PATH_HPP
#define SHASTA_SHORTEST_PATH_HPP

// Function to find the shortest weighted path
// between two vertices of an undirected graph.
// We cannot use boost::dijkstra_shortest_paths
// from the Boost Graph library because we want to be
// able to reuse the priority queue, to
// reduce memory allocation activity.

// Requirements on class Graph:
// - Must support a subset of the boost::graph::adjacency_list API.
// - The vertex must have the following data members:
//   vertex_descriptor predecessor;
//   uint64_t distance;
//   uint8_t color;
// - The edge must have a uint64_t weight data member
//   which defines path lengths.



// shasta
#include "orderPairs.hpp"

// Boost Graph library.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "cstddef.hpp"
#include <queue>
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {

        // The last argument to findShortestPath is a work area with this type.
        template<class Graph> using FindShortestPathQueue  =
            std::priority_queue<
                pair< uint64_t, typename Graph::vertex_descriptor>,
                vector< pair< uint64_t, typename Graph::vertex_descriptor> >,
                OrderPairsByFirstOnlyGreater< uint64_t, typename Graph::vertex_descriptor>
            >;

        template<class Graph> void findShortestPath(
            Graph&,
            typename Graph::vertex_descriptor vSource,
            typename Graph::vertex_descriptor vTarget,
            vector<typename Graph::vertex_descriptor>& path,

            // Work area. Does not need to be initialized.
            // When calling findShortestPath repeatedly,
            // use the same FindShortestPathQueue to reduce memory
            // allocation activity.
            FindShortestPathQueue<Graph>&
            );
    }
}



template<class Graph> inline void ChanZuckerberg::shasta::findShortestPath(
    Graph& graph,
    typename Graph::vertex_descriptor vSource,
    typename Graph::vertex_descriptor vTarget,
    vector<typename Graph::vertex_descriptor>& path,
    FindShortestPathQueue<Graph>& q
    )
{
    // Trivial special case.
    if(vTarget == vSource) {
        path.clear();
        path.push_back(vSource);
        return;
    }

    using vertex_descriptor = typename Graph::vertex_descriptor;

    // Initialize.
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        auto& vertex = graph[v];
        vertex.predecessor = Graph::null_vertex();
        vertex.distance = std::numeric_limits<uint64_t>::max();
        vertex.color = 0;
    }
    graph[vSource].predecessor = vSource;
    graph[vSource].distance = 0;
    while(!q.empty()) {
        q.pop();
    }
    q.push(make_pair(0, vSource));



    // Main loop.
    while(!q.empty()) {

        // Dequeue the closest vertex in the queue.
        const auto p0 = q.top();
        q.pop();
        const uint64_t distance0 = p0.first;
        const vertex_descriptor v0 = p0.second;
        // cout << "Dequeued " << v0.v << " at distance " << distance0 << endl;

        // If already encountered, skip.
        // This can happen because the inner loop uses "lazy deletion".
        // For example, see here:
        // https://stackoverflow.com/questions/9209323/easiest-way-of-using-min-priority-queue-with-key-update-in-c
        auto& vertex0 = graph[v0];
        if(vertex0.color == 1) {
            // cout << "Already encountered, skipped." << endl;
            continue;
        }
        vertex0.color = 1;

        // If we found vLast, construct the path and be done.
        if(v0 == vTarget) {
            path.clear();
            vertex_descriptor v = v0;
            while(true) {
                path.push_back(v);
                if(v == vSource) {
                    break;
                }
                v = graph[v].predecessor;
            }
            std::reverse(path.begin(), path.end());
            while(!q.empty()) {
                q.pop();
            }
            return;
        }

        // Loop over its out-edges.
        BGL_FORALL_OUTEDGES_T(v0, e01, graph, Graph) {
            const vertex_descriptor v1 = target(e01, graph);
            auto& vertex1 = graph[v1];
            if(vertex1.color == 1) {
                // cout << "    " << v1.v << " skipped because of color." << endl;
                continue;
            }
            const uint64_t weight = graph[e01].weight;
            const uint64_t distance1 = distance0 + weight;
            // cout << "    Found " << v1.v << " at distance " << distance1 << endl;

            if(distance1 < vertex1.distance) {
                q.push(make_pair(distance1, v1));
                vertex1.predecessor = v0;
                vertex1.distance = distance1;
                // cout << "        Enqueued." << endl;
            }
        }
    }

    // If getting here, the queue is empty but we have not found vLast.
    // This means that there is no path between vFirst and vlast.
    path.clear();
}




#endif
