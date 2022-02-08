#ifndef SHASTA_SHORTEST_PATH_BOUNDED_DISTANCE_HPP
#define SHASTA_SHORTEST_PATH_BOUNDED_DISTANCE_HPP

// #include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/iteration_macros.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include "vector.hpp"

namespace shasta {
    template<class Graph, class DistanceType>
        DistanceType shortestPathBoundedDistance(
            const Graph&,
            typename Graph::vertex_descriptor source,
            typename Graph::vertex_descriptor destination,
            DistanceType maxPathLength,
            const std::map<typename Graph::edge_descriptor, DistanceType>& edgeLength,
            vector<typename Graph::edge_descriptor>& path
        );

    void testShortestPathBoundedDistance();
}



// Compute a shortest path between two given vertices
// of a directed graph, subject to a maximum path length.
// If there is no such path, an empty path is computed.
// The graph must be a directed graph described
// as a boost::graph::adjacency_list.
// A map with edge lengths must be given.
// Edge lengths must be non-negative.
// All data structures used in this function are of size
// proportional to size of the portion of the graph within a distance maxPathlength
// from the source vertex.
// It may be possible to implement this using boost::dijkstra_shortest_paths
// with a suitable visitor.

// See pseudocode here
// https://en.wikipedia.org/wiki/Dijkstra's_algorithm#Pseudocode

// The code mostly follows that pseudocode with a few differences:
// - Instead of maintaining Q, with maintain its complement notQ.
//   This way all data structures are of size proportional
//   to the size of the graph neighborhood visited.
// - The prev map stores the edge, not the vertex, that led to the
//   shortest path for each vertex encountered.
// - Distances larger than maxPathLength are never stored.
// - Returns the length of the path found.
template<class Graph, class DistanceType> DistanceType
    shasta::shortestPathBoundedDistance(
        const Graph& graph,
        typename Graph::vertex_descriptor source,
        typename Graph::vertex_descriptor destination,
        DistanceType maxPathLength,
        const std::map<typename Graph::edge_descriptor, DistanceType>& edgeLength,
        vector<typename Graph::edge_descriptor>& path
    )
{
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;

    // The Boost MultiIndex library is used to implement the
    // updatable priority queue required by the Dijkstra algorithm.
    using boost::multi_index::indexed_by;
    using boost::multi_index::multi_index_container;
    using boost::multi_index::ordered_unique;
    using boost::multi_index::ordered_non_unique;
    using boost::multi_index::member;
    using boost::multi_index::tag;

    // The complement of the Q vertex set.
    std::set<vertex_descriptor> notQ;

    // The prev map contain for each vertex the edge
    // that gave the minimum distance for that vertex.
    std::map<vertex_descriptor, edge_descriptor> prev;

    // The vertices and their distances, accessible both by vertex descriptor
    // and by distance.
    using VertexInfo = pair<vertex_descriptor, DistanceType>;
    using VertexContainer = multi_index_container< VertexInfo, indexed_by<
        ordered_unique<     member<VertexInfo, vertex_descriptor, &VertexInfo::first > >,
        ordered_non_unique< member<VertexInfo, DistanceType     , &VertexInfo::second> > > >;
    VertexContainer vertexContainer;

    // Vertex indexes by vertex descriptor and by distance.
    // That use of the template keyword is weird but necessary.
    auto& verticesByDescriptor = vertexContainer.template get<0>();
    auto& verticesByDistance = vertexContainer.template get<1>();

    // Insert the source vertex.
    vertexContainer.insert(VertexInfo(source, 0));



    // Main loop.
    while(not vertexContainer.empty()) {

        // Find the vertex non in notQ with minimum distance.
        vertex_descriptor u = Graph::null_vertex();
        DistanceType uDistance = std::numeric_limits<DistanceType>::max();
        for(const auto& uInfo: verticesByDistance) {
            if(notQ.find(uInfo.first) != notQ.end()) {
                continue;
            }
            u = uInfo.first;
            uDistance = uInfo.second;
            break;
        }
        if(u == Graph::null_vertex()) {
            break;
        }

         // Instead of removing u from Q (see the pseudocode reference above),
         // add u to not notQ.
         notQ.insert(u);

         // If done, construct the path back to the source.
         if(u == destination) {
             path.clear();
             DistanceType pathLength = DistanceType(0ULL);
             while(u != source) {
                 const edge_descriptor e = prev[u];
                 path.push_back(e);
                 auto ite = edgeLength.find(e);
                 SHASTA_ASSERT(ite != edgeLength.end());
                 pathLength += ite->second;;
                 u = boost::source(e, graph);
             }
             std::reverse(path.begin(), path.end());
             return pathLength;
         }

         // Loop over children of u that are not in notQ.
         BGL_FORALL_OUTEDGES_T(u, e, graph, Graph) {
             const vertex_descriptor v = target(e, graph);
             if(notQ.find(v) != notQ.end()) {
                 // It is notQ. Skip.
                 continue;
             }

             // Compute the tentative new distance of v.
             auto ite = edgeLength.find(e);
             SHASTA_ASSERT(ite != edgeLength.end());
             const DistanceType alt = uDistance + ite->second;
             if(alt > maxPathLength) {
                 continue;
             }

             // If we found a lower distance, update it.
             auto itv = verticesByDescriptor.find(v);
             if(itv == verticesByDescriptor.end()) {
                 vertexContainer.insert(VertexInfo(v, alt));
                 prev[v] = e;
             } else {
                 const DistanceType vDistance = itv->second;
                 if(alt < vDistance) {
                     verticesByDescriptor.replace(itv, VertexInfo(v, alt));
                     prev[v] = e;
                 }
             }
         }

    }

    path.clear();
    return DistanceType(0ULL);
}


inline void shasta::testShortestPathBoundedDistance()
{
    using namespace boost;
    using Graph = adjacency_list<listS, listS, bidirectionalS, uint64_t>;
    using vertex_descriptor = Graph::vertex_descriptor;
    using edge_descriptor = Graph::edge_descriptor;
    using DistanceType = uint64_t;

    Graph graph;

    // Create the vertices.
    vector<vertex_descriptor> v;
    const uint64_t n = 4;
    for(uint64_t i=0; i<n; i++) {
        v.push_back(add_vertex(i, graph));
    }

    // Edges and their lengths.
    vector< pair< pair<uint64_t, uint64_t>, DistanceType> > edgeTable = {
        {{0, 1}, 1},
        {{0, 2}, 1},
        {{1, 3}, 10},
        {{2, 3}, 4},
        {{0, 3}, 1000},
    };

    // Create the edges.
    vector<edge_descriptor> edges;
    std::map<edge_descriptor, DistanceType> edgeLength;
    for(const auto& p: edgeTable) {
        edge_descriptor e;
        tie(e, ignore) = add_edge(v[p.first.first], v[p.first.second], graph);
        edgeLength.insert(make_pair(e, p.second));
    }


    // Compute the shortest path.
    const uint64_t maxPathLength = 5;
    vector<edge_descriptor> path;
    const uint64_t pathLength = shortestPathBoundedDistance(graph, v[0], v[3], maxPathLength, edgeLength, path);
    cout << "Found a path with " << path.size() << " edges of length " << pathLength << endl;\
    for(const edge_descriptor e: path) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        cout << graph[v0] << " " << graph[v1] << endl;
    }

}


#endif
