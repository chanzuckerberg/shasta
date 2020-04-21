#ifndef SHASTA_SUBGRAPH_HPP
#define SHASTA_SUBGRAPH_HPP

// Function to create a local subgraph of a given Boost directed graph.
// The local subgraph is a copy of the original graph that
// only includes vertices within a specified distance from a
// given set of start vertices, and edges between those vertices.

#include <boost/bimap.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <map>
#include <queue>
#include "vector.hpp"

namespace shasta {

    template<class DirectedGraph> void createLocalSubgraph(

        // The input graph.
        const DirectedGraph& graph,

        // The start vertices. Only vertices within the specified distance
        // from one of these vertices are kept in the the subgraph.
        const vector<typename DirectedGraph::vertex_descriptor>& startVertices,

        // The maximum distance from these start vertices.
        uint64_t maxDistance,

        // The subgraph created by this function. It must
        // initially be empty.
        DirectedGraph& subgraph,

        // Vertex map.
        // vertexMap.left gives the vertex descriptor in the graph for a given
        // vertex descriptor in the subgraph.
        // vertexMap.right gives the vertex descriptor in the subgraph for a given
        // vertex descriptor in the graph.
        boost::bimap<
            typename DirectedGraph::vertex_descriptor,
            typename DirectedGraph::vertex_descriptor>& vertexMap,

        // Edge map.
        // edgeMap.left gives the edge descriptor in the graph for a given
        // edge descriptor in the subgraph.
        // edgeMap.right gives the edge descriptor in the subgraph for a given
        // edge descriptor in the graph.
        boost::bimap<
            typename DirectedGraph::edge_descriptor,
            typename DirectedGraph::edge_descriptor>& edgeMap,

        // Distance map for the subgraph.
        // It gives the distance from the start vertices
        // of each vertex of the subgraph.
        std::map<
            typename DirectedGraph::vertex_descriptor,
            uint64_t>& distanceMap
        );

}



template<class DirectedGraph> inline void shasta::createLocalSubgraph(

    // The input graph.
    const DirectedGraph& graph,

    // The start vertices. Only vertices within the specified distance
    // from one of these vertices are kept in the the subgraph.
    const vector<typename DirectedGraph::vertex_descriptor>& startVertices,

    // The maximum distance from these start vertices.
    uint64_t maxDistance,

    // The subgraph created by this function. It must
    // initially be empty.
    DirectedGraph& subgraph,

    // Vertex map.
    // vertexMap.left gives the vertex descriptor in the graph for a given
    // vertex descriptor in the subgraph.
    // vertexMap.right gives the vertex descriptor in the subgraph for a given
    // vertex descriptor in the graph.
    boost::bimap<
        typename DirectedGraph::vertex_descriptor,
        typename DirectedGraph::vertex_descriptor>& vertexMap,

    // Edge map.
    // edgeMap.left gives the edge descriptor in the graph for a given
    // edge descriptor in the subgraph.
    // edgeMap.right gives the edge descriptor in the subgraph for a given
    // edge descriptor in the graph.
    boost::bimap<
        typename DirectedGraph::edge_descriptor,
        typename DirectedGraph::edge_descriptor>& edgeMap,

    // Distance map for the subgraph.
    // It gives the distance from the start vertices
    // of each vertex of the subgraph.
    std::map<
        typename DirectedGraph::vertex_descriptor,
        uint64_t>& distanceMap
    )
{
    using vertex_descriptor = typename DirectedGraph::vertex_descriptor;
    using edge_descriptor = typename DirectedGraph::edge_descriptor;

    using VertexMapValueType =
        typename boost::bimap<
        typename DirectedGraph::vertex_descriptor,
        typename DirectedGraph::vertex_descriptor>::value_type;
    using EdgeMapValueType =
        typename boost::bimap<
        typename DirectedGraph::edge_descriptor,
        typename DirectedGraph::edge_descriptor>::value_type;

    // Vertex queue used for the BFS.
    // It contains vertices of the full graph.
    std::queue<vertex_descriptor> q;


    // We use "v" to indicate vertices of the full graph and "u"
    // for vertices of the subgraph.

    // Initialize the BFS with our starting vertices.
    for(const vertex_descriptor v: startVertices) {
        if(maxDistance > 0) {
            q.push(v);
        }
        const vertex_descriptor u = boost::add_vertex(graph[v], subgraph);
        vertexMap.insert(VertexMapValueType(u, v));
        distanceMap.insert(make_pair(u, 0));
    }



    // Do the BFS.
    while(not q.empty()) {

        // Dequeue a vertex of the full graph.
        vertex_descriptor v0 = q.front();
        q.pop();

        // Find the corresponding vertex of the subgraph.
        const auto it0 = vertexMap.right.find(v0);
        SHASTA_ASSERT(it0 != vertexMap.right.end());
        const vertex_descriptor u0 = it0->second;

        // Find its distance.
        const auto jt0 = distanceMap.find(u0);
        SHASTA_ASSERT(jt0 != distanceMap.end());
        const uint64_t distance0 = jt0->second;
        const uint64_t distance1 = distance0 + 1;

        // Loop over children.
        BGL_FORALL_OUTEDGES_T(v0, e01, graph, DirectedGraph) {
            const vertex_descriptor v1 = boost::target(e01, graph);
            if(vertexMap.right.find(v1) == vertexMap.right.end()) {
                const vertex_descriptor u1 = boost::add_vertex(graph[v1], subgraph);
                vertexMap.insert(VertexMapValueType(u1, v1));
                distanceMap.insert(make_pair(u1, distance1));
                if(distance1 < maxDistance) {
                    q.push(v1);
                }
            }
        }

        // Loop over parents.
        BGL_FORALL_INEDGES_T(v0, e01, graph, DirectedGraph) {
            const vertex_descriptor v1 = boost::source(e01, graph);
            if(vertexMap.right.find(v1) == vertexMap.right.end()) {
                const vertex_descriptor u1 = boost::add_vertex(graph[v1], subgraph);
                vertexMap.insert(VertexMapValueType(u1, v1));
                distanceMap.insert(make_pair(u1, distance1));
                if(distance1 < maxDistance) {
                    q.push(v1);
                }
            }
        }

    }



    // Add the edges.
    for(const auto& p: vertexMap.right) {
        const vertex_descriptor v0 = p.first;
        const vertex_descriptor u0 = p.second;
        BGL_FORALL_OUTEDGES_T(v0, e01, graph, DirectedGraph) {
            const vertex_descriptor v1 = boost::target(e01, graph);
            const auto it1 = vertexMap.right.find(v1);
            if(it1 != vertexMap.right.end()) {
                const vertex_descriptor u1 = it1->second;
                edge_descriptor e01New;
                bool edgeWasAdded;
                tie(e01New, edgeWasAdded) = boost::add_edge(u0, u1, graph[e01], subgraph);
                SHASTA_ASSERT(edgeWasAdded);
                edgeMap.insert(EdgeMapValueType(e01New, e01));
            }
        }
    }
}


#endif
