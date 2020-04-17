#ifndef SHASTA_SUBGRAPH_HPP
#define SHASTA_SUBGRAPH_HPP

// Function to create a local subgraph of a given Boost directed graph.
// The local subgraph is a copy of the original graph that
// only includes vertices within a specified distance from a
// given set of start vertices, and edges between those vertices.

#include <map>
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
        // Gives the vertex descriptor in the graph for a given
        // vertex descriptor in the subgraph.
        std::map<
            typename DirectedGraph::vertex_descriptor,
            typename DirectedGraph::vertex_descriptor>& vertexMap,

        // Distance map for the subgraph.
        // It gives the distance from the start vertices
        // of each vertex of the subgraph.
        std::map<
            typename DirectedGraph::vertex_descriptor,
            uint64_t>& distanceMap
        );

}

#endif
