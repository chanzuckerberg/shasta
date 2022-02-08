#ifndef SHASTA_MAKE_BICONNECTED_HPP
#define SHASTA_MAKE_BICONNECTED_HPP

// Boost libraries.
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <vector>



// Recursively remove articulation points from a graph
// until there are none left, at which point the graph is biconnected.
// Also computes a map containing the index of the biconnected component
// that each edge belongs to.

namespace shasta {
    template<class Graph> void makeBiconnected(
        Graph& graph,
        std::map<typename Graph::edge_descriptor, uint64_t>& componentMap);
}



template<class Graph> void shasta::makeBiconnected(
    Graph& graph,
    std::map<typename Graph::edge_descriptor, uint64_t>& componentMap)
{
    using boost::vertex_index_map;
    using boost::make_assoc_property_map;
    using vertex_descriptor = typename Graph::vertex_descriptor;

    std::map<vertex_descriptor, uint64_t> vertexMap;
    vector<vertex_descriptor> articulationPoints;
    while(true) {

        vertexMap.clear();
        uint64_t vertexIndex = 0;
        BGL_FORALL_VERTICES_T(v, graph, DynamicConflictReadGraph) {
            vertexMap.insert(make_pair(v, vertexIndex++));
        }
        articulationPoints.clear();
        componentMap.clear();
        boost::biconnected_components(
            graph,
            make_assoc_property_map(componentMap),
            back_inserter(articulationPoints),
            vertex_index_map(make_assoc_property_map(vertexMap)));

        if(articulationPoints.empty()) {
            return;
        }

        // Remove the articulation points.
        for(const vertex_descriptor v: articulationPoints) {
            clear_vertex(v, graph);
            remove_vertex(v, graph);
        }

    }

}


#endif
