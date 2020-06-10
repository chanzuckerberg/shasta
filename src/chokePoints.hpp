#ifndef SHASTA_CHOKE_POINTS_HPP
#define SHASTA_BOTTLENECKS_HPP

#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/dominator_tree.hpp>
#include <boost/graph/reverse_graph.hpp>

// #include "iostream.hpp"
#include <list>
#include <map>
#include "vector.hpp"

namespace shasta {

    // Given a directed graph and a start vertex,
    // find vertices (the "choke points") that are on all paths
    // starting at the start vertex.
    template<class Graph> void findForwardChokePoints(
        const Graph&,
        typename Graph::vertex_descriptor startVertex,
        vector<typename Graph::vertex_descriptor>& chokePoints);

    // Same, but in the opposite direction.
    template<class Graph> void findBackwardChokePoints(
        const Graph&,
        typename Graph::vertex_descriptor startVertex,
        vector<typename Graph::vertex_descriptor>& chokePoints);
}



// Given a directed graph and a start vertex,
// find vertices (the "choke points") that are on all paths
// starting at the start vertex.
// The Graph must have a vertex_index property
// (which in practice means its second template parameter must be vecS).
// They are computed as follows:
// - Compute the dominator tree for the given start vertex.
// - For all vertices with out-degree 0 (in the graph,
//   not the dominator tree), walk up the dominator tree
//   up to the start vertex.
// - The common portion of the walks up the dominator tree
//   gives the choke points.
template<class Graph> void shasta::findForwardChokePoints(
    const Graph& graph,
    typename Graph::vertex_descriptor startVertex,
    vector<typename Graph::vertex_descriptor>& chokePoints)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using boost::lengauer_tarjan_dominator_tree;
    using boost::make_assoc_property_map;

    // Compute the dominator tree.
    std::map<vertex_descriptor, vertex_descriptor> predecessorMap;
    lengauer_tarjan_dominator_tree(
        graph,
        startVertex,
        make_assoc_property_map(predecessorMap));

    /*
    cout << "Dominator tree: " << endl;
    for(const auto& p: predecessorMap) {
        cout << graph[p.first].segmentId << " " << graph[p.second].segmentId << endl;
    }
    */



    // Follow the tree up from all terminal vertices - that is,
    // vertices with out-degree 0.
    // Find the common portion of the path up.
    vector<vertex_descriptor> commonPath;
    vector<vertex_descriptor> path;
    bool firstTime = true;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {

        // If this is not a terminal vertex, skip.
        if(out_degree(v, graph) != 0) {
            continue;
        }

        // If there is no predecessor, it means that this
        // vertex is not accessible from the start vertex. Skip.
        if(predecessorMap.find(v) == predecessorMap.end()) {
            continue;
        }

        // cout << "Walk up the tree from " << graph[v].segmentId << endl;

        // Walk the tree up.
        path.clear();
        while(v != startVertex) {
            path.push_back(v);
            v = predecessorMap[v];
        }
        reverse(path.begin(), path.end());

        /*
        cout << "Path:";
        for(const vertex_descriptor u: path) {
            cout << " " << graph[u].segmentId;
        }
        cout << endl;
        */


        // Find the common portion.
        if(firstTime) {
            commonPath = path;
            firstTime = false;
        } else {
            uint64_t index = 0;
            for(; index<path.size() and index<commonPath.size(); index++) {
                if(path[index] != commonPath[index]) {
                    // cout << "Found a difference at index " << index << endl;
                    break;
                }
            }
            commonPath.resize(index);
        }

    }
    chokePoints = commonPath;
}



template<class Graph> void shasta::findBackwardChokePoints(
    const Graph& graph,
    typename Graph::vertex_descriptor startVertex,
    vector<typename Graph::vertex_descriptor>& chokePoints)
{
    findForwardChokePoints(
        boost::make_reverse_graph(graph),
        startVertex,
        chokePoints);
}

#endif

