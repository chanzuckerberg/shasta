#ifndef SHASTA_DOMINATOR_TREE_HPP
#define SHASTA_DOMINATOR_TREE_HPP

/*******************************************************************************


This provides a bug fix for the 3-argument version of
boost::lengauer_tarjan_dominator_tree.
The fixed version is in the shasta namespace.

The bug fix is as follows.

The original version contains

    std::vector<VerticesSizeType> dfnum(numOfVertices, 0);

The fixed version below has instead

    std::vector<VerticesSizeType> dfnum(numOfVertices, std::numeric_limits<VerticesSizeType>::max());

As a result of the bug, the original version produces incorrect results
if the graph contains vertices that are unreachable from the entrance.

lengauer_tarjan_dominator_tree_without_dfs, which is eventually called,
includes the following comment about initializing dfnum/dfnumMap:

   * @pre Unreachable nodes must be masked as
   *      (std::numeric_limits<VerticesSizeType>::max)() in dfnumMap.

The above fix implements this requirement.
The fixed version below works correctly even if the graph contains unreachable vertices.

*******************************************************************************/

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dominator_tree.hpp>

namespace shasta {

    template<class Graph, class DomTreePredMap>
    void lengauer_tarjan_dominator_tree(const Graph &g,
        const typename boost::graph_traits<Graph>::vertex_descriptor &entry,
        DomTreePredMap domTreePredMap)
    {
        using namespace boost;

        // typedefs
        typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename graph_traits<Graph>::vertices_size_type VerticesSizeType;
        typedef typename property_map<Graph, vertex_index_t>::const_type IndexMap;
        typedef iterator_property_map<
            typename std::vector<VerticesSizeType>::iterator, IndexMap> TimeMap;
        typedef iterator_property_map<typename std::vector<Vertex>::iterator,
            IndexMap> PredMap;

        // Make property maps
        const VerticesSizeType numOfVertices = num_vertices(g);
        if (numOfVertices == 0)
            return;

        const IndexMap indexMap = get(vertex_index, g);

        // ********* THE BUG FIX IS HERE **********************
        // The original versionwas filling the vector with 0.
        std::vector<VerticesSizeType> dfnum(numOfVertices,
            std::numeric_limits < VerticesSizeType > ::max());

        TimeMap dfnumMap(make_iterator_property_map(dfnum.begin(), indexMap));

        std::vector<Vertex> parent(numOfVertices,
            graph_traits < Graph > ::null_vertex());
        PredMap parentMap(make_iterator_property_map(parent.begin(), indexMap));

        std::vector<Vertex> verticesByDFNum(parent);

        // Run main algorithm
        boost::lengauer_tarjan_dominator_tree(g, entry, indexMap, dfnumMap,
            parentMap, verticesByDFNum, domTreePredMap);
    }

}

#endif

