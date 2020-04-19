#ifndef SHASTA_APPROXIMATE_COLORING_HPP
#define SHASTA_APPROXIMATE_COLORING_HPP

/*******************************************************************************

Approximate coloring of a Boost undirected graph using the Dsatur algorithm.

See this Wikipedia article for greedy coloring in general:
https://en.wikipedia.org/wiki/Greedy_coloring

This code uses the Brelaz (1979) Dsatur
algorithm described there in this section
https://en.wikipedia.org/wiki/Greedy_coloring#Adaptive

The main Wikipedia article for this algorithm is here:
https://en.wikipedia.org/wiki/DSatur

The Dsatur algorithm was ingtroduced in this paper:
Brelaz, Daniel (April 1979),
"New methods to color the vertices of a graph",
Communications of the ACM, 22 (4): 251–256, doi:10.1145/359094.359101

Many alternative approximate coloring methods are also available, for example
https://www.gerad.ca/~alainh/RLFPaper.pdf

For a partial coloring of a graph, the Dsatur algorithm defines
the saturation degree of a vertex as the number of distinct colors
to which the vertex is adjacent (only counting, of course,
vertices that have already been colored).
With this definition, the Dsatur algorithm is simply described
as follows:

- Starting with the uncolored graph, iterate over vertices in the order
  defined below.
- At each iteration, assign to each vertex the lowest possible color.
  This is the lowest color that the vertex is not adjacent to.
- At each iteration choose the vertex with maximum saturation degree
  (maximum number of adjacent colors). In case of ties, break
  the ties by selecting the vertex with maximum degree in the uncolored
  subgraph - that is, the vertex with the greatest number of uncolored
  adjacent vertices.

The DSatur algoritm is an approximate coloring algorithm, but it
is exact for 2-colorable graphs.

This function applies the above algorithm to a Boost undirected graph.
It assumes that the graph is defined as a
boost::adjacency_list<
    OutEdgeList,
    VertexList,
    Directed,
    VertexProperties, EdgeProperties, GraphProperties,
    EdgeList>
with:

- VertexList = boost::vecS (so vertex descriptors are integers starting at zero).
- Directed = boost::undirectedS.

It returns the computed coloring in its second argument, indexed
by vertex_descriptor.

*******************************************************************************/

// Shasta.
#include "SHASTA_ASSERT.hpp"

// Boost libraries.
#include <boost/container/flat_set.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "algorithm.hpp"
#include "cstdint.hpp"
#include <map>
#include <numeric>
#include "vector.hpp"

namespace shasta {
    template<class Graph> void approximateColoring(
        const Graph&,
        vector<uint64_t>&
    );
}



template<class Graph> void shasta::approximateColoring(
    const Graph& graph,
    vector<uint64_t>& color
)
{
    // Get the number of vertices and initialize the colors.
    const uint64_t n = boost::num_vertices(graph);
    color.resize(n);
    const uint64_t noColor = std::numeric_limits<uint64_t>::max();
    fill(color.begin(), color.end(), noColor);

    // Vector to contain the number of colored vertices adjacent to
    // each uncolored vertex. This is used to break ties in the Dsatur algorithm.
    // We only maintain this for uncolored vertices.
    // Once a vertex gets colored we stop updating this as it will never be needed.
    vector<uint64_t> adjacentColoredCount(n, 0);

    // The colors adjacent to each uncolored vertex. The number of colors adjacent
    // to each vertex is small, so we use boost::container::flat_set
    // instead of std::set.
    // We only maintain this for uncolored vertices.
    // Once a vertex gets colored we stop updating this as it will never be needed.
    using boost::container::flat_set;
    vector< flat_set<uint64_t> > adjacentColors(n);

    // A multimap that stores uncolored vertices
    // keyed by saturation degree (that is, the number of colors adjacent to
    // each vertex), in decreasing order.
    // This is used to select, at each iteration, the next vertex to be colored.
    using MultiMap = std::multimap<uint64_t, uint64_t, std::greater<uint64_t> >;
    MultiMap verticesBySaturationDegree;
    for(uint64_t v=0; v<n; v++) {
        verticesBySaturationDegree.insert(make_pair(0, v));
    }



    // Main iteration loop. At each iteration we color one vertex.
    for(uint64_t iteration=0; iteration<n; iteration++) {

        // To choose the vertex to color at this iteration, loop over
        // all uncolored vertices with maximum saturation degree.
        SHASTA_ASSERT(not verticesBySaturationDegree.empty());
        const uint64_t highestSaturationDegree = verticesBySaturationDegree.begin()->first;
        MultiMap::iterator begin, end;
        tie(begin, end) = verticesBySaturationDegree.equal_range(highestSaturationDegree);
        SHASTA_ASSERT(begin != end);
        uint64_t v0 = begin->second;
        uint64_t adjacentColoredCount0 = adjacentColoredCount[v0];
        MultiMap::iterator it0 = begin;
        for(MultiMap::iterator it=begin; it!=end; ++it) {
            const uint64_t v1 = it->second;
            if(adjacentColoredCount[v1] > adjacentColoredCount0) {
                v0 = v1;
                adjacentColoredCount0 = adjacentColoredCount[v1];
                it0 = it;
            }
        }

        // Remove v0 from verticesBySaturationDegree.
        verticesBySaturationDegree.erase(it0);

        // This iteration is coloring LocalVertexId v0.
        // We choose the lowest color that is not present in adjacentColors[v0].
        const auto& adjacentColors0 = adjacentColors[v0];
        uint64_t color0;
        for(uint64_t color=0; ; color++) {
            if(adjacentColors0.find(color) == adjacentColors0.end()) {
                color0 = color;
                break;
            }
        }

        // Color v0 with the color chosen in this way.
        color[v0] = color0;

        // Increment adjacentColoredCount for adjacent uncolored vertices.
        BGL_FORALL_ADJ_T(v0, v1, graph, Graph) {
            if(color[v1] == noColor) {
                adjacentColoredCount[v1]++;
            }
        }


        // Update adjacentColors for adjacent uncolored vertices, and
        // update verticesBySaturationDegree accordingly.
        BGL_FORALL_ADJ_T(v0, v1, graph, Graph) {
            if(color[v1] == noColor) {
                auto& adjacentColors1 = adjacentColors[v1];
                if(adjacentColors1.find(color0) == adjacentColors1.end()) {
                    adjacentColors1.insert(color0);

                    // The saturation degree of this vertex has increased.
                    MultiMap::iterator begin, end;
                    tie(begin, end) = verticesBySaturationDegree.equal_range(adjacentColors1.size()-1);
                    SHASTA_ASSERT(begin != end);
                    bool done = false;
                    for(MultiMap::iterator it=begin; it!=end; ++it) {
                        if(it->second == v1) {
                            verticesBySaturationDegree.erase(it);
                            verticesBySaturationDegree.insert(make_pair(adjacentColors1.size(), v1));
                            done = true;
                            break;
                        }
                    }
                    SHASTA_ASSERT(done);
                }
            }
        }
    }

    // Verify that we colored all the vertices.
    SHASTA_ASSERT(find(color.begin(), color.end(), noColor) == color.end());

}

#endif

