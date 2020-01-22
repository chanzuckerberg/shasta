// Shasta.
#include "ConflictReadGraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/container/flat_set.hpp>

// Standard library.
#include <map>
#include "fstream.hpp"



/*******************************************************************************

See this Wikipedia article for greedy coloring in general:
https://en.wikipedia.org/wiki/Greedy_coloring

This code uses the Brelaz (1979) Dsatur
algorithm described there in this section
https://en.wikipedia.org/wiki/Greedy_coloring#Adaptive

The main Wikipedia article for this algorithm is here:
https://en.wikipedia.org/wiki/DSatur

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

This function applies the above algorithm to the specified set of vertices,
which are assumed to be a connected component of the graph
and stored in sorted order.

*******************************************************************************/

void ConflictReadGraph::colorConnectedComponent(const vector<VertexId>& component)
{
    SHASTA_ASSERT(std::is_sorted(component.begin(), component.end()));

    // A LocalVertexId is an index into the component vector.
    using Int = VertexId;

    // The number of vertices in this connected component.
    const Int n = component.size();
    SHASTA_ASSERT(n > 1);

    // Vector to contain the number of colored vertices adjacent to
    // each uncolored vertex. This is used to break ties in the Dsatur algorithm.
    // We only maintain this for uncolored vertices.
    // Once a vertex gets colored we stop updating this as it will never be needed.
    vector<Int> adjacentColoredCount(n, 0);

    // The colors adjacent to each uncolored vertex. The number of colors adjacent
    // to each vertex is small, so we use boost::container::flat_set
    // instead of std::set.
    // We only maintain this for uncolored vertices.
    // Once a vertex gets colored we stop updating this as it will never be needed.
    using boost::container::flat_set;
    vector< flat_set<Int> > adjacentColors(n);

    // A multimap that stores uncolored vertices (as local vertex ids)
    // keyed by saturation degree (that is, the number of colors adjacent to
    // each vertex), in decreasing order.
    // This is used to select, at each iteration, the next vertex to be colored.
    using MultiMap = std::multimap<Int, Int, std::greater<Int> >;
    MultiMap verticesBySaturationDegree;
    for(Int v=0; v<n; v++) {
        verticesBySaturationDegree.insert(make_pair(0, v));
    }

    // Make sure all vertices have no cluster id.
    for(Int v=0; v<n; v++) {
        getVertex(component[v]).clusterId = ConflictReadGraphVertex::invalid;
    }



    // Main iteration loop. At each iteration we color one vertex.
    for(Int iteration=0; iteration<n; iteration++) {
        // cout << "Iteration " << iteration << " begins." << endl;

        // To choose the vertex to color at this iteration, loop over
        // all uncolored vertices with maximum saturation degree.
        SHASTA_ASSERT(not verticesBySaturationDegree.empty());
        const Int highestSaturationDegree = verticesBySaturationDegree.begin()->first;
        MultiMap::iterator begin, end;
        tie(begin, end) = verticesBySaturationDegree.equal_range(highestSaturationDegree);
        SHASTA_ASSERT(begin != end);
        Int v0 = begin->second;
        Int adjacentColoredCount0 = adjacentColoredCount[v0];
        MultiMap::iterator it0 = begin;
        for(MultiMap::iterator it=begin; it!=end; ++it) {
            const Int v1 = it->second;
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
        Int color0 = ConflictReadGraphVertex::invalid;
        for(Int color=0; ; color++) {
            if(adjacentColors0.find(color) == adjacentColors0.end()) {
                color0 = color;
                break;
            }
        }

        // Color v0 with the color chosen in this way.
        const VertexId u0 = component[v0];
        // cout << "Iteration " << iteration << " colored vertex " << v0 << " " << getOrientedReadId(u0) <<
        //     " color " << color0 << endl;
        ConflictReadGraphVertex& vertex0 = getVertex(u0);
        SHASTA_ASSERT(vertex0.clusterId = ConflictReadGraphVertex::invalid);
        vertex0.clusterId = color0;

        // Increment adjacentColoredCount for adjacent uncolored vertices.
        for(EdgeId edgeId: incidentEdges(u0)) {
            const VertexId u1 = otherVertex(edgeId, u0);
            const ConflictReadGraphVertex& vertex1 = getVertex(u1);
            if(vertex1.clusterId == ConflictReadGraphVertex::invalid) {
                const auto it = lower_bound(component.begin(), component.end(), u1);
                SHASTA_ASSERT(*it == u1);
                const Int v1 = it - component.begin();
                adjacentColoredCount[v1]++;
            }
        }
        // cout << "***A" << endl;


        // Update adjacentColors for adjacent uncolored vertices, and
        // update verticesBySaturationDegree accordingly.
        for(EdgeId edgeId: incidentEdges(u0)) {
            const VertexId u1 = otherVertex(edgeId, u0);
            const ConflictReadGraphVertex& vertex1 = getVertex(u1);
            if(vertex1.clusterId == ConflictReadGraphVertex::invalid) {
                const auto it = lower_bound(component.begin(), component.end(), u1);
                SHASTA_ASSERT(*it == u1);
                const Int v1 = it - component.begin();
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
        // cout << "***B" << endl;

    }


    // Find the number of colors used.
    Int maxColor = 0;
    for(Int v=0; v<n; v++) {
        maxColor = max(maxColor, getVertex(component[v]).clusterId);
    }
    cout << "Used " << maxColor+1 << " colors." << endl;

}



void ConflictReadGraph::writeGraphviz(const string& fileName) const
{
    ofstream s(fileName);
    writeGraphviz(s);
}
void ConflictReadGraph::writeGraphviz(ostream& s) const
{
    s <<
        "graph G {\n" <<
        "node [shape=point];\n";

    // Write the vertices.
    for(VertexId v=0; v<vertices.size(); v++) {
        const OrientedReadId orientedReadId = getOrientedReadId(v);
        s << v << "[tooltip=\"" << orientedReadId << "\"];\n";
    }

    // Write the edges.
    for(EdgeId e=0; e<edges.size(); e++) {
        s << v0(e) << "--" << v1(e) << ";\n";
    }

    s << "}\n";
}
