#ifndef SHASTA_ASSEMBLY_PATH_GRAPH_HPP
#define SHASTA_ASSEMBLY_PATH_GRAPH_HPP



/***************************************************************************

The assembly path graph is a directed graph in which edge corresponds
to a path in the assembly graph. The path is described as an ordered
sequence of assembly graph edges.

Each vertex of the assembly path graph corresponds to a vertex
of the assembly graph, but not all assembly graph vertices
have a corresponding vertex in the assembly path graph.

The assembly path graph is used for detangling the assembly graph.

When the assembly path graph is first created, each edge stores
a path consisting of a single assembly graph edge.
However, as detangling proceeds, individual edges are combined
into longer paths.

Note that an assembly graph edhe can appear in multiple
edges of the assembly path graph.

The assembly path graph is represented as a boost::adjacency_list
from the Boost Graph library. This makes it easy to manipulate the
assembly graph path. However this also means that the assembly path
graph cannot be stored in memory mapped files like many other
data structures.

***************************************************************************/

// Shasta.
#include "AssemblyGraph.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "algorithm.hpp"
#include "iosfwd.hpp"
#include "string.hpp"
#include "vector.hpp"


namespace shasta {
    class AssemblyPathGraph;
    class AssemblyPathGraphVertex;
    class AssemblyPathGraphEdge;

    using AssemblyPathGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        AssemblyPathGraphVertex,
        AssemblyPathGraphEdge
        >;
    inline ostream& operator<<(
        ostream&,
        const AssemblyPathGraphEdge&);

    class AssemblyGraph;
}



class shasta::AssemblyPathGraphVertex {
public:
    AssemblyGraph::VertexId vertexId;
    AssemblyPathGraphVertex(AssemblyGraph::VertexId vertexId) :
        vertexId(vertexId) {}
};



class shasta::AssemblyPathGraphEdge {
public:

    // The AsssemblyGraph path corresponding to this edge.
    vector <AssemblyGraph::EdgeId> path;

    // The length of the path, as measured on the marker graph.
    uint64_t pathLength;

    // Initialize the path to a single AssemblyGraph edge.
    AssemblyPathGraphEdge(AssemblyGraph::EdgeId edgeId) :
        path(1, edgeId) {}

    // The OrientedReadId's on this path, sorted.
    vector<OrientedReadId> orientedReadIds;

    // Represent it as a string consisting of the path edge ids,
    // separated by dashes.
    operator string() const
    {
        string s;
        for(uint64_t i=0; i<path.size(); i++) {
            s += to_string(path[i]);
            if(i != path.size()-1) {
                s += "-";
            }
        }
        return s;
    }
};

inline std::ostream& shasta::operator<<(
    ostream& s,
    const AssemblyPathGraphEdge& edge)
{
    s << string(edge);
    return s;
}



class shasta::AssemblyPathGraph : public AssemblyPathGraphBaseClass {
public:
    AssemblyPathGraph(const AssemblyGraph&);
    void detangle();

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

};


#endif
