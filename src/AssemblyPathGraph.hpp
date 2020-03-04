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

Note that an assembly graph edge can appear in multiple
edges of the assembly path graph.

The assembly path graph is represented as a boost::adjacency_list
from the Boost Graph library. This makes it easy to manipulate the
assembly graph path. However this also means that the assembly path
graph cannot be stored in memory mapped files like many other
data structures.



TANGLES AND THEIR PROPERTIES

A tangle in the assembly path graph is defined as an edge v0->v1
with source vertex v0 and target vertex v1 such that:

- In-degree(v0)>1, out-degree(v0)=1
- In-degree(v1)=1, out-degree(v1)>1

In words, a tangle is a "bottleneck" where a number of incoming branches
all converge and from which a number of outgoing branches all diverge.

For example, here is an illustration of a tangle with 2 incoming edges
and 2 outgoing edges:

               v0            v1
________________________________________________
               /              \
______________/                \________________


Some nomenclature used below:

- Edge v0->v1 is called the tangle edge.
- The incoming edges of v0 are called the in-edges of the tangle.
- The outgoing edges of v1 are called the out-edges of the tangle.
- The number of in-edges of the tangle is called the in-degree of the tangle.
- The number of out-edges of the tangle are called the out-degree of the tangle.

In a typical tangle, the in-degree and out-degree are the same,
but this is not necessary.

Some properties:

- The tangle edge of a tangle cannot also be the tangle edge of another tangle,
- The tangle edge of a tangle cannot also be an in-edge or out-edge of another tangle.
- An in-edge of a tangle cannot also be an in-edge of another tangle.
- An out-edge of a tangle cannot also be an out-edge of another tangle.
- However, an in-edge of a tangle can also be an out-edge of another tangle.
- And conversely, an out-edge of a tangle can also be an in-edge of another tangle.

Because of these properties, the following holds:
- Any edge can be the tangle edge of zero or one tangles.
- An edge can be an in-edge of zero or one tangles.
- An edge can be an out-edge of zero or one tangles.

And we use the following additional nomenclature:
- If an edge is an in-edge of a tangle, we call that tangle the out-tangle of the edge.
- If an edge is an out-edge of a tangle, we call that tangle the in-tangle of the edge.

This means that for each edge we can store up to 3 tangles:
- The tangle the edge creates, if the edge is a tangle edge of a tangle.
- The in-tangle.
- The out-tangle.
These are stored as TangleId's and set to invalidTangleId when
there is no such tangle.

***************************************************************************/



// Shasta.
#include "AssemblyGraph.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "algorithm.hpp"
#include "iosfwd.hpp"
#include <map>
#include "string.hpp"
#include "vector.hpp"


namespace shasta {
    class AssemblyPathGraph;
    class AssemblyPathGraphVertex;
    class AssemblyPathGraphEdge;
    class Tangle;

    using TangleId = uint64_t;
    static const TangleId invalidTangleId = std::numeric_limits<TangleId>::max();

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

    AssemblyPathGraphBaseClass::vertex_descriptor reverseComplementVertex =
        AssemblyPathGraphBaseClass::null_vertex();
};



class shasta::AssemblyPathGraphEdge {
public:

    // The AsssemblyGraph path corresponding to this edge.
    vector <AssemblyGraph::EdgeId> path;

    // The length of the path, as measured on the marker graph.
    uint64_t pathLength = 0;

    // The tangles that this edge participates in.
    // These are set to invalidTangleId if missing.
    // See above for nomenclature.
    TangleId tangle = invalidTangleId;
    TangleId inTangle = invalidTangleId;
    TangleId outTangle = invalidTangleId;
    void clearTangles()
    {
        tangle = invalidTangleId;
        inTangle = invalidTangleId;
        outTangle = invalidTangleId;
    }

    // The reverse complement of this edge.
    AssemblyPathGraphBaseClass::edge_descriptor reverseComplementEdge;


    // Initialize the path to a single AssemblyGraph edge.
    AssemblyPathGraphEdge(AssemblyGraph::EdgeId edgeId) :
        path(1, edgeId) {}
    AssemblyPathGraphEdge() {}

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

    void mergeOrientedReadIds(
        const vector<OrientedReadId>&,
        const vector<OrientedReadId>&
        );
};

inline std::ostream& shasta::operator<<(
    ostream& s,
    const AssemblyPathGraphEdge& edge)
{
    s << string(edge);
    return s;
}



class shasta::Tangle {
public:
    using vertex_descriptor = AssemblyPathGraphBaseClass::vertex_descriptor;
    using edge_descriptor = AssemblyPathGraphBaseClass::edge_descriptor;

    TangleId tangleId;
    edge_descriptor edge;
    vector<edge_descriptor> inEdges;
    vector<edge_descriptor> outEdges;
    bool unsolvable = false;

    uint64_t inDegree() const
    {
        return inEdges.size();
    }
    uint64_t outDegree() const
    {
        return outEdges.size();
    }

    // The tangle matrix stores the number of common reads between
    // each pair of inEdges and outEdges.
    // Indexed by [i][j] where i is an index into inEdges and j
    // ins an index into outEdges.
    vector< vector<uint64_t> > matrix;
    bool hasZeroMatrixElements() const;
};



class shasta::AssemblyPathGraph : public AssemblyPathGraphBaseClass {
public:

    // The constructor does not fill in the oriented read ids for each edge.
    // This must be done separately (see Assembler::detangle).
    AssemblyPathGraph(const AssemblyGraph&);

    // The tangles currently present in the graph, keyed by their ids.
    TangleId nextTangleId = 0;
    std::map<TangleId, Tangle> tangles;
    Tangle& getTangle(TangleId);
    const Tangle& getTangle(TangleId) const;
    TangleId getReverseComplementTangle(TangleId) const;
    void removeTangle(TangleId);

    void fillReverseComplementEdge(edge_descriptor);

    // Initial creation of all tangles.
    void createTangles();

    // Create tangles involving a given edge.
    // This can create up to two tangles involving
    // the given edge as an in-edge, out-edge, or tangle edge.
    // This is used for incrementally create new tangles as
    // edges are created during detangling.
    void createTanglesInvolvingEdge(edge_descriptor e);

    // Create a new tangle that has the specified edge
    // as the tangle edge, if such a tangle is valid
    // and does not already exist.
    // Return true if the new tangle was created.
    bool createTangleAtEdge(edge_descriptor e);

    // Return the next tangle to work on.
    TangleId findNextTangle() const;

    // Return true if a tangle collides with its reverse complement.
    bool collidesWithReverseComplement(TangleId) const;

    // Detangle all we can.
    // The average number of bases per marker is only used
    // for GFA output.
    void detangle(double basesPerMarker);

    // Detangle a single tangle.
    // This does not fill in the reverseComplementEdge of newly created edges,
    // and does not create new tangles involving those edges.
    bool detangle(TangleId, vector<edge_descriptor>& newEdges);

    // Output in Graphviz format.
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

    // Html output.
    void writeHtml(const string& fileName) const;
    void writeHtml(ostream&) const;
    void writeVerticesHtml(ostream&) const;
    void writeEdgesHtml(ostream&) const;
    void writeTanglesHtml(ostream&) const;

    // GFA output (without sequence).
    void writeGfa(const string& fileName, double basesPerMarker) const;
    void writeGfa(ostream&, double basesPerMarker) const;

private:
    void removeIsolatedVertices();
};




#endif
