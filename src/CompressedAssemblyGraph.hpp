#ifndef SHASTA_COMPRESSED_ASSEMBLY_GRAPH_HPP
#define SHASTA_COMPRESSED_ASSEMBLY_GRAPH_HPP

/*******************************************************************************

The compressed assembly graph is a representation of the assembly graph
in which each linear sequence of bubbles is compressed to a single edge.

*******************************************************************************/

#include "AssemblyGraph.hpp"
#include <boost/graph/adjacency_list.hpp>

namespace shasta {
    class CompressedAssemblyGraph;
    class CompressedAssemblyGraphEdge;
    class CompressedAssemblyGraphVertex;

    using CompressedAssemblyGraphBaseClass =
        boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        CompressedAssemblyGraphVertex,
        CompressedAssemblyGraphEdge
        >;

    class Assembler;
}


class shasta::CompressedAssemblyGraphVertex {
public:
    AssemblyGraph::VertexId vertexId;

    CompressedAssemblyGraphVertex(AssemblyGraph::VertexId vertexId) :
        vertexId(vertexId) {}
};



class shasta::CompressedAssemblyGraphEdge {
public:

    // The chain of assembly graph vertices associated
    // with this edge.
    // This includes the assembly graph vertices
    // associated with the source and target of this edge.
    vector<AssemblyGraph::VertexId> vertices;

    // The chain of sets of parallel assembly graph edges.
    vector< vector<AssemblyGraph::EdgeId> > edges;

    // An id assigned to this edge of the compressed assembly graph
    // and used in gfa and other output.
    uint64_t id;
    string gfaId() const;

    // The minimum and maximum marker count (path length) in markers.
    uint64_t minMarkerCount;
    uint64_t maxMarkerCount;
    void fillMarkerCounts(const AssemblyGraph&);

    // Find the oriented reads that appear in marker graph vertices
    // internal to this edge of the compressed assembly graph.
    void findOrientedReads(const Assembler&);
    vector<OrientedReadId> orientedReadIds;
    vector<uint64_t> orientedReadIdsFrequency;

private:

    // Append to orientedReadIds the oriented reads that
    // appear in a given marker graph edge.
    void findOrientedReads(
        const Assembler&,
        const MarkerGraph::EdgeId&);
};



class shasta::CompressedAssemblyGraph :
    public CompressedAssemblyGraphBaseClass {
public:
    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::EdgeId;


    // Create the CompressedAssemblyGraph from the AssemblyGraph.
    CompressedAssemblyGraph(
        const Assembler&);

    // GFA output (without sequence).
    void writeGfa(const string& fileName, double basesPerMarker) const;
    void writeGfa(ostream&, double basesPerMarker) const;

    // HTML output.
    void writeHtml(const string& fileName) const;
    void writeHtml(ostream&) const;

private:

    // Create a vertex for each vertex of the assembly graph.
    void createVertices(
        uint64_t vertexCount,
        vector<vertex_descriptor>& vertexTable);

    // Create an edge for each set of parallel edges of the assembly graph.
    void createEdges(
        const AssemblyGraph&,
        const vector<vertex_descriptor>& vertexTable
    );

    // Merge linear chains of edges.
    void mergeLinearChains();

    // Assign an id to each edge.
    void assignEdgeIds();

    // Fill in the assembly graph edges that go into each
    // edge of the compressed assembly graph.
    void fillContributingEdges(const AssemblyGraph&);

    // Fill in minimum and maximum marker counts for each edge.
    void fillMarkerCounts(const AssemblyGraph&);

    // Find the oriented reads that appear in marker graph vertices
    // internal to each edge of the compressed assembly graph.
    void findOrientedReads(const Assembler&);
};


#endif
