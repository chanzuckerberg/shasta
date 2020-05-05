#ifndef SHASTA_COMPRESSED_ASSEMBLY_GRAPH_HPP
#define SHASTA_COMPRESSED_ASSEMBLY_GRAPH_HPP

/*******************************************************************************

The compressed assembly graph is a representation of the assembly graph
in which each linear sequence of bubbles is compressed to a single edge.

*******************************************************************************/

#include "AssemblyGraph.hpp"

#include <boost/bimap.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "array.hpp"
#include <map>

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
    double averageMarkerCount() const
    {
        return double(minMarkerCount + maxMarkerCount) / 2.;
    }
    void fillMarkerCounts(const AssemblyGraph&);

    // Find the oriented reads that appear in marker graph vertices
    // internal to this edge of the compressed assembly graph.
    void findOrientedReads(const Assembler&);
    vector<OrientedReadId> orientedReadIds;
    vector<uint64_t> orientedReadIdsFrequency;

    // The edges that have at least one oriented read in common
    // with this edge.
    vector<CompressedAssemblyGraphBaseClass::edge_descriptor>
        relatedEdges;

    // Return the maximum number of branches in a bubble.
    uint64_t maxPloidy() const;

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
    CompressedAssemblyGraph(const Assembler&);

    // Create a local subgraph.
    // See createLocalSubgraph for argument explanation.
    CompressedAssemblyGraph(
        const CompressedAssemblyGraph& graph,
        const Assembler&,
        const vector<vertex_descriptor>& startVertices,
        uint64_t maxDistance,
        boost::bimap<vertex_descriptor, vertex_descriptor>& vertexMap,
        boost::bimap<edge_descriptor, edge_descriptor>& edgeMap,
        std::map<vertex_descriptor, uint64_t>& distanceMap
        );

    // Return the edge with a given GFA id.
    pair<edge_descriptor, bool> getEdgeFromGfaId(const string&) const;

    // GFA output (without sequence).
    void writeGfa(const string& fileName, double basesPerMarker) const;
    void writeGfa(ostream&, double basesPerMarker) const;

    // Graphviz output.
    void computeVertexLayout(
        uint64_t sizePixels,
        double vertexScalingFactor,
        double edgeLengthPower,
        double edgeLengthScalingFactor,
        double timeout,
        std::map<vertex_descriptor, array<double, 2 > >& vertexPositions
    ) const;
    void writeGraphviz(
        const string& fileName,
        uint64_t sizePixels,
        double vertexScalingFactor,
        double edgeLengthScalingFactor,
        double edgeThicknessScalingFactor,
        double edgeArrowScalingFactor,
        std::map<vertex_descriptor, array<double, 2 > >& vertexPositions) const;
    void writeGraphviz(
        ostream&,
        uint64_t sizePixels,
        double vertexScalingFactor,
        double edgeLengthScalingFactor,
        double edgeThicknessScalingFactor,
        double edgeArrowScalingFactor,
        std::map<vertex_descriptor, array<double, 2 > >& vertexPositions) const;

    // Dump everything to csv files.
    void writeCsv() const;

    // Create a csv file with coloring.
    // If the string passed in is an oriented read,
    // it colors all edges that have that read.
    // If it is the gfaId of an edge, it colors that edge in red
    // and all related edges in green.
    // This can be loaded in Bandage to color the edges.
    void color(
        const string&,
        const string& fileName) const;
    void color(
        const string&,
        ostream&) const;

private:
    void writeCsvEdges() const;
    void writeCsvBubbleChains() const;
    void writeCsvOrientedReadsByEdge() const;
    void writeCsvOrientedReads() const;

    // Create a vertex for each vertex of the assembly graph.
    void createVertices(
        uint64_t vertexCount,
        vector<vertex_descriptor>& vertexTable);

    // Create an edge for each set of parallel edges of the assembly graph.
    void createEdges(
        const AssemblyGraph&,
        const vector<vertex_descriptor>& vertexTable
    );

    // Remove back edges that create reverse bubbles.
    void removeReverseBubbles();

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

    // The edges that each oriented read appears in.
    // Indexed by OrientedRead::getValue().
    vector< vector<edge_descriptor> > orientedReadTable;
    void fillOrientedReadTable(const Assembler&);

    // Find edges that have at least one common oriented read
    // which each edge.
    void findRelatedEdges();
    void findRelatedEdges(edge_descriptor);

private:
    uint64_t maxPloidy() const;

};


#endif
