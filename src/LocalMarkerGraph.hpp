#ifndef CZI_NANOPORE2_LOCAL_MARKER_GRAPH_HPP
#define CZI_NANOPORE2_LOCAL_MARKER_GRAPH_HPP


/*******************************************************************************

The local marker graph is a directed graph in which each vertex
corresponds to a marker in an oriented read.
Consecutive markers in each oriented read are joined by an edge,
so each oriented read generates
a linear sequence of vertices in the graph.

Therefore, the local marker graph initially consists of a
separate linear sequence for each oriented read.

Then, for each marker alignment found between the oriented reads,
we merge pairs of vertices corresponding to the aligned markers.

This results in a graph that is a representation of
the alignment of the oriented reads.

This class can only handle a small number of oriented reads,
and therefore can only be used for local assemblies,
for testing/debugging. Separate code will be necessary to handle
the global marker graph.

The given oriented reads are indexed using a local oriented read id.

*******************************************************************************/

// Nanopore2.
#include "Kmer.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>
#include "tuple.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {

        // Forward declaration of types declared in this file.
        class LocalMarkerGraphVertex;
        class LocalMarkerGraphEdge;
        class LocalMarkerGraph;
        using LocalMarkerGraphBaseClass = boost::adjacency_list<
            boost::setS,
            boost::listS,
            boost::bidirectionalS,
            LocalMarkerGraphVertex,
            LocalMarkerGraphEdge
            >;

        // Forward declarations of types defined elsewhere.
        class LongBaseSequence;
        class Marker;
    }
}



class ChanZuckerberg::Nanopore2::LocalMarkerGraphVertex {
public:

    // The k-mer id of the marker corresponding to this vertex.
    KmerId kmerId;

    // The markers that were merged into this vertex.
    class MarkerId {
    public:
        uint32_t localOrientedReadId;
        uint32_t ordinal;
        MarkerId(uint32_t localOrientedReadId, uint32_t ordinal) :
            localOrientedReadId(localOrientedReadId), ordinal(ordinal) {}
        bool operator<(const MarkerId& that) const
        {
            return tie(localOrientedReadId, ordinal) <
                tie(that.localOrientedReadId, that.ordinal);
        }
    };
    vector<MarkerId> markerIds;

    // A vertex id used for graphviz output.
    uint32_t vertexId;

    // Color, used for topological sorting.
    uint32_t color;
};



class ChanZuckerberg::Nanopore2::LocalMarkerGraphEdge {
public:
    class Data {
    public:
        uint32_t localOrientedReadId;
        array<uint32_t, 2> ordinals;
        Data(uint32_t localOrientedReadId, uint32_t ordinal0, uint32_t ordinal1) :
            localOrientedReadId(localOrientedReadId),
            ordinals(array<uint32_t, 2>{ordinal0, ordinal1}) {}
        int shift() const
        {
            return int(ordinals[1]) - int(ordinals[0]);
        }
        bool operator<(const Data& that) const
        {
            return localOrientedReadId < that.localOrientedReadId;
        }
        bool operator==(const Data& that) const
        {
            return localOrientedReadId == that.localOrientedReadId &&
                ordinals == that.ordinals;
        }
    };
    vector<Data> data;
    bool isInOptimalSpanningTree = false;
};



class ChanZuckerberg::Nanopore2::LocalMarkerGraph :
    public LocalMarkerGraphBaseClass {
public:

    // Construct the graph given the markers of
    // given oriented reads, sorted by position.
    // Each pair in the second input vector is (position, KmerId).
    // This constructs the initial graph without any
    // alignment information, consisting of a
    // separate linear sequence for each input oriented read.
    LocalMarkerGraph(
        size_t k,
        const vector<OrientedReadId>&,
        const vector<LongBaseSequence>& sequences,
        const vector< vector<Marker> >&,
        size_t minCoverage,     // For a vertex to be considered strong.
        size_t minConsensus     // For an edge to be considered strong.
        );

    // Merge two vertices.
    // The vertices to be merged are specified by
    // the local oriented read ids (indexes in the orientedReadIds vector)
    // and the ordinals.
    void mergeVertices(
        uint32_t localOrientedReadId0,
        uint32_t ordinal0,
        uint32_t localOrientedReadId1,
        uint32_t ordinal1
    );

    // Remove vertices with coverage less than minCoverage.
    void removeWeakVertices();

    // Remove edges with consensus less than the minConsensus
    // (consensus is maximum number of oriented reads that agree on the sequence).
    void removeWeakEdges();

    // Recursively prune weak leaves.
    void pruneWeakLeaves();

    // Same as above, but keep all edges marked as spanning tree edges.
    void removeWeakNonSpanningTreeEdges();

    // Sort the markers in each vertex.
    void sort();

    // Sort the markers in each vertex
    // and remove or split ambiguous vertices.
    // A vertex is ambiguous if it has more
    // than one occurrence for at least one oriented read.
    void sortAndRemoveAmbiguousVertices();
    void sortAndSplitAmbiguousVertices();

    // Fill in edge data.
    // This assumes that there are no ambiguous vertices
    // and that the occurrences vectors in the vertices are sorted.
    // This must be done before calling write.
    // This only fills in edge data when the
    // ordinal in the source and target vertices
    // are contiguous.
    void fillEdgeData();

    // Add edges that were removed because of low coverage,
    // skipping if necessary (we can only use surviving vertices).
    void addMissingEdges();

    // Mark edges that are part of an optimal spanning tree.
    void computeOptimalSpanningTree();

    // Extract the longest sequence present in the graph.
    // For each position the returned sequence contains a base and a coverage value.
    using Base = Nanopore2::Base;   // Override boost::graph definition.
    vector< pair<Base, int> > extractLongestSequence();

    // Write in Graphviz format.
    void write(ostream&, bool addEdgeLabels) const;
    void write(const string& fileName, bool addEdgeLabels) const;

    // Check that the graph and the vertex map are consistent.
    void check();

private:

    // The length of the markers we are using.
    size_t k;

    // The OrientedReadId's used by this assembly graph.
    vector<OrientedReadId>  orientedReadIds;

    // Their base sequences.
    const vector<LongBaseSequence> sequences;

    // The markers the input oriented reads, sorted by position.
    vector< vector<Marker> > markers;

    // Minimum coverage for a vertex to be considered strong.
    size_t minCoverage;

    // Minimum consensus for an edge to be considered strong.
    size_t minConsensus;

    // The vertex map is indexed by [localOrientedReadId][ordinal].
    vector< vector<vertex_descriptor> > vertexMap;
    uint32_t nextVertexId = 0;

    // Split a vertex, creating a vertex for each k-mer occurrence.
    void split(vertex_descriptor);

    // Fill in edge data.
    // This assumes that there are no ambiguous vertices
    // and that the occurrences vectors in the vertices are sorted.
    // This must be done before calling write.
    void fillEdgeData(edge_descriptor);

    // Compute the consensus of an edge.
    // It is the maximum number of oriented read ids
    // that agree on the sequence associated with the edge.
    size_t edgeConsensus(edge_descriptor) const;

    // Construct the sequence implied by an LocalMarkerGraphEdge::Data.
    // If the source and target k-mers overlap,
    // the string is a number representing the number of overlap bases.
    // Otherwise it is a base sequence.
    // We should phase this out in favor of getEdgeSequence.
    string getSequence(const LocalMarkerGraphEdge::Data&) const;



    // Class to describe the sequence implied by an LocalMarkerGraphEdge::Data.
    // If the edge sequences overlap, sequence is empty and overlap
    // contains the number of overlapping bases.
    // Otherwise, overlap is zero and sequence contain the sequence
    // between the two vertices.
    // If sequence is empty and overlap is zero, the sequences
    // of the two vertices immediately follow each other
    // without any overlap and without any intervening bases.
    class EdgeSequence {
    public:
        vector<Base> sequence;
        int overlap = 0;
        bool operator<(const EdgeSequence& that) const
        {
            return tie(overlap, sequence) < tie(that.overlap, that.sequence);
        }
    };
    EdgeSequence getEdgeSequence(const LocalMarkerGraphEdge::Data&) const;

    // Get the dominant sequence of an edge and its coverage.
    pair<EdgeSequence, int> getDominantEdgeSequence(edge_descriptor) const;



    // Find the longest path in the graph.
    // This assumes that the graph is acyclic.
    vector<edge_descriptor> findLongestPath();

    // Topological sort of the graph.
    // This assumes that the graph is acyclic.
    vector<vertex_descriptor> topologicalSort();

    // Extract the sequence corresponding to a path.
    // For each position the returned sequence contains a base and a coverage value.
    vector< pair<Base, int> > getPathSequence(const vector<edge_descriptor>&);

    class Writer {
    public:
        Writer(const LocalMarkerGraph&, bool addEdgeLabels);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const LocalMarkerGraph& graph;
        bool addEdgeLabels;
    };
    friend class Writer;
};

#endif
