#ifndef SHASTA_LOCAL_MARKER_GRAPH_HPP
#define SHASTA_LOCAL_MARKER_GRAPH_HPP

/*******************************************************************************

The local marker graph created by class LocalMarkerGraph is a subgraph
of the global marker graph, created by starting at a given vertex,
and extending out to a specified distance in both directions.
Distance is number of edges on the global marker graph.

Like in the global marker graph, each vertex corresponds to
a group of aligned markers.

*******************************************************************************/

// Shasta.
#include "AssemblyGraph.hpp"
#include "Coverage.hpp"
#include "Kmer.hpp"
#include "LocalMarkerGraphRequestParameters.hpp"
#include "MarkerGraph.hpp"
#include "Reads.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

namespace shasta {

    class LocalMarkerGraphVertex;
    class LocalMarkerGraphEdge;
    class LocalMarkerGraph;
    using LocalMarkerGraphBaseClass = boost::adjacency_list<
        boost::listS,   // Permit parallel edges created by createMarkerGraphEdgesStrict
        boost::listS,
        boost::bidirectionalS,
        LocalMarkerGraphVertex,
        LocalMarkerGraphEdge
        >;

    class CompressedMarker;
    class ConsensusCaller;
    class LocalMarkerGraphRequestParameters;
    class LongBaseSequences;
}



class shasta::LocalMarkerGraphVertex {
public:

    // The global vertex id of the vertex of the global marker
    // graph that corresponds to this vertex.
    MarkerGraph::VertexId vertexId;

    // The distance from the start vertex.
    uint64_t distance;

    // The markers of this vertex.
    class MarkerInfo {
    public:
        MarkerId markerId;
        OrientedReadId orientedReadId;
        uint32_t ordinal;
    };
    vector<MarkerInfo> markerInfos;

    LocalMarkerGraphVertex(
        MarkerGraph::VertexId vertexId,
        uint64_t distance) :
        vertexId(vertexId),
        distance(distance)
        {}

    // Look for the ordinal for a given oriented read id.
    // If found, returns pair(true, ordinal).
    // Otherwise, returns pair(false, don't care).
    // If more than an ordinal is found, the first one is returned.
    pair<bool, uint32_t> getOrdinal(OrientedReadId) const;

    // Coverage information at each of the k positions.
    // All reads agree on the bases, which are the marker bases,
    // but the repeat counts can be different.
    vector<Coverage> coverages;

    vector<uint8_t> storedConsensusRepeatCounts;

    // Fields used by approximateTopologicalSort.
    uint32_t color = 0;
    size_t rank = 0;

    // Compute coverage for each strand.
    array<uint64_t, 2> strandCoverage() const
    {
        array<uint64_t, 2> c = {0, 0};
        for(const MarkerInfo& markerInfo: markerInfos) {
            const Strand strand = markerInfo.orientedReadId.getStrand();
            ++c[strand];
        }
        return c;
    }

};



class shasta::LocalMarkerGraphEdge {
public:

    // Class to describe the intervening sequence between
    // the two markers that define the edge.
    class Sequence {
    public:
        // The number of overlapping bases between the
        // marker of the source vertex and the
        // marker of the target vertex.
        uint8_t overlappingBaseCount;

        // The intervening sequence between the two markers.
        vector<Base> sequence;

        // There are three possible cases:
        // overlappingBaseCount>0 && sequence.empty():
        //     The sequence of the two markers overlap by
        //     overlappingBaseCount bases.
        // overlappingBaseCount==0 && !sequence.empty():
        //     Sequence stores the base sequence between the two markers.
        // overlappingBaseCount==0 && sequence.empty():
        //     The two markers immediately follow each other,
        //     without any intervening sequence.

        bool operator<(const Sequence& that) const
        {
            return tie(overlappingBaseCount, sequence) <
                tie(that.overlappingBaseCount, that.sequence);
        }
    };



    // The oriented vertices of this edge, grouped by sequence.
    // Sorted by decreasing number of supporting reads.
    vector< pair<Sequence, vector<MarkerIntervalWithRepeatCounts> > > infos;

    // Coverage is the total number of reads supporting this edge.
    size_t coverage() const
    {
        size_t c = 0;
        for(const auto& p: infos) {
            c += p.second.size();
        }
        return c;
    }

    // Compute coverage for each strand.
    array<uint64_t, 2> strandCoverage() const
    {
        array<uint64_t, 2> c = {0, 0};
        for(const auto& info: infos) {
            for(const MarkerIntervalWithRepeatCounts& interval: info.second) {
                const Strand strand = interval.orientedReadId.getStrand();
                ++c[strand];
            }

        }
        return c;
    }

    // Look for the ordinals for a given oriented read id.
    // If found, returns true.
    // If more than an ordinal pairs is found, the first one is returned.
    bool getOrdinals(OrientedReadId, array<uint32_t, 2>& ordinals) const;

    // Id of the global edge corresponding to this edge.
    // Only filled in when the graph is created using stored connectivity.
    MarkerGraph::EdgeId edgeId = MarkerGraph::invalidEdgeId;

    // The locations of this edge in the assembly graph.
    // Each location is a pair(AssemblyGraph::EdgeId, position),
    // where position is the index of this marker graph edge in the assembly graph edge.
    // Note that, before detangling, or if detangling is not used,
    // each marker graph edge corresponds to at most one location
    // in the assembly graph. However, after detangling a marker
    // graph edge can correspond to multiple locations in the
    // assembly graph.
    vector< pair<AssemblyGraph::EdgeId, uint32_t> > assemblyGraphLocations;

    // Flag that is set if the edge was removed during
    // approximate transitive reduction by flagWeakMarkerGraphEdges.
    uint8_t wasRemovedByTransitiveReduction = 0;

    // Set if this edge was removed during pruning.
    uint8_t wasPruned = 0;

    // Set if this edge belongs to a bubble/superbubble that was removed.
    uint8_t isSuperBubbleEdge = 0;

    // Set if this edge was removed as a low coverage cross edge
    uint8_t isLowCoverageCrossEdge = 0;

    // Flag that is set if this edge corresponds to a global marker graph
    // edge that was assembled.
    uint8_t wasAssembled = 0;

    // Flag for secondary edges in assembly mode 1.
    uint8_t isSecondary = 0;

    // Field used by approximateTopologicalSort.
    bool isDagEdge = true;

    // Consensus sequence for this edge.
    vector<Base> consensusSequence;
    vector<uint32_t> consensusRepeatCounts;
    uint8_t consensusOverlappingBaseCount = 0;

};



class shasta::LocalMarkerGraph :
    public LocalMarkerGraphBaseClass {
public:

    LocalMarkerGraph(
        uint64_t readRepresentation,
        uint32_t k,
        const Reads& reads,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MemoryMapped::Vector<MarkerGraph::CompressedVertexId>& globalMarkerGraphVertex,
        const ConsensusCaller&
        );

    // Override base class Base defined in Boost Graph library.
    // Use shasta::Base instead.
    using Base = shasta::Base;

    // Find out if a vertex with the given MarkerGraph::VertexId exists.
    // If it exists, return make_pair(true, v).
    // Otherwise, return make_pair(false, null_vertex());
    pair<bool, vertex_descriptor> findVertex(MarkerGraph::VertexId) const;

    // Add a vertex with the given MarkerGraph::VertexId
    // and return its vertex descriptor.
    // A vertex with this MarkerGraph::VertexId must not exist.
    vertex_descriptor addVertex(
        MarkerGraph::VertexId,
        uint64_t distance,
        span<MarkerId> markers);

    // Get the KmerId for a vertex.
    KmerId getKmerId(vertex_descriptor) const;

    // Get the repeat counts for a MarkerInfo of a vertex.
    vector<uint8_t> getRepeatCounts(const LocalMarkerGraphVertex::MarkerInfo&) const;

    // Fill in the ConsensusInfo's for each vertex.
    void computeVertexConsensusInfo();
    void computeVertexConsensusInfo(vertex_descriptor);

    // Store sequence information in the edge.
    // Takes as input a vector of the
    // LocalMarkerGraphEdge::Info that caused the edge to be created.
    void storeEdgeInfo(edge_descriptor, const vector<MarkerInterval>&);


    // Write in Graphviz format.
    void write(
        ostream&,
        const LocalMarkerGraphRequestParameters&) const;
    void write(
        const string& fileName,
        const LocalMarkerGraphRequestParameters&) const;


    // Approximate topological sort, adding edges
    // in order of decreasing coverage. The topological sort rank
    // of each vertex is stored in LocalMarkerGrapg2Vertex::rank.
    // In addition, the vertices are stored in topological sort order
    // in vector topologicallySortedVertices.
    void approximateTopologicalSort();
    vector<vertex_descriptor> topologicallySortedVertices;

private:

    // Map a global vertex id to a vertex descriptor for the local graph.
    std::map<MarkerGraph::VertexId, vertex_descriptor> vertexMap;

    // Reads representation: 0 = raw, 1 = RLE.
    uint64_t readRepresentation;

    // The length of k-mers used as markers.
    uint32_t k;

    // Reference to the global data structure containing all reads and markers
    // (not just those in this local marker graph).
    const Reads& reads;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;

    // A reference to the vector containing the global marker graph vertex id
    // corresponding to each marker.
    // Indexed by MarkerId.
    const MemoryMapped::Vector<MarkerGraph::CompressedVertexId>& globalMarkerGraphVertex;

    // Object used to compute consensus bases and repeat counts.
    // This is owned by the caller (the Assembler object).
    const ConsensusCaller& consensusCaller;


    // Class used for graphviz output.
    class Writer : public LocalMarkerGraphRequestParameters {
    public:
        Writer(
            const LocalMarkerGraph&,
            const LocalMarkerGraphRequestParameters&);

        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const LocalMarkerGraph& graph;

        // Vertex and edge colors.
        static const string vertexColorZeroDistance;
        static const string vertexColorIntermediateDistance;
        static const string vertexColorMaxDistance;
        static const string edgeArrowColorRemovedDuringTransitiveReduction;
        static const string edgeArrowColorRemovedDuringPruning;
        static const string edgeArrowColorRemovedDuringSuperBubbleRemoval;
        static const string edgeArrowColorRemovedAsLowCoverageCrossEdge;
        static const string edgeArrowColorNotRemovedNotAssembled;
        static const string edgeArrowColorNotRemovedAssembled;
        static const string edgeLabelColorRemovedDuringTransitiveReduction;
        static const string edgeLabelColorRemovedDuringPruning;
        static const string edgeLabelColorRemovedDuringSuperBubbleRemoval;
        static const string edgeLabelColorRemovedAsLowCoverageCrossEdge;
        static const string edgeLabelColorNotRemovedNotAssembled;
        static const string edgeLabelColorNotRemovedAssembled;
        string vertexColor(const LocalMarkerGraphVertex&) const;
        string edgeArrowColor(const LocalMarkerGraphEdge&) const;
        string edgeLabelColor(const LocalMarkerGraphEdge&) const;
    };
    friend class Writer;

public:
    static void writeColorLegendVerticesByDistance(ostream&);
    static void writeColorLegendEdgeArrowsByFlags(ostream&);
    static void writeColorLegendEdgeLabelsByFlags(ostream&);
};

#endif

