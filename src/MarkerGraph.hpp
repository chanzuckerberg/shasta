#ifndef SHASTA_MARKER_GRAPH_HPP
#define SHASTA_MARKER_GRAPH_HPP

#include "Base.hpp"
#include "Coverage.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "Uint.hpp"
#include "cstdint.hpp"

namespace shasta {

    class MarkerGraph;

    // Type used to globally identify a marker on an oriented read.
    // This is the global index of the marker in Assembler::markers.
    // For a human assembly with coverage 40X the total number
    // of markers is more than 20 billions (counting both strands),
    // so this needs to be uint64_t. There could, however, be situations
    // where uint32_t is sufficient.
    using MarkerId = uint64_t;

}



class shasta::MarkerGraph : public MultithreadedObject<MarkerGraph> {
public:

    using VertexId = MarkerId;
    using EdgeId = MarkerId;
    static const VertexId invalidVertexId;
    static const EdgeId invalidEdgeId;

    // To save memory, store vertex ids using 5 bytes.
    // This allows for up to 2^40 = 1 Ti markers (both strands).
    // A human size run with 40x coverage and 10% markers
    // has around 25 G markers (both strands).
    using CompressedVertexId = Uint40;
    static const CompressedVertexId invalidCompressedVertexId;

    MarkerGraph();

    // The marker ids of the markers corresponding to
    // each vertex of the global marker graph.
    // Indexed by VertexId.
    // For a given vertex, the marker ids are sorted.
    MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId> vertices;

    // The global marker graph vertex corresponding to each marker.
    // Indexed by MarkerId.
    // For markers that don't correspond to a marker graph vertex,
    // this stores invalidCompressedVertexId.
    MemoryMapped::Vector<CompressedVertexId> vertexTable;



    // Remove marker graph vertices and update vertices and vertexTable.
    // After this is called, the only
    // two MarkerGraph field filled in are vertices and vertexTable.
    // Everything else has to be recreated.
    void removeVertices(
        const MemoryMapped::Vector<VertexId>& verticesToBeKept,
        uint64_t pageSize,
        uint64_t threadCount);
private:
    class RemoveVerticesData {
    public:
        const MemoryMapped::Vector<VertexId>* verticesToBeKept;
        MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId> newVertices;
    };
    RemoveVerticesData removeVerticesData;
    void removeVerticesThreadFunction1(size_t threadId);
    void removeVerticesThreadFunction2(size_t threadId);
    void removeVerticesThreadFunction3(size_t threadId);
public:



    // The reverse complement of each vertex.
    // Indexed by VertexId.
    MemoryMapped::Vector<VertexId> reverseComplementVertex;

    // The edges of the marker graph.
    class Edge {
    public:
        Uint40 source;  // The source vertex (index into globalMarkerGraphVertices).
        Uint40 target;  // The target vertex (index into globalMarkerGraphVertices).
        uint8_t coverage;   // (255 indicates 255 or more).

        // Flags used to mark the edge as removed from the marker graph.
        bool wasRemoved() const
        {
            return
                wasRemovedByTransitiveReduction ||
                wasPruned ||
                isLowCoverageCrossEdge ||
                isSuperBubbleEdge;
        }

        // Flag that is set if the edge was removed during
        // approximate transitive reduction by flagWeakMarkerGraphEdges.
        uint8_t wasRemovedByTransitiveReduction : 1;

        // Set if this edge was removed during pruning.
        uint8_t wasPruned : 1;

        // Set if this edge belongs to a bubble/superbubble that was removed.
        uint8_t isSuperBubbleEdge : 1;

        // Flag set if this edge corresponds to a low coverage cross edge
        // of the assembly graph.
        uint8_t isLowCoverageCrossEdge: 1;

        // Flag set if this edge was assembled.
        // If set, edgeConsensusOverlappingBaseCount and edgeConsensus
        // for this edge are set.
        uint8_t wasAssembled : 1;

        // Unused.
        uint8_t flag4 : 1;
        uint8_t flag5 : 1;
        uint8_t flag6 : 1;

        void clearFlags()
        {
            wasRemovedByTransitiveReduction = 0;
            wasPruned = 0;
            isSuperBubbleEdge = 0;
            isLowCoverageCrossEdge = 0;
            wasAssembled = 0;
            flag4 = 0;
            flag5 = 0;
            flag6 = 0;
        }
        Edge() :
            source(MarkerGraph::invalidCompressedVertexId),
            target(MarkerGraph::invalidCompressedVertexId),
            coverage(0)
        {
            clearFlags();
        }
    };
    MemoryMapped::Vector<Edge> edges;
    const Edge* findEdge(Uint40 source, Uint40 target) const;
    EdgeId findEdgeId(Uint40 source, Uint40 target) const;

    // The MarkerIntervals for each of the above edges.
    MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> edgeMarkerIntervals;

    // The edges that each vertex is the source of.
    // Contains indexes into the above edges vector.
    MemoryMapped::VectorOfVectors<Uint40, uint64_t> edgesBySource;

    // The edges that each vertex is the target of.
    // Contains indexes into the above edges vector.
    MemoryMapped::VectorOfVectors<Uint40, uint64_t> edgesByTarget;

    // Compute in-degree or out-degree of a vertex,
    // counting only edges that were not removed.
    uint64_t inDegree(VertexId) const;
    uint64_t outDegree(VertexId) const;

    // The reverse complement of each edge.
    // Indexed by EdgeId.
    MemoryMapped::Vector<EdgeId> reverseComplementEdge;

    // The consensus repeat counts of each vertex of the marker graph.
    // There are assemblerInfo->k entries for each vertex.
    // The first entry for a vertex is at index vertexId*assemblerInfo->k.
    MemoryMapped::Vector<uint8_t> vertexRepeatCounts;

    // Consensus sequence and repeat counts for each marker graph edge.
    // This excludes the sequence of flanking markers and their repeat counts.
    // Indexed by the marker graph edge id.
    // - For edges that were marked as removed,
    //   edgeConsensusOverlappingBaseCount is 0 and edgeConsensus is empty.
    // - For edges that were not marked as removed:
    //   * If the consensus sequence has one or more intervening bases
    //     between the flanking markers,
    //     edgeConsensusOverlappingBaseCount is 0 and edgeConsensus
    //     stores those intervening bases with their repeat count consensus.
    //   * Otherwise, edgeConsensus is empty and
    //     edgeConsensusOverlappingBaseCount stores the number of
    //     overlapping bases (for the consensus sequence)
    //     between the two flanking markers. This can be zero
    //     if the consensus sequence has the flanking markers
    //     exactly adjacent.
    MemoryMapped::VectorOfVectors<pair<Base, uint8_t>, uint64_t> edgeConsensus;
    MemoryMapped::Vector<uint8_t> edgeConsensusOverlappingBaseCount;


    // Details of vertex coverage.
    // These are not stored by default.
    // They can be used to calibrate the Bayesian model for repeat counts
    // and for some types of analyses.
    // Indeed by VertexId. For each vertex, contains pairs (position, CompressedCoverageData),
    // ordered by position.
    // Note that the bases at a given position are all identical by construction.
    MemoryMapped::VectorOfVectors<pair<uint32_t, CompressedCoverageData>, uint64_t>
        vertexCoverageData;

    // Details of edge coverage.
    // These are not stored by default.
    // They can be used to calibrate the Bayesian model for repeat counts
    // and for some types of analyses.
    // Indeed by EdgeId. For each edge, contains pairs (position, CompressedCoverageData),
    // ordered by position.
    MemoryMapped::VectorOfVectors<pair<uint32_t, CompressedCoverageData>, uint64_t>
        edgeCoverageData;
};

#endif
