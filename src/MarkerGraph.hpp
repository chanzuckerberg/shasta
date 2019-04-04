#ifndef CZI_SHASTA_MARKER_GRAPH_HPP
#define CZI_SHASTA_MARKER_GRAPH_HPP

#include "MemoryMappedVectorOfVectors.hpp"
#include "Uint.hpp"
#include "cstdint.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class MarkerGraph;
    }
}



class ChanZuckerberg::shasta::MarkerGraph {
public:

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
                isSuperBubbleEdge;
        }

        // Flag that is set if the edge was removed during
        // approximate transitive reduction by flagWeakMarkerGraphEdges.
        uint8_t wasRemovedByTransitiveReduction : 1;

        // Set if this edge was removed during pruning.
        uint8_t wasPruned : 1;

        // Set if this edge belongs to a bubble/superbubble that was removed.
        uint8_t isSuperBubbleEdge : 1;

        // Unused.
        uint8_t flag3 : 1;
        uint8_t flag4 : 1;
        uint8_t flag5 : 1;
        uint8_t flag6 : 1;
        uint8_t flag7 : 1;

        void clearFlags()
        {
            wasRemovedByTransitiveReduction = 0;
            wasPruned = 0;
            isSuperBubbleEdge = 0;
            flag3 = 0;
            flag4 = 0;
            flag5 = 0;
            flag6 = 0;
            flag7 = 0;
        }
        Edge() :
            source(invalidCompressedGlobalMarkerGraphVertexId),
            target(invalidCompressedGlobalMarkerGraphVertexId),
            coverage(0)
        {
            clearFlags();
        }
    };
    MemoryMapped::Vector<Edge> edges;
    const Edge* findEdge(Uint40 source, Uint40 target) const;

    // The MarkerIntervals for each of the above edges.
    MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> edgeMarkerIntervals;

    // The edges that each vertex is the source of.
    // Contains indexes into the above edges vector.
    MemoryMapped::VectorOfVectors<Uint40, uint64_t> edgesBySource;

    // The edges that each vertex is the target of.
    // Contains indexes into the above edges vector.
    MemoryMapped::VectorOfVectors<Uint40, uint64_t> edgesByTarget;

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
    //     if the consensus sequence has tghe flanking markers
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
