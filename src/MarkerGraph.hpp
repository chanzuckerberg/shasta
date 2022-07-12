#ifndef SHASTA_MARKER_GRAPH_HPP
#define SHASTA_MARKER_GRAPH_HPP

#include "Base.hpp"
#include "MarkerInterval.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "shastaTypes.hpp"
#include "Uint.hpp"

#include "cstdint.hpp"
#include "memory.hpp"

namespace shasta {

    class MarkerGraph;
    class CompressedCoverageData;

    extern template class MultithreadedObject<MarkerGraph>;
}



class shasta::MarkerGraph : public MultithreadedObject<MarkerGraph> {
public:

    using VertexId = MarkerGraphVertexId;
    using EdgeId = MarkerGraphEdgeId;
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
    // Stored as a shared pointer to permit easy replacement of the vertices.
    shared_ptr< MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId> > verticesPointer;
    void constructVertices()
    {
        verticesPointer = make_shared<MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId> >();
    }
    void destructVertices() {
        verticesPointer = 0;
    }



    // Vertices access functions.
    // Return the number of vertices.
    MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId>& vertices()
    {
        return *verticesPointer;
    }
    const MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId>& vertices() const
    {
        return *verticesPointer;
    }
    uint64_t vertexCount() const {
        return verticesPointer->size();
    }
    // Return the number of markers for a given vertex.
    uint64_t vertexCoverage(VertexId vertexId) const
    {
        return verticesPointer->size(vertexId);
    }
    // Return the marker ids for a given vertex.
    span<MarkerId> getVertexMarkerIds(VertexId vertexId) {
        return vertices()[vertexId];
    }
    span<const MarkerId> getVertexMarkerIds(VertexId vertexId) const {
        return vertices()[vertexId];
    }

    void remove();

    // The global marker graph vertex corresponding to each marker.
    // Indexed by MarkerId.
    // For markers that don't correspond to a marker graph vertex,
    // this stores invalidCompressedVertexId.
    MemoryMapped::Vector<CompressedVertexId> vertexTable;



    // This renumbers the vertex table to make sure that
    // vertices are numbered contiguously starting at 0.
    // This must be called after the vertexTable is changed,
    // as in Assembler::cleanupDuplicateMarkers.
    // After this is called, all other data structures
    // are inconsistent and need to be recreated.
    // The second version can be called if the maximum vertex id
    // present in the vertex table is already known, and is faster.
    // Returns the maximmum vertex id after renumbering.
    VertexId renumberVertexTable(size_t threadCount);
    VertexId renumberVertexTable(size_t threadCount, VertexId maxVertexId);
private:
    void renumberVertexTableThreadFunction1(size_t threadId);
    void renumberVertexTableThreadFunction2(size_t threadId);
    class RenumberVertexTableData {
    public:
        // Set to true for VertexId values represented in the starting vertexTable.
        MemoryMapped::Vector<bool> isPresent;

        // The new VertexId corresponding to each old VertexId.
        MemoryMapped::Vector<VertexId> newVertexId;
    };
    RenumberVertexTableData renumberVertexTableData;



    // Find the maximum valid VertexId in the vertex table.
    VertexId findMaxVertexTableEntry(size_t threadCount);
    void findMaxVertexTableEntryThreadFunction(size_t threadId);
    class FindMaxVertexTableEntryData {
    public:
        // The maximum VertexId found by each thread.
        vector<VertexId> threadMaxVertexId;
    };
    FindMaxVertexTableEntryData findMaxVertexTableEntryData;
public:


    // Recreate the vertices from the vertexTable.
    // This assumes that valid VertexId's in the vertex table
    // are numbered contiguously starting at 0 (call renumberVertexTable to ensure that).
    void createVerticesFromVertexTable(size_t threadCount, VertexId maxVertexId);
private:
    void createVerticesFromVertexTableThreadFunction1(size_t threadId);
    void createVerticesFromVertexTableThreadFunction2(size_t threadId);
    void createVerticesFromVertexTableThreadFunction3(size_t threadId);
    void createVerticesFromVertexTableThreadFunction4(size_t threadId);
    class CreateVerticesFromVertexTableData {
    public:
        // Like the vertices, but the second template argument is VertexId
        // instead of CompressedVertexId. This is necessariy to be able to
        // work on it in multilthreaded code efficiently.
        MemoryMapped::VectorOfVectors<MarkerId, VertexId> vertices;
    };
    CreateVerticesFromVertexTableData createVerticesFromVertexTableData;
public:



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
        shared_ptr<MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId> > newVerticesPointer;
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
                isSuperBubbleEdge ||
                wasRemovedWhileSplittingSecondaryEdges
                ;
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

        // Flag for secondary edges in assembly mode 1.
        uint8_t isSecondary;

        // This is set for secondary edges that are created and later split.
        // Assembly mode 2 only.
        uint8_t wasRemovedWhileSplittingSecondaryEdges : 1;

        // Unused.
        uint8_t flag6 : 1;

        void clearFlags()
        {
            wasRemovedByTransitiveReduction = 0;
            wasPruned = 0;
            isSuperBubbleEdge = 0;
            isLowCoverageCrossEdge = 0;
            wasAssembled = 0;
            isSecondary = 0;
            wasRemovedWhileSplittingSecondaryEdges = 0;
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
    EdgeId getFirstNonRemovedOutEdge(VertexId) const;
    EdgeId getFirstNonRemovedInEdge(VertexId) const;

    // The reverse complement of each edge.
    // Indexed by EdgeId.
    MemoryMapped::Vector<EdgeId> reverseComplementEdge;

    // Return total coverage of an edge.
    uint64_t edgeCoverage(EdgeId edgeId) const
    {
        return edgeMarkerIntervals.size(edgeId);
    }

    // Return coverage for each strand for an edge.
    array<uint64_t, 2> edgeStrandCoverage(EdgeId edgeId) const
    {
        array<uint64_t, 2> coverage = {0, 0};
        for(const MarkerInterval& markerInterval: edgeMarkerIntervals[edgeId]) {
            ++coverage[markerInterval.orientedReadId.getStrand()];
        }
        return coverage;
    }

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
