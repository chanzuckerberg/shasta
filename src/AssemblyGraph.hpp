#ifndef SHASTA_ASSEMBLY_GRAPH_HPP
#define SHASTA_ASSEMBLY_GRAPH_HPP



/***************************************************************************

The assembly graph is a "compressed" version of the marker graph,
in which each linear chain of edges is replaced by a single edge.

Each assembly graph vertex  corresponds to a marker graph vertex,
but the reverse is not true, because marker graph vertices
that are internal to a linear chain of edges don't have a corresponding
vertex in the assembly graph.

***************************************************************************/

// Shasta.
#include "LongBaseSequence.hpp"
#include "MarkerGraph.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

// Standard library.
#include <limits>

namespace shasta {
    class AssemblyGraph;
}



class shasta::AssemblyGraph {
public:

    // Use the same vertex and edge ids of the marker graph.
    // We could probably get away with 32 bits.
    using VertexId = MarkerGraph::VertexId;
    using EdgeId = MarkerGraph::EdgeId;
    static const VertexId invalidVertexId = std::numeric_limits<VertexId>::max();
    static const VertexId invalidEdgeId = std::numeric_limits<EdgeId>::max();

    // The vertices of the assembly graph.
    // Each assembly graph vertex  corresponds to a marker graph vertex.
    // Indexed by vertex id in the assembly graph.
    // Contains the corresponding vertex id in the marker graph.
    MemoryMapped::Vector<VertexId> vertices;

    // The reverse complement of each vertex.
    // Indexed by VertexId.
    MemoryMapped::Vector<VertexId> reverseComplementVertex;



    // The edges of the assembly graph.
    // Indexed by edge id in the assembly graph.
    class Edge {
    public:
        VertexId source;
        VertexId target;

        // Minimum, average, and maximum coverage of marker graph vertices and edges
        // corresponding to this assembly graph edge.
        // Only internal vertices contribute to this -
        // the first and last vertex don't contribute.
        // If there is only one marker graph edge, there are
        // no internal vertices, and in that case vertex coverage
        // metrics are set to zero.
        uint32_t minVertexCoverage;
        uint32_t averageVertexCoverage;
        uint32_t maxVertexCoverage;
        uint32_t minEdgeCoverage;
        uint32_t averageEdgeCoverage;
        uint32_t maxEdgeCoverage;

        // Reason for removal.
        enum class RemovalReason : uint8_t {
            NotRemoved = 0,
            LowCoverageCrossEdge = 1
        };
        RemovalReason removalReason;
        bool wasRemoved() const
        {
            return removalReason != RemovalReason::NotRemoved;
        }

        Edge() : removalReason(RemovalReason::NotRemoved) {}
    };
    MemoryMapped::Vector<Edge> edges;

    // Return the number of edges that were not removed.
    EdgeId edgeCount() const;

    // The reverse complement of each edge.
    // Indexed by EdgeId.
    MemoryMapped::Vector<EdgeId> reverseComplementEdge;

    // Return true if this edge is an assembled edge.
    // To avoid assembling both strands, we only assemble
    // one edge in each reverse complemented pair.
    bool isAssembledEdge(EdgeId edgeId) const
    {
        return edgeId <= reverseComplementEdge[edgeId];
    }



    // The edges that each vertex is the source of.
    // Indexed by vertex id in the assembly graph.
    // Contains edge ids that can be used as indexes into edges and edgeLists.
    MemoryMapped::VectorOfVectors<EdgeId, EdgeId> edgesBySource;

    // The edges that each vertex is the target of.
    // Indexed by vertex id in the assembly graph.
    // Contains edge ids that can be used as indexes into edges and edgeLists.
    MemoryMapped::VectorOfVectors<EdgeId, EdgeId> edgesByTarget;

    // Fill in edgesBySource and edgesByTarget.
    void computeConnectivity();

    // The edge ids of global marker graph edges corresponding
    // to each edge of the assembly graph.
    // Indexed by the edge id of the assembly graph edge.
    MemoryMapped::VectorOfVectors<EdgeId, EdgeId> edgeLists;

    // Find the out-degree or in-degree of a vertex.
    // This is not simply the same as counting edgesBySource
    // and edgesByTarget, because we have to skip edges
    // that were removed.
    VertexId inDegree(VertexId) const;
    VertexId outDegree(VertexId) const;

    // Find in-edges/out-edges of a vertex
    // that were not removed.
    // They are returned sorted by edge id.
    void findInEdges(VertexId, vector<EdgeId>&) const;
    void findOutEdges(VertexId, vector<EdgeId>&) const;



    // Bubbles in the assembly graph.
    // A bubble is a set of two vertices v0, v1
    // such that the outgoing edges of v0 are the same
    // as the outgoing edges of v1, and
    // out-degree(v0) = in-degree(v1) > 1.
    // v0 is called the bubble source.
    // v1 is called the bubble target.
    // In defining and detecting bubbles, edges
    // that were removed are considered to not exist.
    class Bubble {
    public:
        VertexId v0;
        VertexId v1;
        Bubble() {}
        Bubble(VertexId v0, VertexId v1) : v0(v0), v1(v1) {}
    };
    MemoryMapped::Vector<Bubble> bubbles;
    void findBubbles(); // Assumes bubble was already initialized.



    // Bubble chains. A bubble chain is a linear sequence of bubbles.
    // Each pair of consecutive bubbles in the sequence may be separated by
    // a homozygous segment, but the bubble can also be
    // immediately adjacent (v1 of thr first bubble is the same
    // as v0 of the second bubble).
    // For each bubble chain, we store the bubble ids (indexes in
    // the bubbles vector above).
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> bubbleChains;
    void findBubbleChains(); // Assumes bubbleChains was already initialized.


    // A table that can be used to find the locations of a marker graph
    // edge in the assembly graph, if any.
    // Note that, before detangling,or if detangling is not used,
    // each marker graph edge corresponds to at most one location
    // in the assembly graph. However, after detangling a marker
    // graph edge can correspond to multiple locations in the
    // assembly graph.
    // Indexed by the edge id in the marker graph, gives for each marker graph
    // edge a vector of pair(EdgeId, position), where:
    // - EdgeId is the id of the assembly graph edge containing the
    //   given marker graph edge.
    // - Position is the index of this marker graph edge in that
    //   assembly graph edge.
    MemoryMapped::VectorOfVectors< pair<EdgeId, uint32_t> , uint64_t> markerToAssemblyTable;
    void createMarkerToAssemblyTable(uint64_t markerGrapEdgeCount);

    // The assembled sequenced and repeat counts for each edge of the
    // assembly graph.
    // Indexed edge id in the assembly graph.
    LongBaseSequences sequences;
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t> repeatCounts;



    // The oriented reads internal to each assembly graph edge.
    // Indexed by EdgeId.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // Number of marker graph  vertices  internal to the
        // assembly graph edge and containing this oriented read.
        uint64_t vertexCount;

        // Number of marker graph  edges  internal to the
        // assembly graph edge and containing this oriented read.
        uint64_t edgeCount;

        OrientedReadInfo() {}
        OrientedReadInfo(
            OrientedReadId orientedReadId,
            uint64_t vertexCount,
            uint64_t edgeCount) :
            orientedReadId(orientedReadId),
            vertexCount(vertexCount),
            edgeCount(edgeCount) {}
    };
    MemoryMapped::VectorOfVectors<OrientedReadInfo, uint64_t> orientedReadsByEdge;

    // Compute the number of oriented reads in common between two segments.
    uint64_t commonOrientedReadCount(
        EdgeId, EdgeId,
        uint64_t minVertexCount,
        uint64_t minEdgeCount) const;




    // Close all open data.
    void close();

    // Close and remove all open data.
    void remove();

    // Basic Graphviz output of the global assembly graph.
    void writeGraphviz(const string& fileName) const;

    // Create a csv file that can be loaded in Bandage to color assembled segments
    // by similarity (number of common oriented reads) with a given assembled segment.
    void colorGfaBySimilarityToSegment(
        EdgeId,
        uint64_t minVertexCount,
        uint64_t minEdgeCount);



    // Assembly graph forks.

    // A fork is the set of outgoing edges of a vertex with out-degree>1,
    // or the set of incoming edges of a vertex with in-degree>1.
    // The edges of the set are called the fork branches.
    // A fork can have any number of branches>1.

    // A bubble generates only one fork (the set of outgoing edges
    // of the bubble source vertex which is also the set of
    // incoming edges of the bubble target vertex).
    class Fork {
    public:
        // The source or target vertex of the fork.
        VertexId vertexId;

        // Flag which is true if the fork consists of the vertex's
        // outgoing edges, false if the fork consists of
        // thevertex's incoming edges.
        // The fork corresponding to a bubble is stored only once.
        bool isForward;
    };
    MemoryMapped::Vector<Fork> forks;
    void createForks();
    span<EdgeId> getForkEdges(uint64_t forkId);

private:
    // Class used by createForks.
    class ForkInfo : public Fork {
    public:
        vector<EdgeId> edgeIds;

        // Compare using only the edgeIds.
        bool operator==(const ForkInfo& that) const
        {
            return edgeIds == that.edgeIds;
        }
        bool operator<(const ForkInfo& that) const
        {
            return edgeIds < that.edgeIds;
        }
    };



};



#endif
