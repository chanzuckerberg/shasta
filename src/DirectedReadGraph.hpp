#ifndef SHASTA_DIRECTED_READ_GRAPH_HPP
#define SHASTA_DIRECTED_READ_GRAPH_HPP


/*******************************************************************************

Directed version of the read graph.

Each vertex correspond to an oriented reads (therefore each read corresponds
to two vertices, one for each orientation). The vertex id of
the read corresponding to each oriented read is orientedRead.getValue().

Edges are directed so that based on the offset between centers of the two reads.
More precisely, consider vertex v0 corresponding to oriented read r0 and
vertex v1 corresponding to oriented read r1.
A directed edge v0->v1 is created if there is a marker alignment with r0 and r1
that places the center of r1 to the right of the center of r0.
In other words, the offset between the centers of r0 and r1 must be non-negative.

To enforce this, when adding edges we compute the offset between centers,
and reverse the edge if necessary.

If the offset between centers is zero, we break the tie in a way that leaves
the read graph invariant under reverse complementing, as follows:
- If strand0 == strand1 ==0, read0 must be less than read1.
- If strand0 == strand1 ==1, read0 must be greater than read1.
- If strand0 != strand1, read0 must be on strand 0 and read1 must be on strand1.

*******************************************************************************/

// Shasta.
#include "Alignment.hpp"
#include "MemoryMappedDirectedGraph.hpp"
#include "ReadId.hpp"

namespace shasta {
    class DirectedReadGraph;
    class DirectedReadGraphEdge;
    class DirectedReadGraphVertex;

    class LocalDirectedReadGraph;

    using DirectedReadGraphBaseClass =
        MemoryMapped::DirectedGraph<DirectedReadGraphVertex, DirectedReadGraphEdge>;
}



// A vertex of the directed read graph.
class shasta::DirectedReadGraphVertex {
public:
    // The number of raw bases and the number of markers
    // of the oriented read corresponding to this vertex.
    uint64_t baseCount;
    uint64_t markerCount;

    DirectedReadGraphBaseClass::VertexId reverseComplementedVertexId =
        DirectedReadGraphBaseClass::invalidVertexId;
};



// An edge of the directed read graph.
class shasta::DirectedReadGraphEdge {
public:

    // Information on the alignment that generated this edge.
    AlignmentInfo alignmentInfo;

    // The reverse complement of this edge.
    DirectedReadGraphBaseClass::EdgeId reverseComplementedEdgeId =
        DirectedReadGraphBaseClass::invalidEdgeId;

    // Flag set if this edge is removed due to transitive reduction.
    uint8_t wasRemovedByTransitiveReduction:1;

    // Transitive coverage begins at 1 and gets set to zero for edges
    // removed during transitive reduction.
    // It is incremented on all edges of the path that causes an edge to be removed.
    uint32_t transitiveCoverage;

    // Constructors.
    DirectedReadGraphEdge(const AlignmentInfo& alignmentInfo) :
        alignmentInfo(alignmentInfo), transitiveCoverage(1)
    {
        wasRemovedByTransitiveReduction = 0;
    }
    DirectedReadGraphEdge() :
        transitiveCoverage(1)
    {
        wasRemovedByTransitiveReduction = 0;
    }

};



class shasta::DirectedReadGraph :
    public DirectedReadGraphBaseClass {
public:
    using BaseClass = DirectedReadGraphBaseClass;
    using Vertex = DirectedReadGraphVertex;
    using Edge = DirectedReadGraphEdge;

    void createVertices(ReadId readCount);

    // Add a pair of edges corresponding to an alignment.
    void addEdgePair(const AlignmentData&);

    // Make sure the graph is invariant under reverse complementing.
    void check();

    // Create a LocalDirectedReadGraph.
    bool extractLocalSubgraph(
        OrientedReadId,
        uint64_t maxDistance,
        uint64_t minAlignedMarkerCount,
        double minAlignedFraction,
        bool allowTransitiveReductionEdges,
        double timeout,
        LocalDirectedReadGraph&);

    void transitiveReduction(
        double offsetTolerance0,
        double offsetTolerance1);

private:

    // Add an edge 0->1, reversing the direction if necessary
    EdgeId addEdge(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        AlignmentInfo);

    // And edge checker that allows only edges not removed by transitive reduction.
    class EdgeFilter : public AbstractEdgeFilter {
    public:
        EdgeFilter(
            uint64_t minAlignedMarkerCount,
            double minAlignedFraction,
            bool allowTransitiveReductionEdges) :
            minAlignedMarkerCount(minAlignedMarkerCount),
            minAlignedFraction(minAlignedFraction),
            allowTransitiveReductionEdges(allowTransitiveReductionEdges) {}

        bool allowEdge(EdgeId edgeId, const Edge& edge) const
        {
            if(not allowTransitiveReductionEdges and edge.wasRemovedByTransitiveReduction) {
                return false;
            }
            return
                edge.alignmentInfo.markerCount >= minAlignedMarkerCount
                and
                min(edge.alignmentInfo.alignedFraction(0), edge.alignmentInfo.alignedFraction(1))
                    >= minAlignedFraction;
        }
        uint64_t minAlignedMarkerCount;
        double minAlignedFraction;
        bool allowTransitiveReductionEdges;
    };


};

#endif

