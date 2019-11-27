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

    // The namber of raw (not RLE) bases and markers
    // for the oriented read corresponding to this vertex.
    uint32_t baseCount;
    uint32_t markerCount;

    // Flag set if there is one alignment in which this oriented read
    // is entirely contained in an another oriented read,
    // except possibly for up to maxTrim markers at each end.
    uint8_t isContained : 1;

    // The VertexId of the reverse complement of this vertex.
    DirectedReadGraphBaseClass::VertexId reverseComplementedVertexId =
        DirectedReadGraphBaseClass::invalidVertexId;

    DirectedReadGraphVertex() :
        baseCount(0), markerCount(0)
    {
        isContained = 0;
    }
};



// An edge of the directed read graph.
class shasta::DirectedReadGraphEdge {
public:

    // Information on the alignment that generated this edge.
    AlignmentInfo alignmentInfo;

    // The EdgeId of the reverse complement of this edge.
    DirectedReadGraphBaseClass::EdgeId reverseComplementedEdgeId =
        DirectedReadGraphBaseClass::invalidEdgeId;

    // Edge flags.
    uint8_t involvesTwoContainedVertices : 1;
    uint8_t involvesOneContainedVertex : 1;

    // Constructors.
    DirectedReadGraphEdge(const AlignmentInfo& alignmentInfo) :
        alignmentInfo(alignmentInfo)
    {
        clearFlags();
    }
    DirectedReadGraphEdge()
    {
        clearFlags();
    }

    void clearFlags()
    {
        involvesTwoContainedVertices = 0;
        involvesOneContainedVertex = 0;
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

    // Flag contained vertices and set edge flags accordingly.
    void flagContainedVertices(uint32_t maxTrim);

    // Create a LocalDirectedReadGraph.
    bool extractLocalSubgraph(
        OrientedReadId,
        uint64_t maxDistance,
        uint64_t minAlignedMarkerCount,
        uint64_t maxOffsetAtCenter,
        double minAlignedFraction,
        bool allowEdgesInvolvingTwoContainedVertices,
        bool allowEdgesInvolvingOneContainedVertex,
        double timeout,
        LocalDirectedReadGraph&);

    void writeEdges();

private:

    // Add an edge 0->1, reversing the direction if necessary
    EdgeId addEdge(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        AlignmentInfo);

    // And edge checker that allows only edges that satisfy specify criteria.
    // Used to create the local directed read graph for display.
    class EdgeFilter : public AbstractEdgeFilter {
    public:
        EdgeFilter(
            uint64_t minAlignedMarkerCount,
            uint64_t maxTwiceOffsetAtCenter,
            double minAlignedFraction,
            bool allowEdgesInvolvingTwoContainedVertices,
            bool allowEdgesInvolvingOneContainedVertex) :

            minAlignedMarkerCount(minAlignedMarkerCount),
            maxTwiceOffsetAtCenter(maxTwiceOffsetAtCenter),
            minAlignedFraction(minAlignedFraction),
            allowEdgesInvolvingTwoContainedVertices(allowEdgesInvolvingTwoContainedVertices),
            allowEdgesInvolvingOneContainedVertex(allowEdgesInvolvingOneContainedVertex)
            {}

        bool allowEdge(EdgeId edgeId, const Edge& edge) const
        {
            if(not allowEdgesInvolvingTwoContainedVertices and edge.involvesTwoContainedVertices) {
                return false;
            }
            if(not allowEdgesInvolvingOneContainedVertex and edge.involvesOneContainedVertex) {
                return false;
            }
            return
                edge.alignmentInfo.markerCount >= minAlignedMarkerCount
                and
                abs(edge.alignmentInfo.twiceOffsetAtCenter()) <= maxTwiceOffsetAtCenter
                and
                edge.alignmentInfo.minAlignedFraction() >= minAlignedFraction
                ;
        }

        uint64_t minAlignedMarkerCount;
        uint64_t maxTwiceOffsetAtCenter;
        double minAlignedFraction;

        bool allowEdgesInvolvingTwoContainedVertices;
        bool allowEdgesInvolvingOneContainedVertex;
    };


};

#endif

