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

A vertex is flagged as contained if there is at least one alignment in which
the oriented read corresponding to the vertex
is entirely contained in an another oriented read,
except possibly for up to maxTrim markers at each end.

The read graph is initially created by adding two edges for each
known alignment. Then, a subset of all the edges are flagged
as "keep" as follows:

- For a contained vertex, the best containedNeighborCount adjacent edges,
as defined by number of aligned markers, are marked as "keep".

- For an uncontained vertex, the best uncontainedNeighborCountPerDirection
out-edges and the best uncontainedNeighborCountPerDirection in-edges
of each vertex are marked as "keep", considering only out-edges and in-edges
to other uncontained vertices.

Only edges marked as "keep" are used to create the marker graph.

*******************************************************************************/

// Shasta.
#include "Alignment.hpp"
#include "MemoryMappedDirectedGraph.hpp"
#include "ReadId.hpp"

namespace shasta {
    class Assembler;
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

    // The number of raw (not RLE) bases and markers
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
    uint8_t keep : 1;

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
        keep = 0;
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

    // Flag as "keep" a subset of all edges.
    // These are the edges that will be used to create the marker graph.
    // See comments at the beginning of this file for more information.
    void flagEdgesToBeKept(
        uint64_t containedNeighborCount,
        uint64_t uncontainedNeighborCountPerDirection);

    // Create a LocalDirectedReadGraph.
    bool extractLocalSubgraph(
        OrientedReadId,
        uint64_t maxDistance,
        uint64_t minAlignedMarkerCount,
        uint64_t maxOffsetAtCenter,
        double minAlignedFraction,
        bool allowEdgesNotKept,
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
            bool allowEdgesNotKept) :

            minAlignedMarkerCount(minAlignedMarkerCount),
            maxTwiceOffsetAtCenter(maxTwiceOffsetAtCenter),
            minAlignedFraction(minAlignedFraction),
            allowEdgesNotKept(allowEdgesNotKept)
            {}

        bool allowEdge(EdgeId edgeId, const Edge& edge) const
        {
            if(not allowEdgesNotKept and not edge.keep) {
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

        bool allowEdgesNotKept;
    };


};

#endif

