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

    using DirectedReadGraphBaseClass =
        MemoryMapped::DirectedGraph<DirectedReadGraphVertex, DirectedReadGraphEdge>;
}



// A vertex of the directed read graph.
class shasta::DirectedReadGraphVertex {
public:
    DirectedReadGraphBaseClass::VertexId reverseComplementedVertexId =
        DirectedReadGraphBaseClass::invalidVertexId;
};



// An edge of the directed read graph.
class shasta::DirectedReadGraphEdge {
public:
    AlignmentInfo alignmentInfo;
    DirectedReadGraphBaseClass::EdgeId reverseComplementedEdgeId =
        DirectedReadGraphBaseClass::invalidEdgeId;;
    DirectedReadGraphEdge(const AlignmentInfo& alignmentInfo) :
        alignmentInfo(alignmentInfo) {}
    DirectedReadGraphEdge() {}

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

private:

    // Add an edge 0->1, reversing the direction if necessary
    EdgeId addEdge(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        AlignmentInfo);

};

#endif

