#ifndef SHASTA_DIRECTED_READ_GRAPH_HPP
#define SHASTA_DIRECTED_READ_GRAPH_HPP


/*******************************************************************************

Directed version of the read graph.

Each vertex correspond to an oriented reads (therefore each read corresponds
to two vertices, one for each orientation).

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
#include "MemoryMappedDirectedGraph.hpp"

namespace shasta {
    class DirectedReadGraph;
    class DirectedReadGraphEdge;
    class DirectedReadGraphVertex;
}



// A vertex of the directed read graph.
// Even though the vertex is empty, is stil occupies one bte.
class shasta::DirectedReadGraphVertex {

};



// An edge of the directed read graph.
class shasta::DirectedReadGraphEdge {

};



class shasta::DirectedReadGraph :
    public MemoryMapped::DirectedGraph<DirectedReadGraphVertex, DirectedReadGraphEdge> {

};

#endif

