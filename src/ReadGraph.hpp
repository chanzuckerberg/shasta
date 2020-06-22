#ifndef SHASTA_READ_GRAPH_HPP
#define SHASTA_READ_GRAPH_HPP



/*******************************************************************************

The read graph is an undirected graph in which each vertex represents
an oriented read. It is similar to the alignment graph,
in which an undirected edge is created if we found an alignment
between the corresponding oriented reads. However,
the read graph only uses a subset of the alignments.

Class ReadGraph is used to store the ReadGraph in permanent
but read-only form using MemoryMapped data structures.

*******************************************************************************/

// Shasta.
#include "MemoryMappedDirectedGraph.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

// Standard library.
#include "cstdint.hpp"
#include <limits>

namespace shasta {
    class ReadGraph;
    class ReadGraphEdge;
}



// An edge of the read graph.
class shasta::ReadGraphEdge {
public:
    array<OrientedReadId, 2> orientedReadIds;

    // The id of the alignment that corresponds to the edge.
    // Note  that if an alignment is used to generate an edge,
    // it is also used to generate the corresponding reverse complemented edge.
    uint64_t alignmentId : 63;

    // Flag set for an edge that jumps across strands.
    uint64_t crossesStrands : 1;

    // Given one of the oriented read ids, get the other.
    OrientedReadId getOther(OrientedReadId orientedReadId) const
    {
        if(orientedReadId == orientedReadIds[0]) {
            return orientedReadIds[1];
        } else if(orientedReadId == orientedReadIds[1]) {
            return orientedReadIds[0];
        } else {
            // The OrientedReadId that was passed in as an argument
            // is neithher of the two OrientedReadId's of this edge.
            SHASTA_ASSERT(0);
        }
    }
};



// Class ReadGraph is used to store the ReadGraph in permanent
// but read-only form using MemoryMapped data structures.
class shasta::ReadGraph {
public:

    MemoryMapped::Vector<ReadGraphEdge> edges;

    // Connectivity of the read graph.
    // Stores, for each OrientedReadId, indexes into the edges vector
    // of the edges that this OrientedReadId is involved in.
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> connectivity;

    // Compute a shortest path, disregarding edges flagged as cross-strand edges.
    void computeShortPath(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        size_t maxDistance,

        // Edge ids of the shortest path starting at orientedReadId0 and
        // ending at orientedReadId1.
        vector<uint32_t>& path,

        // Work areas.
        vector<uint32_t>& distance, // One per vertex, equals infiniteDistance before and after.
        vector<OrientedReadId>& reachedVertices,   // For which distance is not infiniteDistance.
        vector<uint32_t>& parentEdges  // One per vertex

    );

    void unreserve();
    
    static const uint32_t infiniteDistance;
};



#endif
