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
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

// Standard library.
#include "cstdint.hpp"
#include <limits>
#include <random>


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
    uint64_t alignmentId : 62;

    // Flag set for an edge that jumps across strands.
    uint64_t crossesStrands : 1;

    // This is set to indicate that the alignment corresponding to this edge
    // was flagged as inconsistent.
    uint64_t hasInconsistentAlignment : 1;

    ReadGraphEdge()
    {
        crossesStrands = 0;
        hasInconsistentAlignment  = 0;
    }

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

    // The edges are stored with reverse complemented pairs at
    // consecutive positions. That way, to get the reverse complement of
    // an edge id we just reverse its lowest significant bit.
    MemoryMapped::Vector<ReadGraphEdge> edges;
    uint64_t getReverseComplementEdgeId(uint64_t edgeId) const;

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
    void remove();
    
    static const uint32_t infiniteDistance;

    void findNeighbors(OrientedReadId, vector<OrientedReadId>&) const;
    void findNeighbors(OrientedReadId, uint64_t maxDistance, vector<OrientedReadId>&) const;

    // Find "bridges" from the read graph.
    // Takes as input a vector<bool> that says, for each alignmentId,
    // whether that alignment is used in the read graph.
    // Updates that vector to set to false the entries corresponding
    // to read graph "bridges".
    void findBridges(vector<bool>& keepAlignment, uint64_t maxDistance);

    void clustering(
        std::mt19937& randomSource,
        vector<ReadId>& cluster,
        bool debug) const;
};



#endif
