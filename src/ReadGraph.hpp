#ifndef CZI_SHASTA_READ_GRAPH_HPP
#define CZI_SHASTA_READ_GRAPH_HPP



/*******************************************************************************

The read graph is an undirected graph in which each vertex represents
an oriented read. It is similar to the alignment graph,
in which an undirected edge is created if we found an alignment
between the corresponding oriented reads. However,
the read graph only uses a subset of the alignments.

Currently, we keep the best maxAlignmentCount alignments
for each oriented read. Therefore, the read graph is
a k-Nearest-Neighbor (k-NN) version of the alignment graph.

*******************************************************************************/

// Shasta.
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

// Standard library.
#include "cstdint.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class ReadGraph;
    }
}



class ChanZuckerberg::shasta::ReadGraph {
public:

    // An edge of the read graph.
    class Edge {
    public:
        array<OrientedReadId, 2> orientedReadIds;

        // The id of the alignment that corresponds to the edge.
        // Note  that if an alignment is used to generate an edge,
        // it is also used to generate the corresponding reverse complemented edge.
        uint64_t alignmentId;

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
                CZI_ASSERT(0);
            }
        }
    };
    MemoryMapped::Vector<Edge> edges;

    // Connectivity of the read graph.
    // Stores, for each OrientedReadId, indexes into the edges vector
    // of the edges that this OrientedReadId is involved in.
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> connectivity;
};

#endif
