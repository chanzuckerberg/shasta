#ifndef CZI_SHASTA_READ_GRAPH_HPP
#define CZI_SHASTA_READ_GRAPH_HPP



/*******************************************************************************

The read graph is an undirected graph in which each vertex represents
an oriented read. It is similar to the alignment graph,
in which an undirected edge is created if we found an alignment
between the corresponding oriented reads. However,
the read graph only uses a subset of the alignments.

Class DynamicReadGraph is used for construction of the read graph.
It uses the Boost graph library to represent and manipulate the
read graph or a subset of it.

Class ReadGraph is used to store the ReadGraph in permanent
but read-only form using MemoryMapped data structures.

*******************************************************************************/

// Shasta.
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "cstdint.hpp"
#include <unordered_map>

namespace ChanZuckerberg {
    namespace shasta {
        class ReadGraph;
        class DynamicReadGraph;
        class DynamicReadGraphEdge;
        class DynamicReadGraphVertex;
        using DynamicReadGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::undirectedS,
            DynamicReadGraphVertex,
            DynamicReadGraphEdge>;
    }
}


class ChanZuckerberg::shasta::DynamicReadGraphVertex {
public:

    // We store the OrientedReadId as an integer
    // so we can also use it as a vertex id for Graphviz output.
    OrientedReadId::Int orientedReadIdInt;
    OrientedReadId getOrientedReadId() const
    {
        return OrientedReadId(orientedReadIdInt);
    }

    DynamicReadGraphVertex(OrientedReadId orientedReadId) :
        orientedReadIdInt(orientedReadId.getValue()) {}
};



class ChanZuckerberg::shasta::DynamicReadGraphEdge {
public:
    // The id of the alignment that corresponds to the edge.
    uint64_t alignmentId;

    DynamicReadGraphEdge(uint64_t alignmentId) :
        alignmentId(alignmentId) {}
};



// Class DynamicReadGraph is used for construction of the read graph.
// It uses the Boost graph library to represent and manipulate the
// read graph or a subset of it.
class ChanZuckerberg::shasta::DynamicReadGraph : public DynamicReadGraphBaseClass {
public:

    // Create a vertex for each of the two oriented reads
    // corresponding to readCount reads.
    // Used to create a DynamicReadGraph representing
    // the entire global read graph.
    DynamicReadGraph(ReadId readCount);

    // Table that gives the vertex_descriptor corresponding to an
    // OrientedReadId. Keyed by orientedReadId.getValue().
    std::unordered_map<OrientedReadId::Int, vertex_descriptor> vertexMap;

};



// Class ReadGraph is used to store the ReadGraph in permanent
// but read-only form using MemoryMapped data structures.
class ChanZuckerberg::shasta::ReadGraph {
public:

    // An edge of the read graph.
    class Edge {
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
                CZI_ASSERT(0);
            }
        }
    };
    MemoryMapped::Vector<Edge> edges;

    // Connectivity of the read graph.
    // Stores, for each OrientedReadId, indexes into the edges vector
    // of the edges that this OrientedReadId is involved in.
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> connectivity;

    // Count the triangles that have an edge as one of the sides.
    size_t countTriangles(uint32_t edgeId) const;
};

#endif
