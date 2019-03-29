#ifndef CZI_SHASTA_READ_GRAPH_HPP
#define CZI_SHASTA_READ_GRAPH_HPP



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

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "cstdint.hpp"
#include <limits>
#include <unordered_map>

namespace ChanZuckerberg {
    namespace shasta {
        class ReadGraph;


        // Class RawReadGraph is only used inside flagCrossStrandReadGraphEdges.
        // Here, each vertex corresponds to a Read, not an OrientedRead.
        // It has one vertex per read instead of two.
        class RawReadGraph;
        class RawReadGraphEdge;
        class RawReadGraphVertex;
        using RawReadGraphBaseClass = boost::adjacency_list<
            boost::vecS,
            boost::vecS,
            boost::undirectedS,
            RawReadGraphVertex,
            RawReadGraphEdge>;
    }
}




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
    static const uint32_t infiniteDistance = std::numeric_limits<uint32_t>::max();
};



class ChanZuckerberg::shasta::RawReadGraphVertex {
public:
    RawReadGraphVertex() : strand(0) {}
    uint8_t strand;     // Or RawReadGraph::undiscovered
};



class ChanZuckerberg::shasta::RawReadGraphEdge {
public:
    bool isSameStrand;
    RawReadGraphEdge(bool isSameStrand) :
        isSameStrand(isSameStrand) {}
};



// Class RawReadGraph is only used inside flagCrossStrandReadGraphEdges.
// Here, each vertex corresponds to a Read, not an OrientedRead.
// It has one vertex per read instead of two.
class ChanZuckerberg::shasta::RawReadGraph : public RawReadGraphBaseClass {
public:
    using RawReadGraphBaseClass::RawReadGraphBaseClass;
    static const uint8_t undiscovered = 2;

    // Visitor for maximum_adjacency_search.
    class Visitor {
    public:
        Visitor(vector<edge_descriptor>& crossStrandEdges) :
            crossStrandEdges(crossStrandEdges) {}
        void initialize_vertex(vertex_descriptor, const RawReadGraph&)
        {
        }
        void start_vertex(vertex_descriptor, const RawReadGraph&)
        {
        }
        void examine_edge(edge_descriptor, const RawReadGraph&);
        void finish_vertex(vertex_descriptor, const RawReadGraph&)
        {
        }
        vector<edge_descriptor>& crossStrandEdges;
    };

    // Write in Graphviz format.
    void write(ostream&) const;
    void write(const string& fileName) const;

    // Graphviz writer.
    class Writer {
    public:
        Writer(const RawReadGraph&);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const RawReadGraph& graph;
    };
};



#endif
