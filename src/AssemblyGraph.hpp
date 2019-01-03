#ifndef CZI_SHASTA_ASSEMBLY_GRAPH_HPP
#define CZI_SHASTA_ASSEMBLY_GRAPH_HPP

/***************************************************************************

In the assembly graph, each vertex corresponds to a linear chain
of edges in the pruned spanning subgraph of the marker graph.
A directed vertex A->B is created if the last marker graph vertex
of the edge chain corresponding to A coincides with the
first marker graph vertex of the edge chain corresponding to B.

***************************************************************************/

// Shasta.
#include "LongBaseSequence.hpp"
#include "MarkerId.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class AssemblyGraph;
    }
}



class ChanZuckerberg::shasta::AssemblyGraph {
public:

    // Use the same vertex and edge ids of the marker graph.
    // We could probably get away with 32 bits.
    using VertexId = GlobalMarkerGraphVertexId;
    using EdgeId = GlobalMarkerGraphEdgeId;

    // The edge ids of global marker graph edges of each vertex.
    // They describe the chain (path in the marker graph)
    // corresponding to each vertex of the assembly graph.
    // Indexed by the vertex id of the assembly graph vertex.
    MemoryMapped::VectorOfVectors<EdgeId, EdgeId> vertices;

    // The edges of the assembly graph.
    // Each edge stores the vertex ids (in the assembly graph).
    class Edge {
    public:
        VertexId source;
        VertexId target;
        Edge(VertexId source, VertexId target) :
            source(source), target(target) {}
        Edge() {}
    };
    MemoryMapped::Vector<Edge> edges;

    // The edges that each vertex is the source of.
    // Contains indexes into the above edges vector.
    MemoryMapped::VectorOfVectors<VertexId, EdgeId> edgesBySource;

    // The edges that each vertex is the target of.
    // Contains indexes into the above edges vector.
    MemoryMapped::VectorOfVectors<VertexId, EdgeId> edgesByTarget;

    // A table that can be used to find the location of a marker graph
    // edge in the assembly graph, if any.
    // Indexed by the GlobalMarkerGraphEdgeId, gives for each marker graph
    // edge a pair(VertexId, position), where:
    // - VertexId is the id of the assembly graph vertex containing the
    //   given marker graph edge, or std::numeric_limits<VertexId>::max()
    //   if the marker graph edge is not part of any assembly graph vertex.
    // - Position is the index of this marker graph edge in the
    //   chain corresponding to that assembly graph vertex.
    MemoryMapped::Vector< pair<VertexId, uint32_t> > markerToAssemblyTable;

    // The assembled sequenced and repeat counts for each vertex of the
    // assembly graph.
    LongBaseSequences sequences;
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t> repeatCounts;

    // Close and remove all open data.
    void remove();

};



#endif
