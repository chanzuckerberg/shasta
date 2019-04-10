#ifndef CZI_SHASTA_ASSEMBLY_GRAPH_HPP
#define CZI_SHASTA_ASSEMBLY_GRAPH_HPP



/***************************************************************************

The assembly graph is a "compressed" version of the marker graph,
in which each linear chain of edges is replaced by a single edge.

Each assembly graph vertex  corresponds to a marker graph vertex,
but the reverse is not true, because marker graph vertices
that are internal to a linear chain of edges don't have a corresponding
vertex in the assembly graph.

***************************************************************************/

// Shasta.
#include "LongBaseSequence.hpp"
#include "MarkerGraph.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

// Standard library.
#include <limits>

namespace ChanZuckerberg {
    namespace shasta {
        class AssemblyGraph;
    }
}



class ChanZuckerberg::shasta::AssemblyGraph {
public:

    // Use the same vertex and edge ids of the marker graph.
    // We could probably get away with 32 bits.
    using VertexId = MarkerGraph::VertexId;
    using EdgeId = MarkerGraph::EdgeId;
    static const VertexId invalidVertexId = std::numeric_limits<VertexId>::max();
    static const VertexId invalidEdgeId = std::numeric_limits<EdgeId>::max();

    // The vertices of the assembly graph.
    // Each assembly graph vertex  corresponds to a marker graph vertex.
    // Indexed by vertex id in the assembly graph.
    // Contains the corresponding vertex id in the marker graph.
    MemoryMapped::Vector<VertexId> vertices;

    // The reverse complement of each vertex.
    // Indexed by VertexId.
    MemoryMapped::Vector<VertexId> reverseComplementVertex;



    // The edges of the assembly graph.
    // Indexed by edge id in the assembly graph.
    class Edge {
    public:
        VertexId source;
        VertexId target;

        // The average coverage of the marker graph edges
        // corresponding to this assembly graph edge.
        uint32_t averageCoverage;
    };
    MemoryMapped::Vector<Edge> edges;

    // The reverse complement of each edge.
    // Indexed by EdgeId.
    MemoryMapped::Vector<EdgeId> reverseComplementEdge;

    // Return true if this edge is an assembled edge.
    // To avoid assembling both strands, we only assemble
    // one edge in each reverse complemented pair.
    bool isAssembledEdge(EdgeId edgeId) const
    {
        return edgeId < reverseComplementEdge[edgeId];
    }



    // The edges that each vertex is the source of.
    // Indexed by vertex id in the assembly graph.
    // Contains edge ids that can be used as indexes into edges and edgeLists.
    MemoryMapped::VectorOfVectors<EdgeId, EdgeId> edgesBySource;

    // The edges that each vertex is the target of.
    // Indexed by vertex id in the assembly graph.
    // Contains edge ids that can be used as indexes into edges and edgeLists.
    MemoryMapped::VectorOfVectors<EdgeId, EdgeId> edgesByTarget;

    // The edge ids of global marker graph edges corresponding
    // to each edge of the assembly graph.
    // Indexed by the edge id of the assembly graph edge.
    MemoryMapped::VectorOfVectors<EdgeId, EdgeId> edgeLists;

    // A table that can be used to find the location of a marker graph
    // edge in the assembly graph, if any.
    // Indexed by the edge id in the marker graph, gives for each marker graph
    // edge a pair(EdgeId, position), where:
    // - EdgeId is the id of the assembly graph edge containing the
    //   given marker graph edge, or invalidEdgeId
    //   if the marker graph edge is not part of any assembly graph edge.
    // - Position is the index of this marker graph edge in the
    //   chain corresponding to that assembly graph edge.
    MemoryMapped::Vector< pair<EdgeId, uint32_t> > markerToAssemblyTable;

    // The assembled sequenced and repeat counts for each edge of the
    // assembly graph.
    // Indexed edge id in the assembly graph.
    LongBaseSequences sequences;
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t> repeatCounts;

    // Close and remove all open data.
    void remove();

    // Basic Graphviz output of the global assembly graph.
    void writeGraphviz(const string& fileName) const;

};



#endif
