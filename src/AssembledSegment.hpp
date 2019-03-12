#ifndef CZI_SHASTA_ASSEMBLED_SEGMENT_HPP
#define CZI_SHASTA_ASSEMBLED_SEGMENT_HPP

// Shasta.
#include "AssemblyGraph.hpp"
#include "Base.hpp"
#include "MarkerId.hpp"

// Standard library.
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class AssembledSegment;
    }
}



// Class to describe a sequence segment assembled
// from an edge of the assemblygraph.
class ChanZuckerberg::shasta::AssembledSegment {
public:

    // The edge id of the assembly graph edge corresponding to this segment.
    AssemblyGraph::EdgeId assemblyGraphEdgeId;

    // The number of marker graph vertices and edges corresponding to this segment.
    // Since this is a linear chian, the number of vertices equals the number of edges
    // plus one.
    size_t edgeCount;
    size_t vertexCount;


    // The assembled run-length sequence  and repeat counts.
    vector<Base> runLengthSequence;
    vector<uint32_t> repeatCounts;

    // The marker graph vertices of the chain corresponding to this segment.
    vector<GlobalMarkerGraphVertexId> vertexIds;

    // Vertex coverage.
    vector<uint32_t> vertexCoverage;

    // Put back into default-constructed state
    // (except for vector capacities).
    void clear();
};




#endif

