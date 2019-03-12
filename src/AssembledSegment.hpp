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

    // The length of a marker.
    size_t k;

    // The number of marker graph vertices and edges corresponding to this segment.
    // Since this is a linear chian, the number of vertices equals the number of edges
    // plus one.
    size_t edgeCount;
    size_t vertexCount;

    // The marker graph vertices of the chain corresponding to this segment.
    vector<GlobalMarkerGraphVertexId> vertexIds;

    // Vertex coverage.
    vector<uint32_t> vertexCoverage;

    // The consensus sequences and repeat counts for the vertices in the chain.
    vector< vector<Base> > vertexSequences;
    vector< vector<uint32_t> > vertexRepeatCounts;

    // The consensus sequences and repeat counts for the edges in the chain.
    vector< vector<Base> > edgeSequences;
    vector< vector<uint32_t> > edgeRepeatCounts;
    vector<uint8_t> edgeOverlappingBaseCounts;

    // Vertex offsets.
    // A vertex offset is the position of the first base
    // of the vertex consensus sequence (run-length)
    // relative to the first base of assembled sequence (run-length).
    vector<uint32_t> vertexOffsets;
    void computeVertexOffsets();


    // The assembled run-length sequence  and repeat counts.
    vector<Base> runLengthSequence;
    vector<uint32_t> repeatCounts;

    // Put back into default-constructed state
    // (except for vector capacities).
    void clear();
};




#endif

