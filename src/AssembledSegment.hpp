#ifndef CZI_SHASTA_ASSEMBLED_SEGMENT_HPP
#define CZI_SHASTA_ASSEMBLED_SEGMENT_HPP

// Shasta.
#include "AssemblyGraph.hpp"
#include "Base.hpp"
#include "Coverage.hpp"
#include "MarkerGraph.hpp"

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

    // The marker graph vertices and edges of the chain corresponding to this segment.
    vector<MarkerGraph::VertexId> vertexIds;
    vector<MarkerGraph::EdgeId> edgeIds;

    // Vertex and edge coverage.
    vector<uint32_t> vertexCoverage;
    vector<uint32_t> edgeCoverage;

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

    // Compute, for each vertex, the portion of vertex sequence that contributes
    // to the assembly. This is the portion that does not overlap a vertex with greater coverage.
    // (Break ties using vertex ids).
    // An edge with overlapping markers does not contribute to the assembly.
    // An edge with at least one intervening base contributes all of its bases
    // to the assembly.
    vector< pair<uint32_t, uint32_t> > vertexAssembledPortion;
    void computeVertexAssembledPortion();

    // The assembled run-length sequence  and repeat counts.
    vector<Base> runLengthSequence;
    vector<uint32_t> repeatCounts;
    vector<Base> assembledRawSequence;

    // Keep track of the range each vertex and edge contributes.
    vector< pair<uint32_t, uint32_t> > vertexRunLengthRange;
    vector< pair<uint32_t, uint32_t> > vertexRawRange;
    vector< pair<uint32_t, uint32_t> > edgeRunLengthRange;
    vector< pair<uint32_t, uint32_t> > edgeRawRange;
    void assemble();

    // Put back into default-constructed state
    // (except for vector capacities).
    void clear();

    // Write out details in html.
    void writeHtml(ostream&) const;


    // Coverage data is computed optionally under control of the storeCoverageData
    // argument to Assembler::assembleAssemblyGraphEdge.
    // Indexed by [i][position] where:
    // - i is the index for the vertex or edge in the chain.
    // - position is in run-length coordinates.
    vector< vector < vector<CompressedCoverageData> > > vertexCoverageData;
    vector< vector < vector<CompressedCoverageData> > > edgeCoverageData;

    // Coverage data for assembled sequence.
    // Indexed by the assembled position in run-length coordinates.
    vector< vector<CompressedCoverageData> > assembledCoverageData;



    // Python-callable accessors.
    // Unless otherwise specified, positions are in run-length coordinates.
    size_t size() const
    {
        return runLengthSequence.size();
    }
    char getBase(uint32_t position) const
    {
        SHASTA_ASSERT(position < runLengthSequence.size());
        return runLengthSequence[position].character();
    }
    uint32_t getRepeatCount(uint32_t position) const
    {
        SHASTA_ASSERT(position < repeatCounts.size());
        return repeatCounts[position];
    }
    const vector<CompressedCoverageData> getCoverageData(uint32_t position) const
    {
        SHASTA_ASSERT(position < assembledCoverageData.size());
        return assembledCoverageData[position];
    }

};




#endif

