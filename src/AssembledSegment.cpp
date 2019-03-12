#include "AssembledSegment.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



void AssembledSegment::clear()
{
    assemblyGraphEdgeId = AssemblyGraph::invalidEdgeId;
    k = 0;
    vertexCount = 0;
    edgeCount = 0;

    vertexIds.clear();
    vertexCoverage.clear();

    vertexSequences.clear();
    vertexRepeatCounts.clear();

    edgeSequences.clear();
    edgeRepeatCounts.clear();
    edgeOverlappingBaseCounts.clear();

    vertexOffsets.clear();

    runLengthSequence.clear();
    repeatCounts.clear();
}



void AssembledSegment::computeVertexOffsets()
{
    vertexOffsets.resize(vertexCount);
    vertexOffsets[0] = 0;

    for(size_t i=0; i<edgeCount; i++) {
        const uint8_t overlap = edgeOverlappingBaseCounts[i];
        if(overlap > 0) {
            CZI_ASSERT(edgeSequences[i].empty());
            CZI_ASSERT(edgeRepeatCounts[i].empty());
            vertexOffsets[i+1] = uint32_t(vertexOffsets[i] + k - overlap);
        } else {
            vertexOffsets[i+1] = uint32_t(vertexOffsets[i] + k + edgeSequences[i].size());
        }
    }
}

