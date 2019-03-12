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



void AssembledSegment::computeVertexAssembledPortion()
{
    // Compute, for each vertex, the portion of vertex sequence that contributes
    // to the assembly. This is the portion that does not overlap a vertex with greater coverage.
    // (Break ties using vertex ids).
    // An edge with overlapping markers does not contribute to the assembly.
    // An edge with at least one intervening base contributes all of its bases
    // to the assembly.

    vertexAssembledPortion.resize(vertexCount);

    for(int i=0; i<int(vertexCount); i++) {

        // Check previous vertices.
        vertexAssembledPortion[i].first = 0;
        for(int j=i-1; j>=0; j--) {
            if(vertexOffsets[j]+k < vertexOffsets[i]) {
                break;
            }
            if(vertexCoverage[j]>vertexCoverage[i] ||
                (vertexCoverage[j]==vertexCoverage[i] && vertexIds[j]<vertexIds[i])) {
                vertexAssembledPortion[i].first =
                    vertexOffsets[j] + uint32_t(k) - vertexOffsets[i];
                break;
            }
        }

        // Check following vertices.
        vertexAssembledPortion[i].second = uint32_t(k);
        for(int j=i+1; j<int(vertexCount); j++) {
            if(vertexOffsets[i]+k < vertexOffsets[j]) {
                break;
            }
            if(vertexCoverage[j]>vertexCoverage[i] ||
                (vertexCoverage[j]==vertexCoverage[i] && vertexIds[j]<vertexIds[i])) {
                vertexAssembledPortion[i].second = vertexOffsets[j] - vertexOffsets[i];
                break;
            }
        }

        // Handle the case of a vertex that contributes nothing.
        if(vertexAssembledPortion[i].second <= vertexAssembledPortion[i].first) {
            vertexAssembledPortion[i].first = 0;
            vertexAssembledPortion[i].second = 0;
        }
        CZI_ASSERT(vertexAssembledPortion[i].second <= k);
    }
}

