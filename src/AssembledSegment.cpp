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
    edgeCoverage.clear();

    vertexSequences.clear();
    vertexRepeatCounts.clear();

    edgeSequences.clear();
    edgeRepeatCounts.clear();
    edgeOverlappingBaseCounts.clear();

    vertexOffsets.clear();

    runLengthSequence.clear();
    repeatCounts.clear();

    assembledRawSequence.clear();
    vertexRunLengthRange.clear();
    vertexRawRange.clear();
    edgeRunLengthRange.clear();
    edgeRawRange.clear();
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



void AssembledSegment::assemble()
{
    vertexRunLengthRange.resize(vertexCount);
    vertexRawRange.resize(vertexCount);
    edgeRunLengthRange.resize(edgeCount);
    edgeRawRange.resize(edgeCount);

    for(size_t i=0; ; i++) {

        // Vertex.
        vertexRunLengthRange[i].first = uint32_t(runLengthSequence.size());
        vertexRawRange[i].first = uint32_t(assembledRawSequence.size());
        for(uint32_t j=vertexAssembledPortion[i].first; j!=vertexAssembledPortion[i].second; j++) {
            const Base base = vertexSequences[i][j];
            const uint32_t repeatCount = vertexRepeatCounts[i][j];
            CZI_ASSERT(repeatCount > 0);
            runLengthSequence.push_back(base);
            repeatCounts.push_back(repeatCount);
            for(uint32_t k=0; k!=repeatCount; k++) {
                assembledRawSequence.push_back(base);
            }
        }
        vertexRunLengthRange[i].second = uint32_t(runLengthSequence.size());
        vertexRawRange[i].second = uint32_t(assembledRawSequence.size());

        // This was the last vertex.
        if(i == edgeCount) {
            break;
        }

        // Edge.
        edgeRunLengthRange[i].first = uint32_t(runLengthSequence.size());
        edgeRawRange[i].first = uint32_t(assembledRawSequence.size());
        if(edgeSequences[i].size() > 0) {
            for(uint32_t j=0; j!=uint32_t(edgeSequences[i].size()); j++) {
                const Base base = edgeSequences[i][j];
                const uint32_t repeatCount = edgeRepeatCounts[i][j];
                CZI_ASSERT(repeatCount > 0);
                runLengthSequence.push_back(base);
                repeatCounts.push_back(repeatCount);
                for(uint32_t k=0; k!=repeatCount; k++) {
                    assembledRawSequence.push_back(base);
                }
            }
        }
        edgeRunLengthRange[i].second = uint32_t(runLengthSequence.size());
        edgeRawRange[i].second = uint32_t(assembledRawSequence.size());
    }
}
