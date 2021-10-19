#include "assembleMarkerGraphPath.hpp"
#include "AssembledSegment.hpp"
using namespace shasta;



void shasta::assembleMarkerGraphPath(
    uint64_t readRepresentation,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    const span<const MarkerGraph::EdgeId>& markerGraphPath,
    bool storeCoverageData,
    AssembledSegment& assembledSegment)
{
    assembledSegment.clear();
    assembledSegment.k = k;

    assembledSegment.edgeCount = markerGraphPath.size();
    assembledSegment.vertexCount = assembledSegment.edgeCount + 1;
    assembledSegment.edgeIds.resize(assembledSegment.edgeCount);
    copy(markerGraphPath.begin(), markerGraphPath.end(), assembledSegment.edgeIds.begin());

    // Gather the vertices of this path in the marker graph.
    assembledSegment.vertexIds.reserve(assembledSegment.vertexCount);
    for(const MarkerGraph::EdgeId edgeId: assembledSegment.edgeIds) {
        const MarkerGraph::Edge& edge =
            markerGraph.edges[edgeId];
        assembledSegment.vertexIds.push_back(edge.source);
    }
    const MarkerGraph::Edge& lastEdge =
        markerGraph.edges[assembledSegment.edgeIds[assembledSegment.edgeIds.size()-1]];
    assembledSegment.vertexIds.push_back(lastEdge.target);

    // Get vertex coverage.
    assembledSegment.vertexCoverage.resize(assembledSegment.vertexCount);
    for(size_t i=0; i<assembledSegment.vertexCount; i++) {
        assembledSegment.vertexCoverage[i] = uint32_t(markerGraph.vertexCoverage(assembledSegment.vertexIds[i]));
    }

    // Edge coverage.
    assembledSegment.edgeCoverage.resize(assembledSegment.edgeCount);
    for(size_t i=0; i<assembledSegment.edgeCount; i++) {
        assembledSegment.edgeCoverage[i] =
            uint32_t(markerGraph.edgeMarkerIntervals.size(assembledSegment.edgeIds[i]));
    }



    // Extract consensus sequence for the vertices of the chain.
    assembledSegment.vertexSequences.resize(assembledSegment.vertexCount);
    assembledSegment.vertexRepeatCounts.resize(assembledSegment.vertexCount);
    for(size_t i=0; i<assembledSegment.vertexCount; i++) {

        // Get the sequence.
        const MarkerId firstMarkerId = markerGraph.getVertexMarkerIds(assembledSegment.vertexIds[i])[0];
        const CompressedMarker& firstMarker = markers.begin()[firstMarkerId];
        const KmerId kmerId = firstMarker.kmerId;
        const Kmer kmer(kmerId, k);

        // Get the repeat counts.
        const auto& storedConsensus = markerGraph.vertexRepeatCounts.begin() + k * assembledSegment.vertexIds[i];

        // Store in the AssembledSegment.
        assembledSegment.vertexSequences[i].resize(k);
        assembledSegment.vertexRepeatCounts[i].resize(k);
        for(size_t j=0; j<k; j++) {
            assembledSegment.vertexSequences[i][j] = kmer[j];
            assembledSegment.vertexRepeatCounts[i][j] = storedConsensus[j];
        }
    }



    // Extract consensus sequence for the edges of the chain.
    assembledSegment.edgeSequences.resize(assembledSegment.edgeCount);
    assembledSegment.edgeRepeatCounts.resize(assembledSegment.edgeCount);
    assembledSegment.edgeOverlappingBaseCounts.resize(assembledSegment.edgeCount);
    for(size_t i=0; i<assembledSegment.edgeCount; i++) {

        const auto& storedConsensus = markerGraph.edgeConsensus[assembledSegment.edgeIds[i]];
        assembledSegment.edgeSequences[i].resize(storedConsensus.size());
        assembledSegment.edgeRepeatCounts[i].resize(storedConsensus.size());
        for(size_t j=0; j<storedConsensus.size(); j++) {
            assembledSegment.edgeSequences[i][j] = storedConsensus[j].first;
            assembledSegment.edgeRepeatCounts[i][j] = storedConsensus[j].second;
        }
        assembledSegment.edgeOverlappingBaseCounts[i] =
            markerGraph.edgeConsensusOverlappingBaseCount[assembledSegment.edgeIds[i]];
    }



    // Extract coverage data for vertices and edges.
    if(storeCoverageData) {

        // Check that coverage data is available.
        if( !markerGraph.vertexCoverageData.isOpen() ||
            !markerGraph.edgeCoverageData.isOpen()) {
            throw runtime_error("Coverage data is not accessible.");
        }

        // Vertices.
        assembledSegment.vertexCoverageData.resize(assembledSegment.vertexCount);
        for(size_t i=0; i<assembledSegment.vertexCount; i++) {
            const auto& input = markerGraph.vertexCoverageData[assembledSegment.vertexIds[i]];
            auto& output = assembledSegment.vertexCoverageData[i];
            output.resize(k);
            for(const pair<uint32_t, CompressedCoverageData>& p: input) {
                const uint32_t position = p.first;
                SHASTA_ASSERT(position < k);
                const CompressedCoverageData& cd = p.second;
                output[position].push_back(cd);
            }
        }

        // Edges.
        assembledSegment.edgeCoverageData.resize(assembledSegment.edgeCount);
        for(size_t i=0; i<assembledSegment.edgeCount; i++) {
            const auto& input = markerGraph.edgeCoverageData[assembledSegment.edgeIds[i]];
            auto& output = assembledSegment.edgeCoverageData[i];
            for(const pair<uint32_t, CompressedCoverageData>& p: input) {
                const uint32_t position = p.first;
                if(position >= output.size()) {
                    output.resize(position+1);
                }
                const CompressedCoverageData& cd = p.second;
                output[position].push_back(cd);
            }
        }
    }



    // Compute vertex offsets.
    // A vertex offset is the position of the first base
    // of the vertex consensus sequence (run-length)
    // relative to the first base of assembled sequence (run-length).
    assembledSegment.computeVertexOffsets();

    // Compute, for each vertex, the portion of vertex sequence that contributes
    // to the assembly. This is the portion that does not overlap a vertex with greater coverage.
    // (Break ties using vertex ids).
    // An edge with overlapping markers does not contribute to the assembly.
    // An edge with at least one intervening base contributes all of its bases
    // to the assembly.
    assembledSegment.computeVertexAssembledPortion();

    // Assemble run-length sequence and raw sequence.
    // Keep track of the range each vertex and edge contributes.
    assembledSegment.assemble();

}
