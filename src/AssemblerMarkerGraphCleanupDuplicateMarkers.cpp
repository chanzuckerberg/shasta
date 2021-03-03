#include "Assembler.hpp"
using namespace shasta;


// Clean up marker graph vertices that have duplicate markers
// (more than one marker on the same oriented reads).
// Such vertices are only generated when using --MarkerGraph.allowDuplicateMarkers.
void Assembler::cleanupDuplicateMarkers(
    uint64_t threadCount,
    double duplicateCoverageRatioThreshold)
{
    // Check that we have what we need.
    using CompressedVertexId = MarkerGraph::CompressedVertexId;
    MemoryMapped::Vector<CompressedVertexId>& vertexTable = markerGraph.vertexTable;
    SHASTA_ASSERT(vertexTable.isOpenWithWriteAccess);
    const MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId>& vertices = markerGraph.vertices();
    SHASTA_ASSERT(vertices.isOpen());
    const uint64_t vertexCount = vertices.size();

    cout << timestamp << "Cleaning up duplicate markers for " << vertexCount << " marker graph vertices." << endl;

    // Store information that needs to be visible to the threads.
    cleanupDuplicateMarkersData.duplicateCoverageRatioThreshold = duplicateCoverageRatioThreshold;

    // Process each vertex in multithreaded code.
    // For each vertex, we possibly change the vertexTable for the markers in that vertex (only).
    const uint64_t batchSize = 100;
    setupLoadBalancing(vertexCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction1, threadCount);

    cout << timestamp << "Cleaning up duplicate markers completed." << endl;
}



void Assembler::cleanupDuplicateMarkersThreadFunction(size_t threadId)
{
    // const double duplicateCoverageRatioThreshold = cleanupDuplicateMarkersData.duplicateCoverageRatioThreshold;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all vertices in this batch.
        for(uint64_t vertexId=begin; vertexId!=end; ++vertexId) {

        }

    }
}
