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

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store information that needs to be visible to the threads.
    cleanupDuplicateMarkersData.duplicateCoverageRatioThreshold = duplicateCoverageRatioThreshold;
    cleanupDuplicateMarkersData.badVertexCount = 0;

    // Process each vertex in multithreaded code.
    // For each vertex, we possibly change the vertexTable for the markers in that vertex (only).
    const uint64_t batchSize = 100;
    setupLoadBalancing(vertexCount, batchSize);
    runThreads(&Assembler::cleanupDuplicateMarkersThreadFunction, threadCount);

    cout << "Found " << cleanupDuplicateMarkersData.badVertexCount <<
        " vertices with duplicate markers." << endl;
    cout << timestamp << "Cleaning up duplicate markers completed." << endl;
}



void Assembler::cleanupDuplicateMarkersThreadFunction(size_t threadId)
{
    // const double duplicateCoverageRatioThreshold = cleanupDuplicateMarkersData.duplicateCoverageRatioThreshold;
    uint64_t badVertexCount = 0;

    // The pairs (orientedReadId, marker ordinal) for the current vertex.
    using MarkerPair = pair<OrientedReadId, uint32_t>;
    vector<MarkerPair> markerPairs;

    // Vector of flags that tells us which MarkerPair's are duplicate (duplicate in OrientedReadId only).
    vector<bool> isDuplicateOrientedReadId;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all vertices in this batch.
        for(uint64_t vertexId=begin; vertexId!=end; ++vertexId) {

            // If this vertex does not have duplicate markers, skip it.
            if(not isBadMarkerGraphVertex(vertexId)) {
                continue;
            }

            // This vertex has duplicate markers (more than one marker on the
            // same oriented read).
            ++badVertexCount;

            // Get the pairs (orientedReadId, marker ordinal) for this vertex.
            const span<MarkerId> markerIds = markerGraph.getVertexMarkerIds(vertexId);
            SHASTA_ASSERT(markerIds.size() > 1);
            markerPairs.clear();
            for(const MarkerId markerId: markerIds) {
                markerPairs.push_back(findMarkerId(markerId));
            }

            // Find the ones that are duplicate.
            // We take advantage of the fact that the pairs are sorted by OrientedReadId.
            isDuplicateOrientedReadId.resize(markerPairs.size());
            fill(isDuplicateOrientedReadId.begin(), isDuplicateOrientedReadId.end(), false);
            for(uint64_t i=1; i<markerPairs.size(); i++) {
                if(markerPairs[i-1].first == markerPairs[i].first) {
                    isDuplicateOrientedReadId[i-1] = true;
                    isDuplicateOrientedReadId[i] = true;
                }
            }

        }
    }

    // Increment the total number of vertices with duplicate markers.
    __sync_fetch_and_add(&cleanupDuplicateMarkersData.badVertexCount, badVertexCount);
}
