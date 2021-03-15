#include "Assembler.hpp"
using namespace shasta;


// Clean up marker graph vertices that have duplicate markers
// (more than one marker on the same oriented reads).
// Such vertices are only generated when using --MarkerGraph.allowDuplicateMarkers.
void Assembler::cleanupDuplicateMarkers(
    uint64_t threadCount,
    double pattern1Threshold,
    bool pattern1CreateNewVertices)
{
    const bool debug = false;

    // Check that we have what we need.
    SHASTA_ASSERT(markers.isOpen());
    using CompressedVertexId = MarkerGraph::CompressedVertexId;
    MemoryMapped::Vector<CompressedVertexId>& vertexTable = markerGraph.vertexTable;
    SHASTA_ASSERT(vertexTable.isOpenWithWriteAccess);
    const MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId>& vertices = markerGraph.vertices();
    SHASTA_ASSERT(vertices.isOpen());
    const uint64_t vertexCount = vertices.size();
    SHASTA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    SHASTA_ASSERT(markerGraph.reverseComplementVertex.size() == vertexCount);

    cout << timestamp << "Cleaning up duplicate markers for " << vertexCount << " marker graph vertices." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store information that needs to be visible to the threads.
    cleanupDuplicateMarkersData.pattern1Threshold = pattern1Threshold;
    cleanupDuplicateMarkersData.pattern1CreateNewVertices = pattern1CreateNewVertices;
    cleanupDuplicateMarkersData.badVertexCount = 0;
    cleanupDuplicateMarkersData.pattern1Count = 0;
    cleanupDuplicateMarkersData.removedCount = 0;
    cleanupDuplicateMarkersData.nextVertexId = vertexCount;

    // Process each vertex in multithreaded code.
    // For each vertex, we possibly change the vertexTable for the markers in that vertex (only).
    // Some vertices can disappear completely.
    const uint64_t batchSize = 100;
    setupLoadBalancing(vertexCount, batchSize);
    runThreads(&Assembler::cleanupDuplicateMarkersThreadFunction, threadCount);

    cout << "Found " << cleanupDuplicateMarkersData.badVertexCount <<
        " vertices with duplicate markers." << endl;
    cout << "Pattern 1 vertex count: " << cleanupDuplicateMarkersData.pattern1Count << endl;
    cout << "Unprocessed (removed) vertex count: " << cleanupDuplicateMarkersData.removedCount << endl;

    // Renumber the vertex table to make sure vertices are numbered contiguously starting at 0.
    if(debug) {
        cout << "Maximum vertex id before renumbering of the vertex table " << cleanupDuplicateMarkersData.nextVertexId - 1 << endl;
    }
    const MarkerGraph::VertexId maxVertexId =
        markerGraph.renumberVertexTable(threadCount, cleanupDuplicateMarkersData.nextVertexId - 1);
    if(debug) {
        cout << "Maximum vertex id after renumbering of the vertex table " << maxVertexId << endl;
    }

    // Now we can recreate the vertices in the marker graph.
    markerGraph.createVerticesFromVertexTable(
        threadCount, maxVertexId);
    if(debug) {
        cout << "New number of vertices is " << markerGraph.vertices().size() << endl;
    }



    // Sanity check.
    if(debug) {
        for(MarkerGraph::VertexId vertexId=0; vertexId<markerGraph.vertices().size(); vertexId++) {
            if(markerGraph.vertices().size(vertexId) == 0) {
                cout << "Failing vertex id " << vertexId << endl;
            }
            SHASTA_ASSERT(markerGraph.vertices().size(vertexId) > 0);
        }
    }



    // Finally, recreate the reverse complement vertices.
    findMarkerGraphReverseComplementVertices(threadCount);


    cout << timestamp << "Cleaning up duplicate markers completed." << endl;
    cout << "Number of marker graph vertices is now " << markerGraph.vertices().size() << endl;
}



void Assembler::cleanupDuplicateMarkersThreadFunction(size_t threadId)
{
    const bool debug = true;
    ofstream out;
    if(debug) {
        out.open("cleanupDuplicateMarkers-" + to_string(threadId) + ".threadLog");
    }

    const double pattern1Threshold = cleanupDuplicateMarkersData.pattern1Threshold;
    const bool pattern1CreateNewVertices = cleanupDuplicateMarkersData.pattern1CreateNewVertices;

    uint64_t badVertexCount = 0;
    uint64_t pattern1Count = 0;
    uint64_t removedCount = 0;

    // The pairs (orientedReadId, marker ordinal) for the current vertex.
    using MarkerPair = pair<OrientedReadId, uint32_t>;
    vector<MarkerPair> markerPairs;

    // Vector of flags that tells us which MarkerPair's are duplicate (duplicate in OrientedReadId only).
    vector<bool> isDuplicateOrientedReadId;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all vertices in this batch.
        for(MarkerGraph::VertexId vertexId=begin; vertexId!=end; ++vertexId) {

            // Process one pair of reverse complemented vertices at a time.
            const MarkerGraph::VertexId vertexIdRc = markerGraph.reverseComplementVertex[vertexId];
            if(vertexIdRc < vertexId) {
                continue;
            }

            // If this vertex does not have duplicate markers, skip it.
            if(not isBadMarkerGraphVertex(vertexId)) {
                continue;
            }

            // This vertex has duplicate markers (more than one marker on the
            // same oriented read).
            if(vertexId == vertexIdRc) {
                ++badVertexCount;   // Unusual/exceptional case.
            } else {
                badVertexCount += 2;
            }
            if(debug) {
                if(vertexId == vertexIdRc) {
                    out << "Working on self-complementary vertex " <<
                        vertexId << endl;
                } else {
                    out << "Working on vertex " <<
                        vertexId << " and its reverse complement " << vertexIdRc << endl;
                }
            }

            // Get the pairs (orientedReadId, marker ordinal) for this vertex.
            const span<MarkerId> markerIds = markerGraph.getVertexMarkerIds(vertexId);
            const uint64_t markerCount = markerIds.size();
            SHASTA_ASSERT(markerCount > 1);
            markerPairs.clear();
            for(const MarkerId markerId: markerIds) {
                markerPairs.push_back(findMarkerId(markerId));
            }

            // Find the ones that are duplicate.
            // We take advantage of the fact that the pairs are sorted by OrientedReadId.
            isDuplicateOrientedReadId.resize(markerPairs.size());
            fill(isDuplicateOrientedReadId.begin(), isDuplicateOrientedReadId.end(), false);
            for(uint64_t i=1; i<markerCount; i++) {
                if(markerPairs[i-1].first == markerPairs[i].first) {
                    isDuplicateOrientedReadId[i-1] = true;
                    isDuplicateOrientedReadId[i] = true;
                }
            }
            const uint64_t duplicateCount =
                std::count(isDuplicateOrientedReadId.begin(), isDuplicateOrientedReadId.end(), true);
            if(debug) {
                out << duplicateCount << " duplicate markers out of " << markerCount << endl;
                for(uint64_t i=0; i<markerCount; i++) {
                    const auto& p = markerPairs[i];
                    out << p.first << " " << p.second;
                    if(isDuplicateOrientedReadId[i]) {
                        out << " duplicate";
                    }
                    out << endl;
                }
            }
            SHASTA_ASSERT(duplicateCount > 0);



            // Pattern 1: the number of duplicate markers is small.
            const double duplicateRatio = double(duplicateCount) / double(markerCount);
            if(duplicateRatio < pattern1Threshold) {
                if(debug) {
                    out << "Vertex " << vertexId << " processed as pattern 1 vertex." << endl;
                }
                SHASTA_ASSERT(duplicateCount < markerCount);
                if(vertexId == vertexIdRc) {
                    ++pattern1Count;   // Unusual/exceptional case.
                } else {
                    pattern1Count += 2;
                }
                cleanupDuplicateMarkersPattern1(vertexId,
                    pattern1CreateNewVertices, markerPairs, isDuplicateOrientedReadId, debug, out);
                continue;
            }


            // If we get here, for lack of a better solution, remove the vertex.
            if(vertexId == vertexIdRc) {
                ++removedCount;   // Unusual/exceptional case.
            } else {
                removedCount += 2;
            }
            if(debug) {
                out << "Vertex " << vertexId << " not processed, removed instead." << endl;
            }
            for(const auto& p: markerPairs) {
                const MarkerId markerId = getMarkerId(p.first, p.second);
                const MarkerId markerIdRc = getReverseComplementMarkerId(p.first, p.second);
                markerGraph.vertexTable[markerId] = MarkerGraph::invalidCompressedVertexId;
                markerGraph.vertexTable[markerIdRc] = MarkerGraph::invalidCompressedVertexId;
            }
        }
    }

    // Increment global counts.
    __sync_fetch_and_add(&cleanupDuplicateMarkersData.badVertexCount, badVertexCount);
    __sync_fetch_and_add(&cleanupDuplicateMarkersData.pattern1Count, pattern1Count);
    __sync_fetch_and_add(&cleanupDuplicateMarkersData.removedCount, removedCount);
}



void Assembler::cleanupDuplicateMarkersPattern1(
    MarkerGraph::VertexId vertexId,
    bool createNewVertices,
    vector< pair<OrientedReadId, uint32_t> > &markerPairs,
    vector<bool>& isDuplicateOrientedReadId,
    bool debug,
    ostream& out)
{
    if(debug) {
        out << "Processing pattern 1 vertex " << vertexId << endl;
    }
    const uint64_t markerCount = markerPairs.size();
    SHASTA_ASSERT(isDuplicateOrientedReadId.size() == markerCount);

    // Loop over markers on this vertex.
    for(uint64_t i=0; i<markerCount; i++) {

        // If not duplicate, don't do anything.
        if(not isDuplicateOrientedReadId[i]) {
            continue;
        }

        const pair<OrientedReadId, uint32_t>& p = markerPairs[i];
        const MarkerId markerId = getMarkerId(p.first, p.second);
        const MarkerId markerIdRc = getReverseComplementMarkerId(p.first, p.second);

        if(createNewVertices) {

            // Assign it to a new vertex.
            markerGraph.vertexTable[markerId] =
                __sync_fetch_and_add(&cleanupDuplicateMarkersData.nextVertexId, 1);
            markerGraph.vertexTable[markerIdRc] =
                __sync_fetch_and_add(&cleanupDuplicateMarkersData.nextVertexId, 1);
        } else {

            // Take this marker out of the current vertex, without
            // assigning it to a new vertex.
            markerGraph.vertexTable[markerId] = MarkerGraph::invalidCompressedVertexId;
            markerGraph.vertexTable[markerIdRc] = MarkerGraph::invalidCompressedVertexId;
        }
    }
}

