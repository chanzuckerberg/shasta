// Shasta.
#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "ConsensusCaller.hpp"
#ifdef SHASTA_HTTP_SERVER
#include "LocalMarkerGraph.hpp"
#endif
#include "timestamp.hpp"
using namespace shasta;

// Spoa.
#include "spoa/spoa.hpp"

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "chrono.hpp"
#include <map>
#include <queue>

#include "GPU.h"



// Loop over all alignments in the read graph
// to create vertices of the global marker graph.
// Throw away vertices with coverage (number of markers)
// less than minCoverage or more than maxCoverage.
// Also throw away "bad" vertices - that is, vertices
// with more than one marker on the same oriented read.
void Assembler::createMarkerGraphVerticesGpu(

    // The maximum frequency of marker k-mers to be used in
    // computing alignments.
    uint32_t maxMarkerFrequency,

    // The maximum ordinal skip to be tolerated between successive markers
    // in the alignment.
    size_t maxSkip,

    // The maximum ordinal drift to be tolerated between successive markers
    // in the alignment.
    size_t maxDrift,

    // Minimum coverage (number of markers) for a vertex
    // of the marker graph to be kept.
    size_t minCoverage,

    // Maximum coverage (number of markers) for a vertex
    // of the marker graph to be kept.
    size_t maxCoverage,

    // Number of threads. If zero, a number of threads equal to
    // the number of virtual processors is used.
    size_t threadCount
)
{

    // Flag to control debug output.
    // Only turn on for a very small test run with just a few reads.
    const bool debug = false;

    // using VertexId = MarkerGraph::VertexId;
    // using CompressedVertexId = CompressedGlobalMarkerGraphVertexId;

    const auto tBegin = steady_clock::now();
    cout << timestamp << "Begin computing marker graph vertices." << endl;

    // Check that we have what we need.
    checkReadsAreOpen();
    SHASTA_ASSERT(readFlags.isOpen);
    checkKmersAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    checkReadGraphIsOpen();

    
    //Compute unique markers
    std::unordered_map <KmerId, uint32_t> uniqueMarkersDict;

    uint32_t numUniqueMarkers = 0;
    uint32_t kmerTableSize = static_cast<uint32_t> (kmerTable.size());
    for (uint32_t j = 0; j < kmerTableSize; j++) {
        if (kmerTable[j].isMarker && kmerTable[j].isRleKmer) {
            uniqueMarkersDict[j] = numUniqueMarkers++;
        }
    }

    // Store parameters so they are accessible to the threads.
    auto& data = createMarkerGraphVerticesData;
    data.maxSkip = maxSkip;
    data.maxDrift = maxDrift;
    data.maxMarkerFrequency = maxMarkerFrequency;

    data.gpuBatchSize = shasta_getGpuBatchSize();
    data.uniqueMarkersDict = uniqueMarkersDict;
    
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Pick the batch size for computing alignments.
    size_t batchSize = alignmentCandidates.candidates.size()/threadCount;
    if (batchSize%2 > 0) {
        batchSize += 1;
    }
    if(batchSize == 0) {
        batchSize = 1000;
    }

    // Initialize computation of the global marker graph.
    data.orientedMarkerCount = markers.totalSize();
    data.disjointSetsData.createNew(
        largeDataName("tmp-DisjointSetData"),
        largeDataPageSize);
    data.disjointSetsData.reserveAndResize(data.orientedMarkerCount);
    data.disjointSetsPointer = std::make_shared<DisjointSets>(
        data.disjointSetsData.begin(),
        data.orientedMarkerCount
        );



    // Update the disjoint set data structure for each alignment
    // in the read graph.
    cout << "Begin processing " << readGraph.edges.size() << " alignments in the read graph." << endl;
    cout << timestamp << "Disjoint set computation begins." << endl;
    setupLoadBalancing(readGraph.edges.size(), batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction1Gpu, threadCount);
    
    cout << timestamp << "Shutting down processors." << endl;
    shasta_shutdownProcessors();
    
    cout << timestamp << "Disjoint set computation completed." << endl;


    // Find the disjoint set that each oriented marker was assigned to.
    cout << timestamp << "Finding the disjoint set that each oriented marker was assigned to." << endl;
    data.disjointSetTable.createNew(
        largeDataName("tmp-DisjointSetTable"),
        largeDataPageSize);
    data.disjointSetTable.reserveAndResize(data.orientedMarkerCount);
    batchSize = 1000000;
    cout << "Processing " << data.orientedMarkerCount << " oriented markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction2, threadCount);

    // Free the disjoint set data structure.
    data.disjointSetsPointer = 0;
    data.disjointSetsData.remove();


    // Debug output.
    if(debug) {
        ofstream out("DisjointSetTable-initial.csv");
        for(MarkerId markerId=0; markerId<data.orientedMarkerCount; markerId++) {
            out << markerId << "," << data.disjointSetTable[markerId] << "\n";
        }

    }



    // Count the number of markers in each disjoint set
    // and store it in data.workArea.
    // We don't want to combine this with the previous block
    // because it would significantly increase the peak memory usage.
    // This way, we allocate data.workArea only after freeing the
    // disjoint set data structure.
    cout << timestamp << "Counting the number of markers in each disjoint set." << endl;
    data.workArea.createNew(
        largeDataName("tmp-WorkArea"),
        largeDataPageSize);
    data.workArea.reserveAndResize(data.orientedMarkerCount);
    fill(data.workArea.begin(), data.workArea.end(), 0ULL);
    cout << "Processing " << data.orientedMarkerCount << " oriented markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction3, threadCount);



    // Debug output.
    if(debug) {
        ofstream out("WorkArea-initial-count.csv");
        for(MarkerId markerId=0; markerId<data.orientedMarkerCount; markerId++) {
            out << markerId << "," << data.workArea[markerId] << "\n";
        }

    }



    // At this point, data.workArea contains the number of oriented markers in
    // each disjoint set.
    // Replace it with a new numbering, counting only disjoint sets
    // with size not less than minCoverage and not greater than maxCoverage.
    // Note that this numbering is not yet the final vertex numbering,
    // as we will later remove "bad" vertices
    // (vertices with more than one marker on the same read).
    // This block is recursive and cannot be multithreaded.
    cout << timestamp << "Renumbering the disjoint sets." << endl;
    MarkerGraph::VertexId newDisjointSetId = 0ULL;
    for(MarkerGraph::VertexId oldDisjointSetId=0;
        oldDisjointSetId<data.orientedMarkerCount; ++oldDisjointSetId) {
        auto& w = data.workArea[oldDisjointSetId];
        const MarkerGraph::VertexId markerCount = w;
        if(markerCount<minCoverage || markerCount>maxCoverage) {
            w = MarkerGraph::invalidVertexId;
        } else {
            w = newDisjointSetId;
            ++newDisjointSetId;
        }
    }
    const auto disjointSetCount = newDisjointSetId;
    cout << "Kept " << disjointSetCount << " disjoint sets with coverage in the requested range." << endl;



    // Debug output.
    if(debug) {
        ofstream out("WorkArea-initial-renumbering.csv");
        for(MarkerId markerId=0; markerId<data.orientedMarkerCount; markerId++) {
            out << markerId << "," << data.workArea[markerId] << "\n";
        }
    }


    // Reassign vertices to disjoint sets using this new numbering.
    // Vertices assigned to no disjoint set will store MarkerGraph::invalidVertexId.
    // This could be multithreaded if necessary.
    cout << timestamp << "Assigning vertices to renumbered disjoint sets." << endl;
    for(MarkerGraph::VertexId markerId=0;
        markerId<data.orientedMarkerCount; ++markerId) {
        auto& d = data.disjointSetTable[markerId];
        const auto oldId = d;
        const auto newId = data.workArea[oldId];
        d = newId;
    }
    // We no longer need the workArea.
    data.workArea.remove();



    // Debug output.
    if(debug) {
        ofstream out("DisjointSetTable-renumbered.csv");
        for(MarkerId markerId=0; markerId<data.orientedMarkerCount; markerId++) {
            out << markerId << "," << data.disjointSetTable[markerId] << "\n";
        }
    }


    // At this point, data.disjointSetTable contains, for each oriented marker,
    // the disjoint set that the oriented marker is assigned to,
    // and all disjoint sets have a number of markers in the requested range.
    // Vertices not assigned to any disjoint set store a disjoint set
    // equal to MarkerGraph::invalidVertexId.



    // Gather the markers in each disjoint set.
    data.disjointSetMarkers.createNew(
        largeDataName("tmp-DisjointSetMarkers"),
        largeDataPageSize);
    cout << timestamp << "Gathering markers in disjoint sets, pass1." << endl;
    data.disjointSetMarkers.beginPass1(disjointSetCount);
    cout << timestamp << "Processing " << data.orientedMarkerCount << " oriented markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction4, threadCount);
    cout << timestamp << "Gathering markers in disjoint sets, pass2." << endl;
    data.disjointSetMarkers.beginPass2();
    cout << timestamp << "Processing " << data.orientedMarkerCount << " oriented markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction5, threadCount);
    data.disjointSetMarkers.endPass2();



    // Sort the markers in each disjoint set.
    cout << timestamp << "Sorting the markers in each disjoint set." << endl;
    setupLoadBalancing(disjointSetCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction6, threadCount);



    // Flag disjoint sets that contain more than one marker on the same read.
    data.isBadDisjointSet.createNew(
        largeDataName("tmp-IsBadDisjointSet"),
        largeDataPageSize);
    data.isBadDisjointSet.reserveAndResize(disjointSetCount);
    cout << timestamp << "Flagging bad disjoint sets." << endl;
    setupLoadBalancing(disjointSetCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction7, threadCount);
    const size_t badDisjointSetCount = std::count(
        data.isBadDisjointSet.begin(), data.isBadDisjointSet.end(), true);
    cout << "Found " << badDisjointSetCount << " bad disjoint sets "
        "with more than one marker on a single read." << endl;



    // Renumber the disjoint sets again, this time without counting the ones marked as bad.
    cout << timestamp << "Renumbering disjoint sets to remove the bad ones." << endl;
    data.workArea.createNew(
            largeDataName("tmp-WorkArea"),
            largeDataPageSize);
    data.workArea.reserveAndResize(disjointSetCount);
    newDisjointSetId = 0ULL;
    for(MarkerGraph::VertexId oldDisjointSetId=0;
            oldDisjointSetId<disjointSetCount; ++oldDisjointSetId) {
        auto& w = data.workArea[oldDisjointSetId];
        if(data.isBadDisjointSet[oldDisjointSetId]) {
            w = MarkerGraph::invalidVertexId;
        } else {
            w = newDisjointSetId;
            ++newDisjointSetId;
        }
    }
    SHASTA_ASSERT(newDisjointSetId + badDisjointSetCount == disjointSetCount);



    // Debug output.
    if(debug) {
        ofstream out("WorkArea-final-renumbering.csv");
        for(MarkerId markerId=0; markerId<disjointSetCount; markerId++) {
            out << markerId << "," << data.workArea[markerId] << "\n";
        }
    }


    // Compute the final disjoint set number for each marker.
    // That becomes the vertex id assigned to that marker.
    // This could be multithreaded.
    cout << timestamp << "Assigning vertex ids to markers." << endl;
    markerGraph.vertexTable.createNew(
            largeDataName("MarkerGraphVertexTable"),
            largeDataPageSize);
    markerGraph.vertexTable.reserveAndResize(data.orientedMarkerCount);
    for(MarkerGraph::VertexId markerId=0;
            markerId<data.orientedMarkerCount; ++markerId) {
        auto oldValue = data.disjointSetTable[markerId];
        if(oldValue == MarkerGraph::invalidVertexId) {
            markerGraph.vertexTable[markerId] = MarkerGraph::invalidCompressedVertexId;
        } else {
            markerGraph.vertexTable[markerId] = data.workArea[oldValue];
        }
    }



    // Store the disjoint sets that are not marker bad.
    // Each corresponds to a vertex of the global marker graph.
    // This could be multithreaded.
    cout << timestamp << "Gathering the markers of each vertex of the marker graph." << endl;
    markerGraph.vertices.createNew(
            largeDataName("MarkerGraphVertices"),
            largeDataPageSize);
    for(MarkerGraph::VertexId oldDisjointSetId=0;
            oldDisjointSetId<disjointSetCount; ++oldDisjointSetId) {
        if(data.isBadDisjointSet[oldDisjointSetId]) {
            continue;
        }
        markerGraph.vertices.appendVector();
        const auto markers = data.disjointSetMarkers[oldDisjointSetId];
        for(const MarkerId markerId: markers) {
            markerGraph.vertices.append(markerId);
        }
    }
    data.isBadDisjointSet.remove();
    data.workArea.remove();
    data.disjointSetMarkers.remove();
    data.disjointSetTable.remove();


    // Check that the data structures we created are consistent with each other.
    // This could be expensive. Remove when we know this code works.
    // cout << timestamp << "Checking marker graph vertices." << endl;
    // checkMarkerGraphVertices(minCoverage, maxCoverage);



    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "Computation of global marker graph vertices ";
    cout << "completed in " << tTotal << " s." << endl;

}



void Assembler::createMarkerGraphVerticesThreadFunction1Gpu(size_t threadId)
{

    array<OrientedReadId, 2> orientedReadIds;
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;

    alignment.ordinals.reserve(SHASTA_MAX_TB);

    const bool debug = false;
    auto& data = createMarkerGraphVerticesData;
    const size_t maxSkip = data.maxSkip;
    const size_t maxDrift = data.maxDrift;
    const uint32_t maxMarkerFrequency = data.maxMarkerFrequency;
    const size_t gpuBatchSize = data.gpuBatchSize;
    std::unordered_map<KmerId, uint32_t> uniqueMarkersDict = data.uniqueMarkersDict;

    const std::shared_ptr<DisjointSets> disjointSetsPointer = data.disjointSetsPointer;

    uint64_t begin, end;
    if (getNextBatch(begin, end)) {

        // We process read graph edges in pairs.
        // In each pair, the second edge is the reverse complement of the first.
        SHASTA_ASSERT((begin%2) == 0);
        SHASTA_ASSERT((end%2) == 0);
        SHASTA_ASSERT((gpuBatchSize%2) == 0);

        std::map<size_t, uint64_t> readIdLenDict;
        uint64_t numPos = 0, numReads = 0, alignCount = 0;

        // host data structures for GPU
        uint32_t* h_alignments = (uint32_t*) malloc(gpuBatchSize*SHASTA_MAX_TB*sizeof(uint32_t));
        uint32_t* h_num_traceback = (uint32_t*) malloc(gpuBatchSize*sizeof(uint32_t));

        uint64_t* batch_rid_marker_pos;
        uint64_t* batch_read_pairs;

        batch_rid_marker_pos = (uint64_t*) malloc(gpuBatchSize*SHASTA_MAX_MARKERS_PER_READ*sizeof(uint64_t));
        batch_read_pairs = (uint64_t*) malloc(2*gpuBatchSize*sizeof(uint64_t));

        // Alignment candidates that are computed on GPU
        vector<size_t> currAlignmentEdgesId;
        // Alignment candidates that fail on GPU
        vector<size_t> remainingAlignmentEdgesId;

        for (size_t first=begin; first<end; first+=gpuBatchSize) {
            size_t last = std::min(first+gpuBatchSize, end);

            if (debug)
            {
                std::lock_guard<std::mutex> lock(mutex);
                cout << "\tThreadid " << threadId << " start time: " << timestamp << endl;
            }

            // Clear vectors
            numPos = 0;
            numReads = 0;
            alignCount = 0;

            readIdLenDict.clear();


            currAlignmentEdgesId.clear();
            remainingAlignmentEdgesId.clear();

            for (size_t i=first; i<last; i+=2) {
                const ReadGraph::Edge& readGraphEdge = readGraph.edges[i];

                // Check that the next edge is the reverse complement of
                // this edge.
                {
                    const ReadGraph::Edge& readGraphNextEdge = readGraph.edges[i + 1];
                    array<OrientedReadId, 2> nextEdgeOrientedReadIds = readGraphNextEdge.orientedReadIds;
                    nextEdgeOrientedReadIds[0].flipStrand();
                    nextEdgeOrientedReadIds[1].flipStrand();
                    SHASTA_ASSERT(nextEdgeOrientedReadIds == readGraphEdge.orientedReadIds);
                }


                if(readGraphEdge.crossesStrands) {
                    continue;
                }
                orientedReadIds = readGraphEdge.orientedReadIds;
                SHASTA_ASSERT(orientedReadIds[0] < orientedReadIds[1]);

                // If either of the reads is flagged chimeric, skip it.
                if( readFlags[orientedReadIds[0].getReadId()].isChimeric ||
                        readFlags[orientedReadIds[1].getReadId()].isChimeric) {
                    continue;
                }

                ReadId rid1, rid2;
                uint64_t l1, l2;
                uint64_t batch_rid1, batch_rid2;

                rid1 = orientedReadIds[0].getValue();
                rid2 = orientedReadIds[1].getValue(); 

                l1 = getNumMarkersFromOrientedReadId(orientedReadIds[0]); 
                l2 = getNumMarkersFromOrientedReadId(orientedReadIds[1]);

                if (l1 >= SHASTA_MAX_MARKERS_PER_READ) {
                    l1 = 0;
                }
                if (l2 >= SHASTA_MAX_MARKERS_PER_READ) {
                    l2 = 0;
                }

                if ((l1 > 0) && (l2 > 0)) {
                    if (readIdLenDict.find(rid1) == readIdLenDict.end()) {
                        vector<KmerId> vec_m1 = getMarkersFromOrientedReadId(orientedReadIds[0]);
                        batch_rid1 = numReads;
                        uint64_t v, val, m;
                        v = (batch_rid1 << (32+SHASTA_LOG_MAX_MARKERS_PER_READ));
                        for (size_t l = 0; l < l1; l++) {
                            m = uniqueMarkersDict[vec_m1[l]];
                            val = v + (m << SHASTA_LOG_MAX_MARKERS_PER_READ);
                            val = val + l;
                            batch_rid_marker_pos[numPos++] = val;
                        }
                        readIdLenDict[rid1] = (batch_rid1 << 32) + l1;
                        numReads++;
                    }

                    if (readIdLenDict.find(rid2) == readIdLenDict.end()) {
                        vector<KmerId> vec_m2 = getMarkersFromOrientedReadId(orientedReadIds[1]);
                        batch_rid2 = numReads;
                        uint64_t v, val, m;
                        v = (batch_rid2 << (32+SHASTA_LOG_MAX_MARKERS_PER_READ));
                        for (size_t l = 0; l < l2; l++) {
                            m = uniqueMarkersDict[vec_m2[l]];
                            val = v + (m << SHASTA_LOG_MAX_MARKERS_PER_READ);
                            val = val + l;
                            batch_rid_marker_pos[numPos++] = val;
                        }
                        readIdLenDict[rid2] = (batch_rid2 << 32) + l2;
                        numReads++;
                    }

                    // Send read pair
                    uint64_t v1, v2;
                    v1 = readIdLenDict[rid1];
                    v2 = readIdLenDict[rid2];
                    batch_read_pairs[2*alignCount] = v1;
                    batch_read_pairs[2*alignCount+1] = v2;
                }
                else {
                    batch_read_pairs[2*alignCount] = 0;
                    batch_read_pairs[2*alignCount+1] = 0;
                    remainingAlignmentEdgesId.push_back(i);
                }
                currAlignmentEdgesId.push_back(i);
                alignCount++;
            }

            if (debug)
            {
                std::lock_guard<std::mutex> lock(mutex);
                cout << "\tThreadid " << threadId << " end time: " << timestamp << endl;
            }

            if (debug) {
                std::lock_guard<std::mutex> lock(mutex);
                fprintf(stdout, "Batchsize: %zu, Number of markers: %zu\n", (last-first), numPos); 
            }

            // find alignments on GPU
            shasta_alignBatchGPU (maxMarkerFrequency, maxSkip, maxDrift, alignCount, numPos, numReads, batch_rid_marker_pos, batch_read_pairs, h_alignments, h_num_traceback);

            for (size_t i=0; i<alignCount; i++) {
                auto edgeId = currAlignmentEdgesId[i];
                alignment.ordinals.clear();
                size_t last_addr = i*SHASTA_MAX_TB+h_num_traceback[i];
                for (size_t j=1; j<=h_num_traceback[i]; j++) {
                    uint32_t v = h_alignments[last_addr-j];
                    uint32_t l, u;
                    l = ((v << 16) >> 16);
                    u = (v >> 16);
                    alignment.ordinals.push_back(
                            array<uint32_t, 2>({(u-1), (l-1)}));
                }
                if (alignment.ordinals.size() == 0) {
                    remainingAlignmentEdgesId.push_back(edgeId);
                }
                else {
                    const ReadGraph::Edge& readGraphEdge = readGraph.edges[edgeId];
                    orientedReadIds = readGraphEdge.orientedReadIds;
                    // In the global marker graph, merge pairs
                    // of aligned markers.
                    for(const auto& p: alignment.ordinals) {
                        const uint32_t ordinal0 = p[0];
                        const uint32_t ordinal1 = p[1];
                        const MarkerId markerId0 = getMarkerId(orientedReadIds[0], ordinal0);
                        const MarkerId markerId1 = getMarkerId(orientedReadIds[1], ordinal1);
                        SHASTA_ASSERT(markers.begin()[markerId0].kmerId == markers.begin()[markerId1].kmerId);
                        disjointSetsPointer->unite(markerId0, markerId1);

                        // Also merge the reverse complemented markers.
                        // This guarantees that the marker graph remains invariant
                        // under strand swap.
                        disjointSetsPointer->unite(
                                findReverseComplement(markerId0),
                                findReverseComplement(markerId1));
                    }
                }
            }

            // Evaluate alignments that failed on GPU
            for (auto& edgeId: remainingAlignmentEdgesId) {
                const ReadGraph::Edge& readGraphEdge = readGraph.edges[edgeId];
                orientedReadIds = readGraphEdge.orientedReadIds;
                // Get the markers for the two oriented reads.
                for(size_t j=0; j<2; j++) {
                    getMarkersSortedByKmerId(orientedReadIds[j], markersSortedByKmerId[j]);
                }

                // Compute the Alignment.
                // We already know that this is a good alignment, otherwise we
                // would not have stored it.
                alignOrientedReads(
                        markersSortedByKmerId,
                        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);


                // In the global marker graph, merge pairs
                // of aligned markers.
                for(const auto& p: alignment.ordinals) {
                    const uint32_t ordinal0 = p[0];
                    const uint32_t ordinal1 = p[1];
                    const MarkerId markerId0 = getMarkerId(orientedReadIds[0], ordinal0);
                    const MarkerId markerId1 = getMarkerId(orientedReadIds[1], ordinal1);
                    SHASTA_ASSERT(markers.begin()[markerId0].kmerId == markers.begin()[markerId1].kmerId);
                    disjointSetsPointer->unite(markerId0, markerId1);

                    // Also merge the reverse complemented markers.
                    // This guarantees that the marker graph remains invariant
                    // under strand swap.
                    disjointSetsPointer->unite(
                            findReverseComplement(markerId0),
                            findReverseComplement(markerId1));
                }
            }
        }
        // free host data structures
        free(h_alignments);
        free(h_num_traceback);
        free(batch_rid_marker_pos);
        free(batch_read_pairs);
    }
}

