// For the moment this only contains a copy of
// Assembler::computeAlignments renamed
// Assembler::computeAlignmentsGpu.
// The rest is under ifdef and is a verbatim copy from src/AssemblerAlign.cpp.

// Shasta.
#include "../Assembler.hpp"
#include "../AlignmentGraph.hpp"
#include "../timestamp.hpp"
using namespace shasta;

// Standard libraries.
#include "../chrono.hpp"
#include "../iterator.hpp"
#include "../tuple.hpp"
#include <map>
#include <atomic>
//#include <cuda_runtime.h>
#include "GPU.h"

#ifdef SHASTA_BUILD_FOR_GPU

size_t numComputedAlignments=0;
size_t numGoodAlignments=0;
size_t numBadTrim=0;
size_t numFewMarkers=0;
size_t numCpuAlignments=0;

// Compute an alignment for each alignment candidate.
// Store summary information for the ones that are good enough,
// without storing details of the alignment.
void Assembler::computeAlignmentsGpu(

        // Marker frequency threshold.
        // When computing an alignment between two oriented reads,
        // marker kmers that appear more than this number of times
        // in either of the two oriented reads are discarded
        // (in both oriented reads).
        uint32_t maxMarkerFrequency,

        // The maximum ordinal skip to be tolerated between successive markers
        // in the alignment.
        size_t maxSkip,

        // The maximum relative ordinal drift to be tolerated between successive markers
        // in the alignment.
        size_t maxDrift,

        // Minimum number of alignment markers for an alignment to be used.
        size_t minAlignedMarkerCount,

        // Maximum left/right trim (in bases) for an alignment to be used.
        size_t maxTrim,

        // Number of threads. If zero, a number of threads equal to
        // the number of virtual processors is used.
        size_t threadCount
        )
{
    const auto tBegin = steady_clock::now();
    cout << timestamp << "Begin computing alignments for ";
    cout << alignmentCandidates.candidates.size() << " alignment candidates." << endl;

    // Check that we have what we need.
    checkReadsAreOpen();
    checkKmersAreOpen();
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    //Compute unique markers
    std::unordered_map <KmerId, uint32_t> uniqueMarkersDict;

    uint32_t numUniqueMarkers = 0;
    uint32_t kmerTableSize = static_cast<uint32_t> (kmerTable.size());
    for (uint32_t j = 0; j < kmerTableSize; j++) {
        if (kmerTable[j].isMarker && kmerTable[j].isRleKmer) {
            uniqueMarkersDict[j] = numUniqueMarkers++;
        }
    }

    // Initialize GPUs
    cout << timestamp << "Initializing GPU device(s)" << endl;
    int nDevices;
    size_t gpuBatchSize;
    std::tie(nDevices, gpuBatchSize) = shasta_initializeProcessors(numUniqueMarkers);
    cout << timestamp << "Initialized " << nDevices << " GPU devices" << endl;

    // Store parameters so they are accessible to the threads.
    auto& data = computeAlignmentsData;
    data.maxMarkerFrequency = maxMarkerFrequency;
    data.maxSkip = maxSkip;
    data.maxDrift = maxDrift;
    data.minAlignedMarkerCount = minAlignedMarkerCount;
    data.maxTrim = maxTrim;
    data.nDevices = nDevices;
    data.gpuBatchSize = gpuBatchSize;
    data.uniqueMarkersDict = uniqueMarkersDict;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Pick the batch size for computing alignments.
    size_t batchSize = alignmentCandidates.candidates.size()/threadCount;
    if(batchSize == 0) {
        batchSize = 1;
    }

    // Compute the alignments.
    data.threadAlignmentData.resize(threadCount);
    cout << timestamp << "Alignment computation begins." << endl;
    setupLoadBalancing(alignmentCandidates.candidates.size(), batchSize);
    runThreads(&Assembler::computeAlignmentsThreadFunctionGPU, threadCount);
    cout << timestamp << "Alignment computation completed." << endl;

    // Store alignmentInfos found by each thread in the global alignmentInfos.
    cout << timestamp << "Storing the alignment info objects." << endl;
    alignmentData.createNew(largeDataName("AlignmentData"), largeDataPageSize);
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        const vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];
        for(const AlignmentData& ad: threadAlignmentData) {
            alignmentData.push_back(ad);
        }
    }
    cout << timestamp << "Creating alignment table." << endl;
    computeAlignmentTable();

    cout << timestamp << "Offloaded " << numCpuAlignments << " alignments to CPU." << endl;
    cout << timestamp << "Computed " << numGoodAlignments << " alignments." << endl;
    cout << timestamp << "Rejected " << numBadTrim << " alignments (bad trim)." << endl;
    cout << timestamp << "Rejected " << numFewMarkers << " alignments (few markers)." << endl;

    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "Computation of alignments ";
    cout << "completed in " << tTotal << " s." << endl;
}



void Assembler::computeAlignmentsThreadFunctionGPU(size_t threadId)
{

    array<OrientedReadId, 2> orientedReadIds;
    array<OrientedReadId, 2> orientedReadIdsOppositeStrand;
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;

    alignment.ordinals.reserve(SHASTA_MAX_TB);

    const bool debug = false;
    auto& data = computeAlignmentsData;
    const uint32_t maxMarkerFrequency = data.maxMarkerFrequency;
    const size_t maxSkip = data.maxSkip;
    const size_t maxDrift = data.maxDrift;
    const size_t minAlignedMarkerCount = data.minAlignedMarkerCount;
    const size_t maxTrim = data.maxTrim;
    const size_t gpuBatchSize = data.gpuBatchSize;
    std::unordered_map<KmerId, uint32_t> uniqueMarkersDict = data.uniqueMarkersDict;

    vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];

    uint64_t begin, end;
    size_t printEvery = 100000;
    if (getNextBatch(begin, end)) {
        std::map<size_t, uint64_t> readIdLenDict;
        uint64_t numPos = 0, numReads = 0;

        // host data structures for GPU
        uint32_t* h_alignments = (uint32_t*) malloc(gpuBatchSize*SHASTA_MAX_TB*sizeof(uint32_t));
        uint32_t* h_num_traceback = (uint32_t*) malloc(gpuBatchSize*sizeof(uint32_t));

        uint64_t* batch_rid_marker_pos;
        uint64_t* batch_read_pairs;

        batch_rid_marker_pos = (uint64_t*) malloc(gpuBatchSize*SHASTA_MAX_MARKERS_PER_READ*sizeof(uint64_t));
        batch_read_pairs = (uint64_t*) malloc(2*gpuBatchSize*sizeof(uint64_t));

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

            readIdLenDict.clear();

            // Compute alignments on GPU
            for (size_t i=first; i<last; i++) {
                auto candidate = alignmentCandidates.candidates[i];
                ReadId rid1, rid2;
                uint64_t l1, l2;
                uint64_t batch_rid1, batch_rid2;
                Strand s1 = 0;
                Strand s2 = candidate.isSameStrand ? 0 : 1;
                rid1 = candidate.readIds[0];
                rid2 = candidate.readIds[1];

                l1 = getNumMarkers(rid1, s1); 
                l2 = getNumMarkers(rid2, s2);

                if (l1 >= SHASTA_MAX_MARKERS_PER_READ) {
                    l1 = 0;
                }
                if (l2 >= SHASTA_MAX_MARKERS_PER_READ) {
                    l2 = 0;
                }

                if ((l1 > 0) && (l2 > 0)) {
                    if (readIdLenDict.find(2*rid1+s1) == readIdLenDict.end()) {
                        vector<KmerId> vec_m1 = getMarkers(rid1, s1);
                        batch_rid1 = numReads;
                        uint64_t v, val, m;
                        v = (batch_rid1 << (32+SHASTA_LOG_MAX_MARKERS_PER_READ));
                        for (size_t l = 0; l < l1; l++) {
                            m = uniqueMarkersDict[vec_m1[l]];
                            val = v + (m << SHASTA_LOG_MAX_MARKERS_PER_READ);
                            val = val + l;
                            batch_rid_marker_pos[numPos++] = val;
                        }
                        readIdLenDict[2*rid1+s1] = (batch_rid1 << 32) + l1;
                        numReads++;
                    }

                    if (readIdLenDict.find(2*rid2+s2) == readIdLenDict.end()) {
                        vector<KmerId> vec_m2 = getMarkers(rid2, s2);
                        batch_rid2 = numReads;
                        uint64_t v, val, m;
                        v = (batch_rid2 << (32+SHASTA_LOG_MAX_MARKERS_PER_READ));
                        for (size_t l = 0; l < l2; l++) {
                            m = uniqueMarkersDict[vec_m2[l]];
                            val = v + (m << SHASTA_LOG_MAX_MARKERS_PER_READ);
                            val = val + l;
                            batch_rid_marker_pos[numPos++] = val;
                        }
                        readIdLenDict[2*rid2+s2] = (batch_rid2 << 32) + l2;
                        numReads++;
                    }

                    // Send read pair
                    uint64_t v1, v2;
                    v1 = readIdLenDict[2*rid1+s1];
                    v2 = readIdLenDict[2*rid2+s2];
                    batch_read_pairs[2*(i-first)] = v1;
                    batch_read_pairs[2*(i-first)+1] = v2;
                }
                else {
                    batch_read_pairs[2*(i-first)] = 0;
                    batch_read_pairs[2*(i-first)+1] = 0;
                }
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
            shasta_alignBatchGPU (maxMarkerFrequency, maxSkip, maxDrift, (last-first), numPos, numReads, batch_rid_marker_pos, batch_read_pairs, h_alignments, h_num_traceback, false);

            // Print progress
            size_t currNumComputedAlignments = __sync_fetch_and_add(&numComputedAlignments, last-first);
            if ((currNumComputedAlignments == 0) || \
                    ((currNumComputedAlignments / printEvery) != ((currNumComputedAlignments+last-first) / printEvery))) {
                std::lock_guard<std::mutex> lock(mutex);
                cout << timestamp << "Working on alignment " << printEvery*((currNumComputedAlignments+last-first)/printEvery);
                cout << " of " << alignmentCandidates.candidates.size() << endl;
            }

            // Evaluate alignments and push good alignments
            for (size_t i=first; i<last; i++) {
                auto candidate = alignmentCandidates.candidates[i];
                alignment.ordinals.clear();
                size_t last_addr = (i-first)*2+2;
                if (h_num_traceback[i-first] >= 2) {
                    {
                        uint32_t v = h_alignments[last_addr-1];
                        uint32_t l, u;
                        l = ((v << 16) >> 16);
                        u = (v >> 16);
                        alignment.ordinals.push_back(
                                array<uint32_t, 2>({(u-1), (l-1)}));
                    }
                    for (size_t j=0; j<h_num_traceback[i-first]-2; j++) {
                        alignment.ordinals.push_back(
                                array<uint32_t, 2>({0, 0}));
                    }
                    {
                        uint32_t v = h_alignments[last_addr-2];
                        uint32_t l, u;
                        l = ((v << 16) >> 16);
                        u = (v >> 16);
                        alignment.ordinals.push_back(
                                array<uint32_t, 2>({(u-1), (l-1)}));
                    }
                }
                // Re-evaluate this alignment on CPU which likely failed on GPU
                if (alignment.ordinals.size() == 0) {
                    SHASTA_ASSERT(candidate.readIds[0] < candidate.readIds[1]);

                    // Get the oriented read ids, with the first one on strand 0.
                    orientedReadIds[0] = OrientedReadId(candidate.readIds[0], 0);
                    orientedReadIds[1] = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);

                    // Get the oriented read ids for the opposite strand.
                    orientedReadIdsOppositeStrand = orientedReadIds;
                    orientedReadIdsOppositeStrand[0].flipStrand();
                    orientedReadIdsOppositeStrand[1].flipStrand();

                    // Get the markers for the two oriented reads in this candidate.
                    for(size_t j=0; j<2; j++) {
                        getMarkersSortedByKmerId(orientedReadIds[j], markersSortedByKmerId[j]);
                    }

                    // Compute the Alignment.
                    const auto t0 = std::chrono::steady_clock::now();
                    alignOrientedReads(
                            markersSortedByKmerId,
                            maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);

                    __sync_fetch_and_add(&numCpuAlignments, 1);

                    const auto t1 = std::chrono::steady_clock::now();
                    const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
                    if(t01 > 1.) {
                        std::lock_guard<std::mutex> lock(mutex);
                        cout << timestamp << "Slow alignment computation for oriented reads ";
                        cout << orientedReadIds[0] << " ";
                        cout << orientedReadIds[1] << ": ";
                        cout << t01 << " s.\n";
                    }

                    // If the alignment has too few markers skip it.
                    if(alignment.ordinals.size() < minAlignedMarkerCount) {
                        __sync_fetch_and_add(&numFewMarkers, 1);
                        continue;
                    }

                    // If the alignment has too much trim, skip it.
                    uint32_t leftTrim;
                    uint32_t rightTrim;
                    tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();
                    if(leftTrim>maxTrim || rightTrim>maxTrim) {
                        __sync_fetch_and_add(&numBadTrim, 1);
                        continue;
                    }

                    if (debug)
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        // Compute alignment info.
                        Strand s1 = 0;
                        Strand s2 = candidate.isSameStrand ? 0 : 1;
                        uint32_t rid1 = candidate.readIds[0];
                        uint32_t rid2 = candidate.readIds[1];
                        vector<KmerId> vec_m1, vec_m2;
                        vec_m1 = getMarkers(rid1, s1);
                        vec_m2 = getMarkers(rid2, s2);

                        fprintf(stdout, "[%u,%u] ", rid1, rid2); 

                        size_t numPrint = alignment.ordinals.size();
                        if (numPrint > 3) {
                            numPrint = 3;
                        }
                        for (size_t i = 0; i < numPrint; i++) {
                            fprintf(stdout, " at (%u, %u) markers (%u, %u), ", alignment.ordinals[i][0], alignment.ordinals[i][1],
                                    vec_m1[alignment.ordinals[i][0]], vec_m2[alignment.ordinals[i][1]]);
                        }
                        fprintf(stdout, "\n");
                    }

                    __sync_fetch_and_add(&numGoodAlignments, 1);
                    // If getting here, this is a good alignment.
                    threadAlignmentData.push_back(AlignmentData(candidate, alignmentInfo));
                }
                else {
                    // Compute alignment info.
                    Strand s1 = 0;
                    Strand s2 = candidate.isSameStrand ? 0 : 1;
                    uint32_t rid1 = candidate.readIds[0];
                    uint32_t rid2 = candidate.readIds[1];
                    size_t l1, l2;
                    l1 = getNumMarkers(rid1, s1);
                    l2 = getNumMarkers(rid2, s2);

                    //std::reverse(alignment.ordinals.begin(), alignment.ordinals.end());
                    alignmentInfo.create(alignment, uint32_t(l1), uint32_t(l2));

                    // If the alignment has too few markers skip it.
                    if(alignment.ordinals.size() < minAlignedMarkerCount) {
                        __sync_fetch_and_add(&numFewMarkers, 1);
                        continue;
                    }

                    // If the alignment has too much trim, skip it.
                    uint32_t leftTrim;
                    uint32_t rightTrim;
                    tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();
                    if(leftTrim>maxTrim || rightTrim>maxTrim) {
                        __sync_fetch_and_add(&numBadTrim, 1);
                        continue;
                    }

                    if (debug)
                    {
                        std::lock_guard<std::mutex> lock(mutex);

                        vector<KmerId> vec_m1 = getMarkers(rid1, s1);
                        vector<KmerId> vec_m2 = getMarkers(rid2, s2);

                        fprintf(stdout, "[%u,%u] ", rid1, rid2); 

                        size_t numPrint = alignment.ordinals.size();
                        if (numPrint > 3) {
                            numPrint = 3;
                        }
                        for (size_t i = 0; i < numPrint; i++) {
                            fprintf(stdout, " at (%u, %u) markers (%u, %u), ", alignment.ordinals[i][0], alignment.ordinals[i][1],
                                    vec_m1[alignment.ordinals[i][0]], vec_m2[alignment.ordinals[i][1]]);
                        }
                        fprintf(stdout, "\n");
                    }

                    __sync_fetch_and_add(&numGoodAlignments, 1);
                    // If getting here, this is a good alignment.
                    threadAlignmentData.push_back(AlignmentData(candidate, alignmentInfo));
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

#endif
