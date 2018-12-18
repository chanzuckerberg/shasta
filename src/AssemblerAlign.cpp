// shasta.
#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Standard libraries.
#include "chrono.hpp"
#include "iterator.hpp"
#include "tuple.hpp"



// Compute a marker alignment of two oriented reads.
void Assembler::alignOrientedReads(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    size_t maxSkip, // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer
)
{
    alignOrientedReads(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1),
        maxSkip, maxVertexCountPerKmer
        );
}
void Assembler::alignOrientedReads(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    size_t maxSkip, // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer
)
{
    checkReadsAreOpen();
    checkReadNamesAreOpen();
    checkMarkersAreOpen();


    // Get the markers sorted by kmerId.
    vector<MarkerWithOrdinal> markers0SortedByKmerId;
    vector<MarkerWithOrdinal> markers1SortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markers0SortedByKmerId);
    getMarkersSortedByKmerId(orientedReadId1, markers1SortedByKmerId);

    // Call the lower level function.
    AlignmentGraph graph;
    Alignment alignment;
    const bool debug = true;
    alignOrientedReads(
        markers0SortedByKmerId,
        markers1SortedByKmerId,
        maxSkip, maxVertexCountPerKmer, debug, graph, alignment);

    // Compute the AlignmentInfo.
    const AlignmentInfo alignmentInfo(alignment);
    uint32_t leftTrim;
    uint32_t rightTrim;
    tie(leftTrim, rightTrim) = computeTrim(
        orientedReadId0,
        orientedReadId1,
        alignmentInfo);
    cout << orientedReadId0 << " has " << reads[orientedReadId0.getReadId()].baseCount;
    cout << " bases and " << markers0SortedByKmerId.size() << " markers." << endl;
    cout << orientedReadId1 << " has " << reads[orientedReadId1.getReadId()].baseCount;
    cout << " bases and " << markers1SortedByKmerId.size() << " markers." << endl;
    cout << "The alignment has " << alignmentInfo.markerCount;
    cout << " markers. Left trim " << leftTrim;
    cout << " markers, right trim " << rightTrim << " markers." << endl;

    // For convenience, also write the two oriented reads.
    ofstream fasta("AlignedOrientedReads.fasta");
    writeOrientedRead(orientedReadId0, fasta);
    writeOrientedRead(orientedReadId1, fasta);
}



// This lower level version takes as input vectors of
// markers already sorted by kmerId.
void Assembler::alignOrientedReads(
    const vector<MarkerWithOrdinal>& markers0SortedByKmerId,
    const vector<MarkerWithOrdinal>& markers1SortedByKmerId,
    size_t maxSkip,             // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer
)
{
    // Compute the alignment.
    AlignmentGraph graph;
    Alignment alignment;
    const bool debug = true;
    alignOrientedReads(
        markers0SortedByKmerId,
        markers1SortedByKmerId,
        maxSkip, maxVertexCountPerKmer, debug, graph, alignment);
}



void Assembler::alignOrientedReads(
    const vector<MarkerWithOrdinal>& markers0SortedByKmerId,
    const vector<MarkerWithOrdinal>& markers1SortedByKmerId,
    size_t maxSkip,             // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer,
    bool debug,
    AlignmentGraph& graph,
    Alignment& alignment
)
{
    align(markers0SortedByKmerId, markers1SortedByKmerId,
        maxSkip, maxVertexCountPerKmer, graph, debug, alignment);
}



// Compute marker alignments of an oriented read with all reads
// for which we have an alignment.
void Assembler::alignOverlappingOrientedReads(
    ReadId readId, Strand strand,
    size_t maxSkip,                 // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer,
    size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
    size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
)
{
    alignOverlappingOrientedReads(
        OrientedReadId(readId, strand),
        maxSkip, maxVertexCountPerKmer, minAlignedMarkerCount, maxTrim);
}



void Assembler::alignOverlappingOrientedReads(
    OrientedReadId orientedReadId0,
    size_t maxSkip,                 // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer,
    size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
    size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
)
{
    // Check that we have what we need.
    checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    // Get the markers for orientedReadId0.
    vector<MarkerWithOrdinal> markers0SortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markers0SortedByKmerId);

    // Loop over all alignments involving this oriented read.
    vector<MarkerWithOrdinal> markers1SortedByKmerId;
    size_t goodAlignmentCount = 0;
    for(const uint64_t i: alignmentTable[orientedReadId0.getValue()]) {
        const AlignmentData& ad = alignmentData[i];

        // Get the other oriented read involved in this alignment.
        const OrientedReadId orientedReadId1 = ad.getOther(orientedReadId0);

        // Get the markers for orientedReadId1.
        getMarkersSortedByKmerId(orientedReadId1, markers1SortedByKmerId);

        // Compute the alignment.
        AlignmentGraph graph;
        Alignment alignment;
        const bool debug = false;
        alignOrientedReads(
            markers0SortedByKmerId,
            markers1SortedByKmerId,
            maxSkip, maxVertexCountPerKmer, debug, graph, alignment);

        // Compute the AlignmentInfo.
        const AlignmentInfo alignmentInfo(alignment);
        uint32_t leftTrim;
        uint32_t rightTrim;
        tie(leftTrim, rightTrim) = computeTrim(
            orientedReadId0,
            orientedReadId1,
            alignmentInfo);

        cout << orientedReadId0 << " " << orientedReadId1 << " " << alignmentInfo.markerCount;
        if(alignmentInfo.markerCount) {
            cout << " " << leftTrim << " " << rightTrim;
            if(alignmentInfo.markerCount >= minAlignedMarkerCount && leftTrim<=maxTrim && rightTrim<=maxTrim) {
                cout << " good";
                ++goodAlignmentCount;
            }
        }
        cout << endl;

    }
    cout << "Found " << goodAlignmentCount << " alignments out of ";
    cout << alignmentTable[orientedReadId0.getValue()].size() << "." << endl;

}



// Given two oriented reads and their computed AlignmentInfo,
// compute the left and right trim, expressed in markers.
// This is the minimum number of markers (over the two reads)
// that are excluded from the alignment on each side.
// If the trim is too high, the alignment is suspicious.
pair<uint32_t, uint32_t> Assembler::computeTrim(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const AlignmentInfo& alignmentInfo) const
{

    const auto markers0 = markers[orientedReadId0.getValue()];
    const auto markers1 = markers[orientedReadId1.getValue()];
    const uint32_t markerCount0 = uint32_t(markers0.size());
    const uint32_t markerCount1 = uint32_t(markers1.size());
    return alignmentInfo.computeTrim(markerCount0, markerCount1);
}



#if 0
// Compute an alignment for each alignment candidate, but  only store the AlignmentInfo.
// Optionally, the alignments are used for creation of the global marker graph.
void Assembler::computeAllAlignments(

    // The  maximum number of vertices in the alignment graph
    // that we allow a single k-mer to generate.
    size_t alignmentMaxVertexCountPerKmer,

    // The maximum ordinal skip to be tolerated between successive markers
    // in the alignment.
    size_t maxSkip,

    // Minimum number of alignment markers for an alignment to be used.
    size_t minAlignedMarkerCount,

    // Maximum left/right trim (in bases) for an alignment to be used.
    size_t maxTrim,

    // Minimum coverage (number of markers) for a vertex
    // of the marker graph to be kept.
    size_t minCoverage,

    // Number of threads. If zero, a number of threads equal to
    // the number of virtual processors is used.
    size_t threadCount
)
{
    using VertexId = GlobalMarkerGraphVertexId;
    using CompressedVertexId = CompressedGlobalMarkerGraphVertexId;

    const auto tBegin = steady_clock::now();
    cout << timestamp << "Begin computing alignments for ";
    cout << alignmentCandidates.size() << " alignment candidates." << endl;

    // Check that we have what we need.
    checkReadsAreOpen();
    checkKmersAreOpen();
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    // Store parameters so they are accessible to the threads.
    auto& data = computeAllAlignmentsData;
    data.maxSkip = maxSkip;
    data.maxVertexCountPerKmer = alignmentMaxVertexCountPerKmer;
    data.minAlignedMarkerCount = minAlignedMarkerCount;
    data.maxTrim = maxTrim;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;


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



    // Compute the alignments and update the disjoint set data structure
    // for each good alignment.
    data.threadAlignmentData.resize(threadCount);
    cout << timestamp << "Alignment and disjoint set computation begins." << endl;
    size_t batchSize = 10000;
    setupLoadBalancing(alignmentCandidates.size(), batchSize);
    runThreads(&Assembler::computeAllAlignmentsThreadFunction1,
        threadCount, "threadLogs/computeAllAlignments1");
    cout << timestamp << "Alignment and disjoint set computation completed." << endl;



    // Store alignmentInfos found by each thread in the global alignmentInfos.
    cout << "Storing the alignment info objects." << endl;
    alignmentData.createNew(largeDataName("AlignmentData"), largeDataPageSize);
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        const vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];
        for(const AlignmentData& ad: threadAlignmentData) {
            alignmentData.push_back(ad);
        }
    }
    cout << timestamp << "Creating alignment table." << endl;
    computeAlignmentTable();



    // Use the disjoint sets to find the vertex of the global
    // marker graph that each oriented marker belongs to.

    // Compute the global marker graph vertex corresponding
    // to each MarkerId.
    cout << timestamp << "Storing the global marker graph vertex each marker belongs to." << endl;
    globalMarkerGraphVertex.createNew(
        largeDataName("GlobalMarkerGraphVertex"),
        largeDataPageSize);
    globalMarkerGraphVertex.reserveAndResize(data.orientedMarkerCount);
    batchSize = 1000000;
    setupLoadBalancing(markers.totalSize(), batchSize);
    runThreads(&Assembler::computeAllAlignmentsThreadFunction2, threadCount);

    // Clean up.
    data.disjointSetsPointer = 0;
    data.disjointSetsData.remove();


    // Using the disjoint sets data structure we computed a
    // "raw" vertex id that each marker belongs to.
    // These "raw" vertex ids are in [0, number of markers),
    // but are not contiguous, because the number of vertices is
    // less than the number of markers.
    // To compute final vertex ids that are contiguous
    // in [0, numberOfVertices), we begin by counting
    // the number of markers for each raw vertex id.
    // This is stored in workArea.
    cout << timestamp << "Counting the number of markers in each vertex." << endl;
    data.workArea.createNew(
        largeDataName("tmp-ComputeAllAlignmentsWorkArea"),
        largeDataPageSize);
    data.workArea.reserveAndResize(data.orientedMarkerCount);
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::computeAllAlignmentsThreadFunction3, threadCount);



    // Now we can loop over the work area, and replace each
    // entry at least equal to minCoverage with an increasing id
    // which will be the final vertex id.
    // Entries that are 0 are set to maxValueMinus1.
    // Entries that are greater than 0 and less than minCoverage
    // are set to maxValue, so markers that don't belong to
    // any vertex receive a vertex id of all one bits.
    const uint64_t maxValue = std::numeric_limits<uint64_t>::max();
    const uint64_t maxValueMinus1 = maxValue - 1ULL;
    cout << timestamp << "Generating contiguous vertex ids." << endl;
    GlobalMarkerGraphVertexId vertexId = 0;
    for(GlobalMarkerGraphVertexId& w: data.workArea) {
        if(w == 0ULL) {
            w = maxValueMinus1;
        } else if(w <minCoverage) {
            w = maxValue;
        } else {
            w = vertexId++;
        }
    }
    const GlobalMarkerGraphVertexId vertexCount = vertexId;
    cout << timestamp << "The global marker graph has " << vertexCount;
    cout << " vertices for " << data.orientedMarkerCount;
    cout << " oriented markers." << endl;

    // Now we can use the workArea to convert raw vertex ids to final vertex ids.
    cout << timestamp << "Storing final vertex ids for all markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::computeAllAlignmentsThreadFunction4, threadCount);
    data.workArea.remove();



    // Gather the oriented marker ids of each vertex of the global marker graph.
    cout << timestamp << "Gathering the oriented markers of each vertex." << endl;
    globalMarkerGraphVertices.createNew(
        largeDataName("GlobalMarkerGraphVertices"),
        largeDataPageSize);
    cout << timestamp << "... pass 1" << endl;
    globalMarkerGraphVertices.beginPass1(vertexCount);
    const CompressedVertexId maxValue40Bits = maxValue;
    for(const CompressedVertexId& v: globalMarkerGraphVertex) {
        if(v != maxValue40Bits) {
            globalMarkerGraphVertices.incrementCount(v);
        }
    }
    cout << timestamp << "... pass 2" << endl;
    globalMarkerGraphVertices.beginPass2();
    for(VertexId i=0; i<globalMarkerGraphVertex.size(); i++) {
        const CompressedVertexId v = globalMarkerGraphVertex[i];
        if(v != maxValue40Bits) {
            globalMarkerGraphVertices.store(globalMarkerGraphVertex[i], i);
        }
    }
    globalMarkerGraphVertices.endPass2();

    // Sort the markers in each vertex.
    cout << timestamp << "Sorting the oriented markers of each vertex." << endl;
    for(VertexId i=0; i<vertexCount; i++) {
        if(globalMarkerGraphVertices.size(i) > 1) {
            sort(globalMarkerGraphVertices.begin(i), globalMarkerGraphVertices.end(i));
        }
    }


    // Check that all the markers of a vertex have the same kmer id.
    if(true) {
        for(VertexId i=0; i<globalMarkerGraphVertices.size(); i++) {
            const MemoryAsContainer<MarkerId> vertexMarkers = globalMarkerGraphVertices[i];
            KmerId kmerId = 0;
            for(size_t j=0; j<vertexMarkers.size(); j++) {
                const MarkerId markerId = vertexMarkers[j];
                const CompressedMarker& marker = markers.begin()[markerId];
                if(j == 0) {
                    kmerId = marker.kmerId;
                } else {
                    CZI_ASSERT(kmerId == marker.kmerId);
                }
            }
        }
    }

    // Create a histogram of number of markers per vertex.
    cout << "Creating marker graph vertex coverage histogram." << endl;
    vector<uint64_t> histogram;
    for(size_t i=0; i<globalMarkerGraphVertices.size(); i++) {
        const size_t count = globalMarkerGraphVertices.size(i);
        if(histogram.size()<= count) {
            histogram.resize(count+1, 0ULL);
        }
        ++histogram[count];
    }
    ofstream csv("GlobalMarkerGraphHistogram.csv");
    csv << "MarkerCount,Frequency\n";
    for(size_t i=0; i<histogram.size(); i++) {
        const auto n = histogram[i];
        if(n) {
            csv << i << "," << n << "\n";
        }
    }


    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "Computation of alignments and global marker graph ";
    cout << "completed in " << tTotal << " s." << endl;
}



void Assembler::computeAllAlignmentsThreadFunction1(size_t threadId)
{
    ostream& out = getLog(threadId);

    array<OrientedReadId, 2> orientedReadIds;
    array<OrientedReadId, 2> orientedReadIdsOppositeStrand;
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    AlignmentGraph graph;
    Alignment alignment;

    const bool debug = false;
    auto& data = computeAllAlignmentsData;
    const size_t maxSkip = data.maxSkip;
    const size_t maxVertexCountPerKmer = data.maxVertexCountPerKmer;
    const size_t minAlignedMarkerCount = data.minAlignedMarkerCount;
    const size_t maxTrim = data.maxTrim;

    const std::shared_ptr<DisjointSets> disjointSetsPointer = data.disjointSetsPointer;
    vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << "Working on batch " << begin << " " << end << endl;

        for(size_t i=begin; i!=end; i++) {
            const OrientedReadPair& candidate = alignmentCandidates[i];
            CZI_ASSERT(candidate.readIds[0] < candidate.readIds[1]);

            // Get the oriented read ids, with the first one on strand 0.
            orientedReadIds[0] = OrientedReadId(candidate.readIds[0], 0);
            orientedReadIds[1] = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);

            // Get the oriented read ids for the opposite strand.
            orientedReadIdsOppositeStrand = orientedReadIds;
            orientedReadIdsOppositeStrand[0].flipStrand();
            orientedReadIdsOppositeStrand[1].flipStrand();


            // out << timestamp << "Working on " << i << " " << orientedReadIds[0] << " " << orientedReadIds[1] << endl;

            // Get the markers for the two oriented reads in this candidate.
            for(size_t j=0; j<2; j++) {
                getMarkersSortedByKmerId(orientedReadIds[j], markersSortedByKmerId[j]);
            }

            // Compute the Alignment.
            alignOrientedReads(
                markersSortedByKmerId[0],
                markersSortedByKmerId[1],
                maxSkip, maxVertexCountPerKmer, debug, graph, alignment);

            // If the alignment has too few markers skip it.
            if(alignment.ordinals.size() < minAlignedMarkerCount) {
                continue;
            }

            // Compute the AlignmentInfo.
            AlignmentInfo alignmentInfo;
            alignmentInfo.create(alignment);

            // If the alignment has too much trim, skip it.
            uint32_t leftTrim;
            uint32_t rightTrim;
            tie(leftTrim, rightTrim) = computeTrim(
                orientedReadIds[0],
                orientedReadIds[1],
                alignmentInfo);
            if(leftTrim>maxTrim || rightTrim>maxTrim) {
                continue;
            }

            // If getting here, this is a good alignment.
            threadAlignmentData.push_back(AlignmentData(candidate, alignmentInfo));

            // In the global marker graph, merge pairs
            // of aligned markers.
            for(const auto& p: alignment.ordinals) {
                const uint32_t ordinal0 = p.first;
                const uint32_t ordinal1 = p.second;
                const MarkerId markerId0 = getMarkerId(orientedReadIds[0], ordinal0);
                const MarkerId markerId1 = getMarkerId(orientedReadIds[1], ordinal1);
                CZI_ASSERT(markers.begin()[markerId0].kmerId == markers.begin()[markerId1].kmerId);
                disjointSetsPointer->unite(markerId0, markerId1);

                // Also do it for the corresponding oriented markers on
                // the opposite strand.
                const uint32_t ordinal0OppositeStrand = uint32_t(markersSortedByKmerId[0].size()) - 1 - ordinal0;
                const uint32_t ordinal1OppositeStrand = uint32_t(markersSortedByKmerId[1].size()) - 1 - ordinal1;
                const MarkerId markerId0OppositeStrand =
                    getMarkerId(orientedReadIdsOppositeStrand[0], ordinal0OppositeStrand);
                const MarkerId markerId1OppositeStrand =
                    getMarkerId(orientedReadIdsOppositeStrand[1], ordinal1OppositeStrand);
                CZI_ASSERT(markers.begin()[markerId0OppositeStrand].kmerId == markers.begin()[markerId1OppositeStrand].kmerId);
                disjointSetsPointer->unite(
                    markerId0OppositeStrand,
                    markerId1OppositeStrand);
            }
        }
    }

}



void Assembler::computeAllAlignmentsThreadFunction2(size_t threadId)
{
    DisjointSets& disjointSets = *computeAllAlignmentsData.disjointSetsPointer;

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            globalMarkerGraphVertex[i] = disjointSets.find(i);
        }
    }
}



void Assembler::computeAllAlignmentsThreadFunction3(size_t threadId)
{
    GlobalMarkerGraphVertexId* workArea =
        computeAllAlignmentsData.workArea.begin();

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            const MarkerId rawVertexId = globalMarkerGraphVertex[i];
            __sync_fetch_and_add(workArea + rawVertexId, 1ULL);
        }
    }
}



void Assembler::computeAllAlignmentsThreadFunction4(size_t threadId)
{
    GlobalMarkerGraphVertexId* workArea =
        computeAllAlignmentsData.workArea.begin();
    const uint64_t maxValue = std::numeric_limits<uint64_t>::max();
    const uint64_t maxValueMinus1 = maxValue - 1ULL;

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            const MarkerId rawVertexId = globalMarkerGraphVertex[i];
            CZI_ASSERT(rawVertexId != maxValueMinus1);
            const MarkerId finalVertexId = workArea[rawVertexId];
            globalMarkerGraphVertex[i] = finalVertexId;
        }
    }
}
#endif



// Compute an alignment for each alignment candidate.
// Store summary information for the ones that are good enough,
// without storing details of the alignment.
void Assembler::computeAlignments(

    // The  maximum number of vertices in the alignment graph
    // that we allow a single k-mer to generate.
    size_t alignmentMaxVertexCountPerKmer,

    // The maximum ordinal skip to be tolerated between successive markers
    // in the alignment.
    size_t maxSkip,

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
    cout << alignmentCandidates.size() << " alignment candidates." << endl;

    // Check that we have what we need.
    checkReadsAreOpen();
    checkKmersAreOpen();
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    // Store parameters so they are accessible to the threads.
    auto& data = computeAlignmentsData;
    data.maxSkip = maxSkip;
    data.maxVertexCountPerKmer = alignmentMaxVertexCountPerKmer;
    data.minAlignedMarkerCount = minAlignedMarkerCount;
    data.maxTrim = maxTrim;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Compute the alignments.
    data.threadAlignmentData.resize(threadCount);
    cout << timestamp << "Alignment computation begins." << endl;
    size_t batchSize = 10000;
    setupLoadBalancing(alignmentCandidates.size(), batchSize);
    runThreads(&Assembler::computeAlignmentsThreadFunction,
        threadCount, "threadLogs/computeAlignmentsThreadFunction");
    cout << timestamp << "Alignment computation completed." << endl;



    // Store alignmentInfos found by each thread in the global alignmentInfos.
    cout << "Storing the alignment info objects." << endl;
    alignmentData.createNew(largeDataName("AlignmentData"), largeDataPageSize);
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        const vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];
        for(const AlignmentData& ad: threadAlignmentData) {
            alignmentData.push_back(ad);
        }
    }
    cout << timestamp << "Creating alignment table." << endl;
    computeAlignmentTable();

    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "Computation of alignments ";
    cout << "completed in " << tTotal << " s." << endl;
}



void Assembler::computeAlignmentsThreadFunction(size_t threadId)
{
    ostream& out = getLog(threadId);

    array<OrientedReadId, 2> orientedReadIds;
    array<OrientedReadId, 2> orientedReadIdsOppositeStrand;
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    AlignmentGraph graph;
    Alignment alignment;

    const bool debug = false;
    auto& data = computeAlignmentsData;
    const size_t maxSkip = data.maxSkip;
    const size_t maxVertexCountPerKmer = data.maxVertexCountPerKmer;
    const size_t minAlignedMarkerCount = data.minAlignedMarkerCount;
    const size_t maxTrim = data.maxTrim;

    vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << "Working on batch " << begin << " " << end << endl;

        for(size_t i=begin; i!=end; i++) {
            const OrientedReadPair& candidate = alignmentCandidates[i];
            CZI_ASSERT(candidate.readIds[0] < candidate.readIds[1]);

            // Get the oriented read ids, with the first one on strand 0.
            orientedReadIds[0] = OrientedReadId(candidate.readIds[0], 0);
            orientedReadIds[1] = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);

            // Get the oriented read ids for the opposite strand.
            orientedReadIdsOppositeStrand = orientedReadIds;
            orientedReadIdsOppositeStrand[0].flipStrand();
            orientedReadIdsOppositeStrand[1].flipStrand();


            // out << timestamp << "Working on " << i << " " << orientedReadIds[0] << " " << orientedReadIds[1] << endl;

            // Get the markers for the two oriented reads in this candidate.
            for(size_t j=0; j<2; j++) {
                getMarkersSortedByKmerId(orientedReadIds[j], markersSortedByKmerId[j]);
            }

            // Compute the Alignment.
            alignOrientedReads(
                markersSortedByKmerId[0],
                markersSortedByKmerId[1],
                maxSkip, maxVertexCountPerKmer, debug, graph, alignment);

            // If the alignment has too few markers skip it.
            if(alignment.ordinals.size() < minAlignedMarkerCount) {
                continue;
            }

            // Compute the AlignmentInfo.
            AlignmentInfo alignmentInfo;
            alignmentInfo.create(alignment);

            // If the alignment has too much trim, skip it.
            uint32_t leftTrim;
            uint32_t rightTrim;
            tie(leftTrim, rightTrim) = computeTrim(
                orientedReadIds[0],
                orientedReadIds[1],
                alignmentInfo);
            if(leftTrim>maxTrim || rightTrim>maxTrim) {
                continue;
            }

            // If getting here, this is a good alignment.
            threadAlignmentData.push_back(AlignmentData(candidate, alignmentInfo));
        }
    }
}



// Compute alignmentTable from alignmentData.
// This could be made multithreaded if it becomes a bottleneck.
void Assembler::computeAlignmentTable()
{
    alignmentTable.createNew(largeDataName("AlignmentTable"), largeDataPageSize);
    alignmentTable.beginPass1(ReadId(2 * reads.size()));
    for(const AlignmentData& ad: alignmentData) {
        const auto& readIds = ad.readIds;
        OrientedReadId orientedReadId0(readIds[0], 0);
        OrientedReadId orientedReadId1(readIds[1], ad.isSameStrand ? 0 : 1);
        alignmentTable.incrementCount(orientedReadId0.getValue());
        alignmentTable.incrementCount(orientedReadId1.getValue());
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        alignmentTable.incrementCount(orientedReadId0.getValue());
        alignmentTable.incrementCount(orientedReadId1.getValue());
    }
    alignmentTable.beginPass2();
    for(uint32_t i=0; i<alignmentData.size(); i++) {
        const AlignmentData& ad = alignmentData[i];
        const auto& readIds = ad.readIds;
        OrientedReadId orientedReadId0(readIds[0], 0);
        OrientedReadId orientedReadId1(readIds[1], ad.isSameStrand ? 0 : 1);
        alignmentTable.store(orientedReadId0.getValue(), i);
        alignmentTable.store(orientedReadId1.getValue(), i);
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        alignmentTable.store(orientedReadId0.getValue(), i);
        alignmentTable.store(orientedReadId1.getValue(), i);
    }
    alignmentTable.endPass2();



    // Sort each section of the alignment table by OrientedReadId.
    vector< pair<OrientedReadId, uint32_t> > v;
    for(ReadId readId0=0; readId0<reads.size(); readId0++) {
        for(Strand strand0=0; strand0<2; strand0++) {
            const OrientedReadId orientedReadId0(readId0, strand0);

            // Access the section of the alignment table for this oriented read.
            const MemoryAsContainer<uint32_t> alignmentTableSection =
                alignmentTable[orientedReadId0.getValue()];

            // Store pairs(OrientedReadId, alignmentIndex).
            v.clear();
            for(uint32_t alignmentIndex: alignmentTableSection) {
                const AlignmentData& alignment = alignmentData[alignmentIndex];
                const OrientedReadId orientedReadId1 = alignment.getOther(orientedReadId0);
                v.push_back(make_pair(orientedReadId1, alignmentIndex));
            }

            // Sort.
            sort(v.begin(), v.end());

            // Store the sorted alignmentIndex.
            for(size_t i=0; i<v.size(); i++) {
                alignmentTableSection[i] = v[i].second;
            }
        }
    }


}



void Assembler::accessAlignmentData()
{
    alignmentData.accessExistingReadOnly(largeDataName("AlignmentData"));
    alignmentTable.accessExistingReadOnly(largeDataName("AlignmentTable"));
}



void Assembler::checkAlignmentDataAreOpen()
{
    if(!alignmentData.isOpen || !alignmentTable.isOpen()) {
        throw runtime_error("Alignment data are not accessible.");
    }
}



// Find in the alignment table the alignments involving
// a given oriented read, and return them with the correct
// orientation (this may involve a swap and/or reverse complement
// of the AlignmentInfo stored in the alignmentTable).
vector< pair<OrientedReadId, AlignmentInfo> >
    Assembler::findOrientedAlignments(OrientedReadId orientedReadId0Argument) const
{
    const ReadId readId0 = orientedReadId0Argument.getReadId();
    const ReadId strand0 = orientedReadId0Argument.getStrand();

    vector< pair<OrientedReadId, AlignmentInfo> > result;

    // Loop over alignment involving this read, as stored in the
    // alignment table.
    const auto alignmentTable0 = alignmentTable[orientedReadId0Argument.getValue()];
    for(const auto i: alignmentTable0) {
        const AlignmentData& ad = alignmentData[i];

        // Get the oriented read ids that the AlignmentData refers to.
        OrientedReadId orientedReadId0(ad.readIds[0], 0);
        OrientedReadId orientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
        AlignmentInfo alignmentInfo = ad.info;

        // Swap oriented reads, if necessary.
        if(orientedReadId0.getReadId() != readId0) {
            swap(orientedReadId0, orientedReadId1);
            alignmentInfo.swap();
        }
        CZI_ASSERT(orientedReadId0.getReadId() == readId0);

        // Get the number of markers for each of the two reads.
        // We will need it below.
        const uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());
        const uint32_t markerCount1 = uint32_t(markers[orientedReadId1.getValue()].size());

        // Reverse complement, if necessary.
        if(orientedReadId0.getStrand() != strand0) {
            orientedReadId0.flipStrand();
            orientedReadId1.flipStrand();
            alignmentInfo.reverseComplement(markerCount0, markerCount1);
        }
        CZI_ASSERT(orientedReadId0.getStrand() == strand0);
        CZI_ASSERT(orientedReadId0 == orientedReadId0Argument);

        result.push_back(make_pair(orientedReadId1, alignmentInfo));
    }
    return result;
}



// Given the oriented alignments computed by findOrientedAlignments,
// compute a vector of alignment coverage for the given oriented read.
// This is a vector that contains, for each marker of the given
// oriented read, the number of alignments "covering" that marker.
// An alignment "covers" as marker if it starts at or before
// the marker and it ends at or after the marker.
// It does not matter whether the marker is in the alignment.
// This is used to detect portions of the read that may contain
// unreliable alignments due to repeats (most commonly for a human genome,
// LINE repeats).
void Assembler::computeAlignmentCoverage(
    OrientedReadId orientedReadId0,
    const vector< pair<OrientedReadId, AlignmentInfo> >& alignments,   // Computed by findOrientedAlignments
    vector<uint32_t>& coverage
    ) const
{
    // Begin with zero coverage.
    const uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());
    coverage.clear();
    coverage.resize(markerCount0, 0);

    // Loop over the alignments.
    for(const auto& p: alignments) {
        const AlignmentInfo& alignment = p.second;
        const auto begin = alignment.firstOrdinals.first;
        const auto end = alignment.lastOrdinals.first;
        for(size_t ordinal=begin; ordinal<=end; ordinal++) {
            (coverage[ordinal])++;
        }
    }


    // For debugging, write to csv file.
    const bool debug = false;
    if(debug) {
        ofstream csv("AlignmentCoverage.csv");
        csv << "Ordinal,Coverage\n";
        for(uint32_t ordinal=0; ordinal<markerCount0; ordinal++) {
            csv << ordinal << "," << coverage[ordinal] << "\n";
        }
    }
}



// Find out if an alignment stored in alignmentData is a containment alignment.
// The alignment is a containment alignment if one (or both)
// of the two reads are entirely contained in the alignment,
// except possibly for up to maxTrim markers at each end.
// The alignment to be checked is specified by its alignmentId, an index
// into the alignmentData vector.
bool Assembler::isContainmentAlignment(
    uint64_t alignmentId,
    size_t maxTrim) const
{
    // Access the specified alignment.
    const AlignmentData& alignment = alignmentData[alignmentId];

    // Get the oriented reads for the alignment as stored
    // (it is stored with the first read on strand 0).
    const OrientedReadId orientedReadId0(alignment.readIds[0], 0);
    const OrientedReadId orientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);

    // Get the number of markers for each of the oriented reads.
    const uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());
    const uint32_t markerCount1 = uint32_t(markers[orientedReadId1.getValue()].size());

    // Compute the number of markers excluded from the alignment on each side.
    const uint32_t leftUnaligned0 = alignment.info.firstOrdinals.first;
    const uint32_t leftUnaligned1 = alignment.info.firstOrdinals.second;
    const uint32_t rightUnaligned0 = markerCount0 - alignment.info.lastOrdinals.first - 1;
    const uint32_t rightUnaligned1 = markerCount1 - alignment.info.lastOrdinals.second - 1;

    // Check for containment.
    if(leftUnaligned0<=maxTrim && rightUnaligned0<=maxTrim) {
        return true;    // The alignment covers the entire oriented read 0.
    }
    if(leftUnaligned1<=maxTrim && rightUnaligned1<=maxTrim) {
        return true;    // The alignment covers the entire oriented read 1.
    }
    return false;   // Neither of the above was the case.
}


