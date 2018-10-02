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
    cout << " bases, right trim " << rightTrim << " bases." << endl;

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
    size_t maxTrim                  // Maximum trim allowed in an alignment.
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
    size_t maxTrim                  // Maximum trim allowed in an alignment.
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
// compute the left and right trim.
// This is the minimum number of bases (over the two reads)
// that are excluded from the alignment on each size.
// This takes as input markers sorted by position.
pair<uint32_t, uint32_t> Assembler::computeTrim(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const AlignmentInfo& alignmentInfo)
{
    const uint32_t baseCount0 = uint32_t(reads[orientedReadId0.getReadId()].baseCount);
    const uint32_t baseCount1 = uint32_t(reads[orientedReadId1.getReadId()].baseCount);

    const auto markers0 = markers[orientedReadId0.getValue()];
    const auto markers1 = markers[orientedReadId1.getValue()];

    if(alignmentInfo.markerCount) {
        const CompressedMarker& firstMarker0 = markers0[alignmentInfo.firstOrdinals.first];
        const CompressedMarker& firstMarker1 = markers1[alignmentInfo.firstOrdinals.second];
        const CompressedMarker& lastMarker0 = markers0[alignmentInfo.lastOrdinals.first];
        const CompressedMarker& lastMarker1 = markers1[alignmentInfo.lastOrdinals.second];
        return make_pair(
            min(firstMarker0.position, firstMarker1.position),
            min(baseCount0-1-lastMarker0.position, baseCount1-1-lastMarker1.position)
            );
    } else {
        // The alignment is empty.
        // The trim equal the entire lenght of the shortest read.
        const uint32_t trim = min(baseCount0, baseCount1);
        return make_pair(trim, trim);

    }

}



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
    for(uint32_t i=0; i<alignmentTable.size(); i++) {
        sort(alignmentTable.begin(i), alignmentTable.end(i));
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

