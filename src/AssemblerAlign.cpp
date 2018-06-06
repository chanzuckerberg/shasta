// Nanopore2.
#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Standard libraries.
#include "chrono.hpp"
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

    // Get the markers sorted by position.
    vector<Marker> markers0SortedByPosition;
    vector<Marker> markers1SortedByPosition;
    getMarkers(orientedReadId0, markers0SortedByPosition);
    getMarkers(orientedReadId1, markers1SortedByPosition);

    // Get the markers sorted by kmerId.
    vector<Marker> markers0SortedByKmerId = markers0SortedByPosition;
    vector<Marker> markers1SortedByKmerId = markers1SortedByPosition;
    sort(markers0SortedByKmerId.begin(), markers0SortedByKmerId.end(), OrderMarkersByKmerId());
    sort(markers1SortedByKmerId.begin(), markers1SortedByKmerId.end(), OrderMarkersByKmerId());

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
        markers0SortedByPosition,
        markers1SortedByPosition,
        alignmentInfo);
    cout << orientedReadId0 << " has " << reads[orientedReadId0.getReadId()].baseCount;
    cout << " bases and " << markers0SortedByPosition.size() << " markers." << endl;
    cout << orientedReadId1 << " has " << reads[orientedReadId1.getReadId()].baseCount;
    cout << " bases and " << markers1SortedByPosition.size() << " markers." << endl;
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
    const vector<Marker>& markers0SortedByKmerId,
    const vector<Marker>& markers1SortedByKmerId,
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
    const vector<Marker>& markers0SortedByKmerId,
    const vector<Marker>& markers1SortedByKmerId,
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
// for which we have an Overlap.
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
    checkOverlapsAreOpen();

    // Get the markers for orientedReadId0.
    vector<Marker> markers0SortedByPosition;
    getMarkers(orientedReadId0, markers0SortedByPosition);
    vector<Marker> markers0SortedByKmerId = markers0SortedByPosition;
    sort(markers0SortedByKmerId.begin(), markers0SortedByKmerId.end(), OrderMarkersByKmerId());

    // Loop over all overlaps involving this oriented read.
    vector<Marker> markers1SortedByPosition;
    vector<Marker> markers1SortedByKmerId;
    size_t goodAlignmentCount = 0;
    for(const uint64_t i: overlapTable[orientedReadId0.getValue()]) {
        const Overlap& overlap = overlaps[i];

        // Get the other oriented read involved in this overlap.
        const OrientedReadId orientedReadId1 = overlap.getOther(orientedReadId0);

        // Get the markers for orientedReadId1.
        getMarkers(orientedReadId1, markers1SortedByPosition);
        markers1SortedByKmerId = markers1SortedByPosition;
        sort(markers1SortedByKmerId.begin(), markers1SortedByKmerId.end(), OrderMarkersByKmerId());

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
            markers0SortedByPosition,
            markers1SortedByPosition,
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
    cout << "Found " << goodAlignmentCount << " good alignments among ";
    cout << overlapTable[orientedReadId0.getValue()].size() << " overlaps." << endl;

}



// Given two oriented reads and their computed AlignmentInfo,
// compute the left and right trim.
// This is the minimum number of bases (over the two reads)
// that are excluded from the alignment on each size.
// This takes as input markers sorted by position.
pair<uint32_t, uint32_t> Assembler::computeTrim(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    vector<Marker>& markers0,
    vector<Marker>& markers1,
    const AlignmentInfo& alignmentInfo)
{
    const uint32_t baseCount0 = uint32_t(reads[orientedReadId0.getReadId()].baseCount);
    const uint32_t baseCount1 = uint32_t(reads[orientedReadId1.getReadId()].baseCount);

    if(alignmentInfo.markerCount) {
        const Marker& firstMarker0 = markers0[alignmentInfo.firstOrdinals.first];
        const Marker& firstMarker1 = markers1[alignmentInfo.firstOrdinals.second];
        const Marker& lastMarker0 = markers0[alignmentInfo.lastOrdinals.first];
        const Marker& lastMarker1 = markers1[alignmentInfo.lastOrdinals.second];
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



// Compute an Alignment for each Overlap, but  only store the AlignmentInfo.
// Optionally, the alignments are used for creation of the global marker graph.
void Assembler::computeAllAlignments(

    // Minimum number of MinHash hits for an alignment to be computed.
    size_t minFrequency,

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
    size_t threadCount,

    // Flag to control computation of the global marker graph.
    // This should normally be true.
    bool computeGlobalMarkerGraph
)
{
    const auto tBegin = steady_clock::now();
    cout << timestamp << "Begin computing alignments for ";
    cout << overlaps.size() << " overlaps." << endl;

    // Check that we have what we need.
    checkReadsAreOpen();
    checkKmersAreOpen();
    checkMarkersAreOpen();
    checkOverlapsAreOpen();

    // Store alignment parameters so they are accessible to the threads.
    computeAllAlignmentsData.maxSkip = maxSkip;
    computeAllAlignmentsData.maxVertexCountPerKmer = alignmentMaxVertexCountPerKmer;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Create the alignmentInfos.
    alignmentInfos.createNew(largeDataName("AlignmentInfos"), largeDataPageSize);
    alignmentInfos.resize(overlaps.size());

    // Do the computation.
    const size_t batchSize = 10000;
    setupLoadBalancing(overlaps.size(), batchSize);
    runThreads(&Assembler::computeAllAlignmentsThreadFunction,
        threadCount, "threadLogs/computeAllAlignments");

    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "Computation of alignments completed in " << tTotal << " s." << endl;
}



void Assembler::computeAllAlignmentsThreadFunction(size_t threadId)
{
    ostream& out = getLog(threadId);

    array<OrientedReadId, 2> orientedReadIds;
    array<vector<Marker>, 2> markersSortedByPosition;
    array<vector<Marker>, 2> markersSortedByKmerId;
    AlignmentGraph graph;
    Alignment alignment;
    const bool debug = false;
    const size_t maxSkip = computeAllAlignmentsData.maxSkip;
    const size_t maxVertexCountPerKmer = computeAllAlignmentsData.maxVertexCountPerKmer;

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << "Working on batch " << begin << " " << end << endl;

        for(size_t i=begin; i!=end; i++) {
            const Overlap& overlap = overlaps[i];
            AlignmentInfo& alignmentInfo = alignmentInfos[i];

            orientedReadIds[0] = OrientedReadId(overlap.readIds[0], 0);
            orientedReadIds[1] = OrientedReadId(overlap.readIds[1], overlap.isSameStrand ? 0 : 1);

            // out << timestamp << "Working on " << i << " " << orientedReadIds[0] << " " << orientedReadIds[1] << endl;

            // Get the markers for the two oriented reads in this Overlap.
            for(size_t j=0; j<2; j++) {
                getMarkers(orientedReadIds[j], markersSortedByPosition[j]);
                markersSortedByKmerId[j] = markersSortedByPosition[j];
                sort(markersSortedByKmerId[j].begin(), markersSortedByKmerId[j].end(),
                    OrderMarkersByKmerId());
            }

            // Compute the alignment.
            alignOrientedReads(
                markersSortedByKmerId[0],
                markersSortedByKmerId[1],
                maxSkip, maxVertexCountPerKmer, debug, graph, alignment);

            // Store the AlignmentInfo.
            alignmentInfo.create(alignment);
        }
    }

}



void Assembler::accessAlignmentInfos()
{
    alignmentInfos.accessExistingReadOnly(
        largeDataName("AlignmentInfos"));
}



void Assembler::checkAlignmentInfosAreOpen()
{
    if(!alignmentInfos.isOpen) {
        throw runtime_error("Alignment infos are not accessible.");
    }
}

