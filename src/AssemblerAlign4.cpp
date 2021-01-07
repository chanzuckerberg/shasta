// Shasta.
#include "Assembler.hpp"
#include "Align4.hpp"
#include "MemoryMappedAllocator.hpp"
#include "orderPairs.hpp"
using namespace shasta;

// Standard library.
#include "array.hpp"


// Python-callable version.
void Assembler::alignOrientedReads4(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    uint64_t deltaX,
    uint64_t deltaY,
    uint64_t minEntryCountPerCell,
    uint64_t maxDistanceFromBoundary,
    uint64_t minAlignedMarkerCount,
    double minAlignedFraction,
    uint64_t maxSkip,
    uint64_t maxDrift,
    uint64_t maxTrim,
    int64_t matchScore,
    int64_t mismatchScore,
    int64_t gapScore) const
{
    // Fill in the options.
    Align4::Options options;
    options.deltaX = deltaX;
    options.deltaY = deltaY;
    options.minEntryCountPerCell = minEntryCountPerCell;
    options.maxDistanceFromBoundary = maxDistanceFromBoundary;
    options.minAlignedMarkerCount = minAlignedMarkerCount;
    options.minAlignedFraction = minAlignedFraction;
    options.maxSkip = maxSkip;
    options.maxDrift = maxDrift;
    options.maxTrim = maxTrim;
    options.matchScore = matchScore;
    options.mismatchScore = mismatchScore;
    options.gapScore = gapScore;

    // Set up the memory allocator.
    MemoryMapped::ByteAllocator byteAllocator(
        largeDataName("tmp-ByteAllocator"), largeDataPageSize, 1024 * 1024 * 1024);


    // Compute the alignment.
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    const bool debug = true;
    alignOrientedReads4(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1),
        options, byteAllocator, alignment, alignmentInfo, debug);
    cout << "The alignment has " << alignmentInfo.markerCount << " markers." << endl;
}



// Align two reads using alignment method 4.
void Assembler::alignOrientedReads4(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const Align4::Options& options,
    MemoryMapped::ByteAllocator& byteAllocator,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug) const
{
    // Access the markers for the two oriented reads.
    array<span<const CompressedMarker>, 2> orientedReadMarkers;
    orientedReadMarkers[0] = markers[orientedReadId0.getValue()];
    orientedReadMarkers[1] = markers[orientedReadId1.getValue()];



    // Compute markers sorted by KmerId.
    array<vector< pair<KmerId, uint32_t> >, 2> orientedReadSortedMarkers;
    for(uint64_t i=0; i<2; i++) {

        // Unsorted markers for this oriented read.
        const span<const CompressedMarker>& um = orientedReadMarkers[i];

        // Sorted markers for this oriented read.
        vector<pair<KmerId, uint32_t> >& sm = orientedReadSortedMarkers[i];

        // Copy the unsorted markers.
        const uint64_t n = um.size();
        sm.resize(n);
        for(uint64_t ordinal=0; ordinal<n; ordinal++) {
            const CompressedMarker& cm = um[ordinal];
            sm[ordinal] = make_pair(cm.kmerId, uint32_t(ordinal));
        }

        // Sort them.
        sort(sm.begin(), sm.end(), OrderPairsByFirstOnly<KmerId, uint32_t>());
    }
    array<span< pair<KmerId, uint32_t> >, 2> orientedReadSortedMarkersSpans =
        {orientedReadSortedMarkers[0], orientedReadSortedMarkers[1]};



    // Compute the alignment.
    Align4::align(orientedReadMarkers, orientedReadSortedMarkersSpans,
        options, byteAllocator, alignment, alignmentInfo, debug);
}





