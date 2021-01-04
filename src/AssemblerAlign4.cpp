#include "Assembler.hpp"
#include "Align4.hpp"
#include "MemoryMappedAllocator.hpp"
#include "html.hpp"
using namespace shasta;


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

    MemoryMapped::ByteAllocator byteAllocator(
        largeDataName("tmp-ByteAllocator"), largeDataPageSize, 1024 * 1024 * 1024);

    Alignment alignment;
    AlignmentInfo alignmentInfo;

    const bool debug = true;

    // Compute the alignment.
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
    const auto markers0 = markers[orientedReadId0.getValue()];
    const auto markers1 = markers[orientedReadId1.getValue()];

    align4(markers0, markers1,
        options, byteAllocator, alignment, alignmentInfo, debug);
}





