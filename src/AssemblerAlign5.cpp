#include "Assembler.hpp"
#include "Align5.hpp"
#include "html.hpp"
using namespace shasta;


// Python-callable version.
void Assembler::alignOrientedReads5(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    uint64_t m,
    uint64_t deltaX,
    uint64_t deltaY,
    int64_t matchScore,
    int64_t mismatchScore,
    int64_t gapScore) const
{
    Align5::Options options;
    options.m = m;
    options.deltaX = deltaX;
    options.deltaY = deltaY;
    options.matchScore = matchScore;
    options.mismatchScore = mismatchScore;
    options.gapScore = gapScore;

    Alignment alignment;
    AlignmentInfo alignmentInfo;

    const bool debug = true;

    // Compute the alignment.
    alignOrientedReads5(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1),
        options, alignment, alignmentInfo, debug);
}



// Align two reads using alignment method 5.
void Assembler::alignOrientedReads5(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const Align5::Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug) const
{
    const auto markers0 = markers[orientedReadId0.getValue()];
    const auto markers1 = markers[orientedReadId1.getValue()];

    align5(markers0, markers1,
        options, alignment, alignmentInfo, debug);
}





