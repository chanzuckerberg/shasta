#include "Assembler.hpp"
using namespace shasta;


// Python-callable version.
void Assembler::alignOrientedReads4(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    uint64_t m,
    uint64_t maxSkip,
    uint64_t maxDrift,
    int64_t matchScore,
    int64_t mismatchScore,
    int64_t gapScore) const
{
    AlignmentGraph4Options options;
    options.m = m;
    options.maxSkip = maxSkip;
    options.maxDrift = maxDrift;
    options.matchScore = matchScore;
    options.mismatchScore = mismatchScore;
    options.gapScore = gapScore;

    Alignment alignment;
    AlignmentInfo alignmentInfo;

    const bool debug = true;
    ofstream html("Align4.html");

    alignOrientedReads4(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1),
        options, alignment, alignmentInfo, debug, html);
}



// Align two reads using alignment method 4.
// If debug is true, detailed output to html is produced.
// Otherwise, html is not used.
void Assembler::alignOrientedReads4(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const AlignmentGraph4Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug,
    ostream& html) const
{
    const auto markers0 = markers[orientedReadId0.getValue()];
    const auto markers1 = markers[orientedReadId1.getValue()];

    align4(markers0, markers1,
        options, alignment, alignmentInfo, debug, html);
}





