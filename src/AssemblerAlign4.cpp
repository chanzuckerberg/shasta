#include "Assembler.hpp"
#include "html.hpp"
using namespace shasta;


// Python-callable version.
void Assembler::alignOrientedReads4(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    uint64_t m,
    uint64_t deltaX,
    uint64_t deltaY,
    int64_t matchScore,
    int64_t mismatchScore,
    int64_t gapScore) const
{
    Align4Options options;
    options.m = m;
    options.deltaX = deltaX;
    options.deltaY = deltaY;
    options.matchScore = matchScore;
    options.mismatchScore = mismatchScore;
    options.gapScore = gapScore;

    Alignment alignment;
    AlignmentInfo alignmentInfo;

    // Debug output in html format.
    const bool debug = true;
    ofstream html("Align4.html");
    shasta::writeHtmlBegin(html, "Align4");
    writeMakeAllTablesCopyable(html);
    html << "<body onload='makeAllTablesCopyable()'>\n"
        "<h1>Alignment method 4 for oriented reads " <<
        OrientedReadId(readId0, strand0) << " and " <<
        OrientedReadId(readId1, strand1) << "</h1>\n";

    // Compute the alignment.
    alignOrientedReads4(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1),
        options, alignment, alignmentInfo, debug, html);

    // Finish the html.
    html << "</body>\n";
    writeHtmlEnd(html);
}



// Align two reads using alignment method 4.
// If debug is true, detailed output to html is produced.
// Otherwise, html is not used.
void Assembler::alignOrientedReads4(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const Align4Options& options,
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





