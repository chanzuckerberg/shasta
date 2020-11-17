#include "AlignmentGraph4.hpp"
#include "html.hpp"
using namespace shasta;

// Align two arbitrary sequences  using alignment method 4.
// If debug is true, detailed output to html is produced.
// Otherwise, html is not used.
void AlignmentGraph4::align(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const AlignmentGraph4::Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug,
    ostream& html)
{
    switch(options.m) {
    case 1:
        align<1>(markers0, markers1, options, alignment, alignmentInfo, debug, html);
        return;
    case 2:
        align<2>(markers0, markers1, options, alignment, alignmentInfo, debug, html);
        return;
    case 3:
        align<3>(markers0, markers1, options, alignment, alignmentInfo, debug, html);
        return;
    case 4:
        align<4>(markers0, markers1, options, alignment, alignmentInfo, debug, html);
        return;
    default:
        SHASTA_ASSERT(0);
    }
}



// Version templated on m, the number of markers that define
// a "feature" used in the alignment.
template<uint64_t m> void AlignmentGraph4::align(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const AlignmentGraph4::Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug,
    ostream& html)
{
    // Check that we are in the templated version consistent with
    /// the options.
    SHASTA_ASSERT(options.m == m);

    if(debug) {
        writeHtmlBegin(html, "");
        writeStyle(html);
        html << "<body>";
        html << "<pre>";
        html << "Computing a marker alignment of two sequences with " <<
            markers0.size() << " and " << markers1.size() << " markers." << endl;
    }


    if(debug) {
        html << "</pre>";
        html << "</body>";
        writeHtmlEnd(html);
    }
}
