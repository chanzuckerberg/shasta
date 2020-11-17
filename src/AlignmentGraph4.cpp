#include "AlignmentGraph4.hpp"
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
    if(debug) {
        html << "Computing a marker alignment of two sequences with " <<
            markers0.size() << " and " << markers1.size() << " markers." << endl;
    }
    SHASTA_ASSERT(0);
}
