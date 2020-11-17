#ifndef SHASTA_ALIGNMENT_GRAPH4_HPP
#define SHASTA_ALIGNMENT_GRAPH4_HPP

// Class AlignmentGraph4 is used to implement alignment method 4.

#include "Marker.hpp"
#include "span.hpp"

namespace shasta {
    class AlignmentGraph4;
    class Alignment;
    class AlignmentInfo;
}


class shasta::AlignmentGraph4 {
public:

    class Options {
    public:
        uint64_t m;
        uint64_t maxSkip;
        uint64_t maxDrift;
        int64_t matchScore;
        int64_t mismatchScore;
        int64_t gapScore;
    };

    // Align two arbitrary sequences  using alignment method 4.
    // If debug is true, detailed output to html is produced.
    // Otherwise, html is not used.
    static void align(
        const span<const CompressedMarker>&,
        const span<const CompressedMarker>&,
        const AlignmentGraph4::Options&,
        Alignment&,
        AlignmentInfo&,
        bool debug,
        ostream& html);
};



#endif
