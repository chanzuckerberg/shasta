#ifndef CZI_NANOPORE2_ALIGNMENT_HPP
#define CZI_NANOPORE2_ALIGNMENT_HPP

#include "cstdint.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {

        class Alignment;
        class AlignmentInfo;
    }
}



class ChanZuckerberg::Nanopore2::Alignment {
public:

    // The ordinals in each of the two oriented reads of the
    // markers in the alignment.
    vector< pair<uint32_t, uint32_t> > ordinals;
};



class ChanZuckerberg::Nanopore2::AlignmentInfo {
public:

    // The first and last ordinals in each of the two oriented reads.
    // These are not filled in if markerCount is zero.
    pair<uint32_t, uint32_t> firstOrdinals;
    pair<uint32_t, uint32_t> lastOrdinals;

    // The number of markers in the alignment.
    uint32_t markerCount;

    AlignmentInfo(const Alignment& alignment)
    {
        create(alignment);
    }
    void create(const Alignment& alignment)
    {
        markerCount = uint32_t(alignment.ordinals.size());
        if(markerCount) {
            firstOrdinals = alignment.ordinals.front();
            lastOrdinals  = alignment.ordinals.back();
        }
    }
    AlignmentInfo() {}

};
#endif
