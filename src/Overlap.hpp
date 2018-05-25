#ifndef CZI_NANOPORE2_OVERLAP_HPP
#define CZI_NANOPORE2_OVERLAP_HPP

#include "ReadId.hpp"

#include "array.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {
        class Overlap;
    }
}


// Class describing the overlap between a pair of oriented reads.
class ChanZuckerberg::Nanopore2::Overlap {
public:

    // The id of the overlapping oriented reads.
    array<OrientedReadId, 2> orientedReadIds;

    // The number of minHash iteration that found this pair.
    uint32_t minHashFrequency;

    Overlap() : minHashFrequency(0) {}
    Overlap(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        uint32_t minHashFrequency
    ) :
    orientedReadIds(array<OrientedReadId, 2>{orientedReadId0, orientedReadId1}),
    minHashFrequency(minHashFrequency)
    {}
};

#endif
