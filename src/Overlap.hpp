#ifndef CZI_SHASTA_OVERLAP_HPP
#define CZI_SHASTA_OVERLAP_HPP

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

    // The ids of the overlapping oriented reads.
    // They are guaranteed to be distinct.
    // That is, overlapping with self is not stored,
    // even on opposite strands.
    array<ReadId, 2> readIds;

    // The number of minHash iteration that found this pair.
    uint32_t minHashFrequency;

    // Flag that is true if the alignment is obtained with the
    // two reads on the same strand.
    bool isSameStrand;

    Overlap() : minHashFrequency(0) {}
    Overlap(
        ReadId readId0,
        ReadId readId1,
        bool isSameStrand,
        uint32_t minHashFrequency
    ) :
    readIds(array<ReadId, 2>{readId0, readId1}),
    minHashFrequency(minHashFrequency),
    isSameStrand(isSameStrand)
    {
        CZI_ASSERT(readId0 != readId1);
    }



    // Given one of the oriented read ids involved in this overlap,
    // return the other.
    OrientedReadId getOther(OrientedReadId orientedReadIdA) const
    {
        // Get the read id and strand of the given oriented read.
        const ReadId readIdA = orientedReadIdA.getReadId();
        const Strand strandA = orientedReadIdA.getStrand();

        // Find out if readIdA occurs at position 0 or 1.
        int iA;
        if(readIdA == readIds[0]) {
            iA = 0;
        } else if(readIdA == readIds[1]) {
            iA = 1;
        } else {
            CZI_ASSERT(0);
        }

        // Find the desired read id and strand.
        const int iB = 1 - iA;
        const ReadId readIdB = readIds[iB];
        const Strand strandB = isSameStrand ? strandA : 1-strandA;

        return OrientedReadId(readIdB, strandB);
    }
};

#endif
