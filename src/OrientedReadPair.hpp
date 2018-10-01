#ifndef CZI_SHASTA_ORIENTED_READ_PAIR_HPP
#define CZI_SHASTA_ORIENTED_READ_PAIR_HPP

// Shasta.
#include "ReadId.hpp"

// Standard library.
#include "array.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class OrientedReadPair;
    }
}



// Class describing a pair of reads with a relative orientation
// (same strand or different strands).
class ChanZuckerberg::shasta::OrientedReadPair {
public:

    // The read ids are guaranteed to be distinct.
    array<ReadId, 2> readIds;

    // Flag that is true if the two reads are on the same strand.
    bool isSameStrand;

    // Constructors.
    OrientedReadPair() {}
    OrientedReadPair(
        ReadId readId0,
        ReadId readId1,
        bool isSameStrand
        ) :
        readIds(array<ReadId, 2>{readId0, readId1}),
        isSameStrand(isSameStrand)
    {
        CZI_ASSERT(readId0 != readId1);
    }



    // Given one of the oriented read ids, return the other,
    // taking into account the relative orientation.
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
