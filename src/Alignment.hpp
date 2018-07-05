#ifndef CZI_SHASTA_ALIGNMENT_HPP
#define CZI_SHASTA_ALIGNMENT_HPP

#include "ReadId.hpp"

#include "array.hpp"
#include "cstdint.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {

        class Alignment;
        class AlignmentData;
        class AlignmentInfo;
    }
}



class ChanZuckerberg::shasta::Alignment {
public:

    // The ordinals in each of the two oriented reads of the
    // markers in the alignment.
    vector< pair<uint32_t, uint32_t> > ordinals;
};



class ChanZuckerberg::shasta::AlignmentInfo {
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



class ChanZuckerberg::shasta::AlignmentData {
public:

    // The ids of the overlapping oriented reads.
    // They are guaranteed to be distinct.
    // That is, overlapping with self is not stored,
    // even on opposite strands.
    array<ReadId, 2> readIds;

    // Flag that is true if the alignment is obtained with the
    // two reads on the same strand.
    bool isSameStrand;

    // The AlignmentInfo computed with the first read on strand 0.
    AlignmentInfo info;

    AlignmentData() {}
    AlignmentData(
        const array<ReadId, 2>& readIds,
        bool isSameStrand,
        const AlignmentInfo& info) :
        readIds(readIds),
        isSameStrand(isSameStrand),
        info(info)
    {}

    // Given one of the oriented read ids involved in this alignment,
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
