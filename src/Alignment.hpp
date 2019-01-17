#ifndef CZI_SHASTA_ALIGNMENT_HPP
#define CZI_SHASTA_ALIGNMENT_HPP

#include "OrientedReadPair.hpp"
#include "ReadId.hpp"

#include "algorithm.hpp"
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
    vector< array<uint32_t, 2> > ordinals;
};



class ChanZuckerberg::shasta::AlignmentInfo {
public:

    // Alignment information for each of the oriented reads in the alignment.
    class Data {
    public:

        // The total number of markers in this oriented read.
        uint32_t markerCount;

        // The ordinal of the first and last marker of this oriented read
        // involved in the alignment.
        uint32_t firstOrdinal;
        uint32_t lastOrdinal;

        // Sanity check.
        void check() const
        {
            CZI_ASSERT(firstOrdinal < markerCount);
            CZI_ASSERT(lastOrdinal < markerCount);
        }
    };
    array<Data, 2> data;



    // The number of markers in the alignment.
    // This is the same for both reads!
    // It is guaranteed to never be zero for a valid alignment.
    uint32_t markerCount;



    // Constructors.
    AlignmentInfo(
        const Alignment& alignment,
        const array<uint32_t, 2>& markerCounts)
    {
        create(alignment, markerCounts);
    }
    AlignmentInfo(
        const Alignment& alignment,
        uint32_t markerCount0,
        uint32_t markerCount1)
    {
        create(alignment, array<uint32_t, 2>({markerCount0, markerCount1}));
    }
    void create(
        const Alignment& alignment,
        const array<uint32_t, 2>& markerCounts)
    {
        // Store the number of markers in the alignment.
        markerCount = uint32_t(alignment.ordinals.size());
        CZI_ASSERT(markerCount > 0);

        // Store alignment information for each of the two oriented reads.
        for(size_t i=0; i<2; i++) {
            Data& d = data[i];
            d.markerCount = markerCounts[i];
            d.firstOrdinal = alignment.ordinals.front()[i];
            d.lastOrdinal  = alignment.ordinals.back() [i];
            d.check();
        }
    }
    AlignmentInfo() : markerCount(0) {}



    // Update the alignment to reflect a swap the two oriented reads.
    void swap()
    {
        std::swap(data[0], data[1]);
    }

    // Update the alignment to reflect reverse complementing of the two reads.
    // Takes as input the total number of markers in each read.
    void reverseComplement(
        uint32_t markerCount0,
        uint32_t markerCount1)
    {
        for(size_t i=0; i<2; i++) {
            const uint32_t& markerCount = (i==0 ? markerCount0 : markerCount1);
            Data& d = data[i];
            std::swap(d.firstOrdinal, d.lastOrdinal);
            d.firstOrdinal = markerCount - 1 - d.firstOrdinal;
            d.lastOrdinal  = markerCount - 1 - d.lastOrdinal;
        }
    }



    // Compute the left and right trim, expressed in markers.
    // This is the minimum number of markers (over the two oriented reads)
    // that are excluded from the alignment on each side.
    // If the trim is too high, the alignment is suspicious.
    pair<uint32_t, uint32_t> computeTrim(
        uint32_t markerCount0,
        uint32_t markerCount1) const
    {
        CZI_ASSERT(markerCount > 0);

        const uint32_t leftTrim = min(
            data[0].firstOrdinal,
            data[1].firstOrdinal);
        const uint32_t rightTrim = min(
            markerCount0 - 1 - data[0].lastOrdinal,
            markerCount1 - 1 - data[1].lastOrdinal);
        return make_pair(leftTrim, rightTrim);
    }
};



class ChanZuckerberg::shasta::AlignmentData :
    public ChanZuckerberg::shasta::OrientedReadPair{
public:

    // The AlignmentInfo computed with the first read on strand 0.
    AlignmentInfo info;

    AlignmentData() {}
    AlignmentData(
        const array<ReadId, 2>& readIds,
        bool isSameStrand,
        const AlignmentInfo& info) :
        OrientedReadPair(readIds, isSameStrand),
        info(info)
    {}
    AlignmentData(
        const OrientedReadPair& orientedReadPair,
        const AlignmentInfo& info) :
        OrientedReadPair(orientedReadPair),
        info(info)
    {}

};

#endif
