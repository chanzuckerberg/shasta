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

        Data() : markerCount(0) {}
        Data(
            uint32_t markerCount,
            uint32_t firstOrdinal,
            uint32_t lastOrdinal) :
            markerCount(markerCount),
            firstOrdinal(firstOrdinal),
            lastOrdinal(lastOrdinal) {}


        // Return the number of markers in this oriented read
        // that are on the left of the alignment, before the alignment begins.
        uint32_t leftTrim() const
        {
            return firstOrdinal;
        }

        // Return the number of markers in this oriented read
        // that are on the right of the alignment, after the alignment ends.
        uint32_t rightTrim() const
        {
            return markerCount - 1 - lastOrdinal;
        }

        // Return the number of markers in the range covered by the alignment.
        uint32_t range() const
        {
            return lastOrdinal + 1 - firstOrdinal;
        }


        // Sanity check.
        void check() const
        {
            CZI_ASSERT(firstOrdinal < markerCount);
            CZI_ASSERT(lastOrdinal < markerCount);
        }

        // Update to reflect reverse complementing of the oriented read.
        void reverseComplement()
        {
            std::swap(firstOrdinal, lastOrdinal);
            firstOrdinal = markerCount - 1 - firstOrdinal;
            lastOrdinal  = markerCount - 1 - lastOrdinal;
        }

    private:

        // The total number of markers in this oriented read.
        uint32_t markerCount;

        // The ordinal of the first and last marker of this oriented read
        // involved in the alignment.
        uint32_t firstOrdinal;
        uint32_t lastOrdinal;

    };
    array<Data, 2> data;



    // The number of markers in the alignment.
    // This is the same for both oriented reads!
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
            data[i] = Data(
                markerCounts[i],
                alignment.ordinals.front()[i],
                alignment.ordinals.back() [i]);
            data[i].check();
        }
    }
    void create(
        const Alignment& alignment,
        uint32_t markerCount0,
        uint32_t markerCount1)
    {
        create(alignment, array<uint32_t, 2>({markerCount0, markerCount1}));
    }
    AlignmentInfo() : markerCount(0) {}



    // Update to reflect a swap the two oriented reads.
    void swap()
    {
        std::swap(data[0], data[1]);
    }

    // Update to reflect reverse complementing of the two oriented reads.
    void reverseComplement()
    {
        for(size_t i=0; i<2; i++) {
            data[i].reverseComplement();
        }
    }

    // Some accessors.
    uint32_t leftTrim(size_t i) const {
        return data[i].leftTrim();
    }
    uint32_t rightTrim(size_t i) const {
        return data[i].rightTrim();
    }
    uint32_t range(size_t i) const {
        return data[i].range();
    }

    // Return the ratio of aligned markers over the alignment range.
    double alignedFraction(size_t i) const
    {
        return double(markerCount) / double(range(i));
    }

    // Compute the left and right trim, expressed in markers.
    // This is the minimum number of markers (over the two oriented reads)
    // that are excluded from the alignment on each side.
    // If the trim is too high, the alignment is suspicious.
    pair<uint32_t, uint32_t> computeTrim() const
    {
        const uint32_t leftTrim  = min(data[0].leftTrim() , data[1].leftTrim() );
        const uint32_t rightTrim = min(data[0].rightTrim(), data[1].rightTrim());
        return make_pair(leftTrim, rightTrim);
    }
};



class ChanZuckerberg::shasta::AlignmentData :
    public ChanZuckerberg::shasta::OrientedReadPair {
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
