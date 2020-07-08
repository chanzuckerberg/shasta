#ifndef SHASTA_ALIGNMENT_HPP
#define SHASTA_ALIGNMENT_HPP

#include "OrientedReadPair.hpp"
#include "ReadId.hpp"

#include "algorithm.hpp"
#include "array.hpp"
#include <cmath>
#include "cstdint.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {

    class Alignment;
    class AlignmentData;
    class AlignmentInfo;
    enum class AlignmentType;
    void reverse(AlignmentType&);

    inline uint32_t computeOverlappingMarkerCount(
        uint32_t markerCount0,
        uint32_t markerCount1,
        int32_t ordinalOffset
    );
}



class shasta::Alignment {
public:

    // The ordinals in each of the two oriented reads of the
    // markers in the alignment.
    vector< array<uint32_t, 2> > ordinals;
    
    void clear() {
        ordinals.clear();
    }

    uint32_t maxSkip() const;
    uint32_t maxDrift() const;

    void swap();
    void reverseComplement(uint32_t markerCount0, uint32_t markerCount1);

    void checkStrictlyIncreasing() const;

};



// Enum used to classify an alignment.
enum class shasta::AlignmentType {
    read0IsContained,   // 0 is contained in 1. Draw as 0tee--1.
    read1IsContained,   // 1 is contained in 0. Draw as 1tee--0.
    read0IsBackward,    // No containement, 0 is backward of 1 at both ends. Draw as 0->1.
    read1IsBackward,    // No containement, 1 is backward of 0 at both ends. Draw as 1->0.
    ambiguous           // Draw as 0diamond--diamond1
};
inline void shasta::reverse(AlignmentType& alignmentType)
{
    switch(alignmentType) {
    case AlignmentType::read0IsContained:
        alignmentType = AlignmentType::read1IsContained;
        return;
    case AlignmentType::read1IsContained:
        alignmentType = AlignmentType::read0IsContained;
        return;
    case AlignmentType::read0IsBackward:
        alignmentType = AlignmentType::read1IsBackward;
        return;
    case AlignmentType::read1IsBackward:
        alignmentType = AlignmentType::read0IsBackward;
        return;
    case AlignmentType::ambiguous:
        return;
    default:
        SHASTA_ASSERT(0);
    }
}



class shasta::AlignmentInfo {
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
            // Use <= to allow for case where the alignment has no markers,
            // in which case everyting is set to zero.
            SHASTA_ASSERT(firstOrdinal <= markerCount);
            SHASTA_ASSERT(lastOrdinal <= markerCount);
        }

        // Update to reflect reverse complementing of the oriented read.
        void reverseComplement()
        {
            std::swap(firstOrdinal, lastOrdinal);
            firstOrdinal = markerCount - 1 - firstOrdinal;
            lastOrdinal  = markerCount - 1 - lastOrdinal;
        }

        double alignmentCenter() const
        {
            return double(firstOrdinal + lastOrdinal) / 2.;
        }
        uint32_t twiceAlignmentCenter() const
        {
            return firstOrdinal + lastOrdinal;
        }

        double centerPosition() const
        {
            return double(markerCount) / 2.;
        }
        uint32_t twiceCenterPosition() const
        {
            return markerCount;
        }

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

    // The minimum, maximum, and average ordinal offset,
    // computed over all aligned markers.
    int32_t minOrdinalOffset;
    int32_t maxOrdinalOffset;
    int32_t averageOrdinalOffset;



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

        // Store alignment information for each of the two oriented reads.
        for(size_t i=0; i<2; i++) {
            data[i] = Data(
                markerCounts[i],
                (markerCount == 0) ? 0 : alignment.ordinals.front()[i],
                (markerCount == 0) ? 0 : alignment.ordinals.back()[i]);
            data[i].check();
        }

        // Compute minimum, maximum, and average ordinal offset.
        minOrdinalOffset = std::numeric_limits<int32_t>::max();
        maxOrdinalOffset = std::numeric_limits<int32_t>::min();
        double sum = 0.;
        for(const auto& ordinals : alignment.ordinals) {
            const int32_t offset =
                int32_t(ordinals[0]) - int32_t(ordinals[1]);
            minOrdinalOffset = min(minOrdinalOffset, offset);
            maxOrdinalOffset = max(maxOrdinalOffset, offset);
            sum += double(offset);
        }
        averageOrdinalOffset = int32_t(std::round(sum / double(markerCount)));
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
        minOrdinalOffset = -minOrdinalOffset;
        maxOrdinalOffset = -maxOrdinalOffset;
        averageOrdinalOffset = -averageOrdinalOffset;
    }

    // Update to reflect reverse complementing of the two oriented reads.
    void reverseComplement()
    {
        for(size_t i=0; i<2; i++) {
            data[i].reverseComplement();
        }
        const int32_t delta = int32_t(data[0].markerCount) - int32_t(data[1].markerCount);
        minOrdinalOffset =  delta - minOrdinalOffset;
        maxOrdinalOffset =  delta - maxOrdinalOffset;
        averageOrdinalOffset =  delta - averageOrdinalOffset;
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
    double minAlignedFraction() const
    {
        return min(alignedFraction(0), alignedFraction(1));
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

    // Find out if this is a containing alignment,
    // that is, if the alignment covers one read
    // entirely, except possibly for up to maxTim
    // markers on each side.
    bool isContaining(uint32_t maxTrim) const {
        for(size_t i=0; i<2; i++) {
            if(leftTrim(i)<=maxTrim && rightTrim(i)<=maxTrim) {
                return true;
            }
        }
        return false;
    }



    // Return true if this alignment has read i unambiguously contained
    // in read 1-i.
    bool isContained(int i, uint32_t maxTrim) const {

        // Figure out if the two reads are fully covered by the alignment
        // (except possibly up to maxTrim markers at eah end).
        array<bool, 2> coversFullRead;
        coversFullRead[0] = (leftTrim(0)<=maxTrim && rightTrim(0)<=maxTrim);
        coversFullRead[1] = (leftTrim(1)<=maxTrim && rightTrim(1)<=maxTrim);

        // Return true only if the first read is fully covered
        // and the second read is not.
        return coversFullRead[i] and not coversFullRead[1-i];
    }



    // Jaccard similarity of the alignment.
    double jaccard() const
    {
        const uint32_t range0 = range(0);
        const uint32_t range1 = range(1);
        const double intersectionRange = double(range0 + range1) / 2.;
        const double unionRange = double(data[0].markerCount + data[1].markerCount) - intersectionRange;
        return intersectionRange / unionRange;
    }


    // Return the offset between the centers of the two oriented reads,
    // as estimated form the alignment.
    // The offset is positive if the center of the second read
    // is to the right of the center of the first read.
    double offsetAtCenter() const
    {
        const double alignmentCenter0 = data[0].alignmentCenter();
        const double alignmentCenter1 = data[1].alignmentCenter();

        const double center0 = data[0].centerPosition();
        const double center1 = data[1].centerPosition();

        return
            (center1 - alignmentCenter1) -
            (center0 - alignmentCenter0);
    }
    int32_t twiceOffsetAtCenter() const
    {
        const uint32_t twiceAlignmentCenter0 = data[0].twiceAlignmentCenter();
        const uint32_t twiceAlignmentCenter1 = data[1].twiceAlignmentCenter();

        const uint32_t twiceCenter0 = data[0].twiceCenterPosition();
        const uint32_t twiceCenter1 = data[1].twiceCenterPosition();

        return
            (int(twiceCenter1) - int(twiceAlignmentCenter1)) -
            (int(twiceCenter0) - int(twiceAlignmentCenter0));
    }


    // Classify this alignment.
    AlignmentType classify(uint32_t maxTrim) const
    {
        // Compute trim.
        const uint32_t leftTrim0  = leftTrim(0);
        const uint32_t leftTrim1  = leftTrim(1);
        const uint32_t rightTrim0 = rightTrim(0);
        const uint32_t rightTrim1 = rightTrim(1);

        // Check for containment.
        const bool isContained0 = (leftTrim0<=maxTrim) && (rightTrim0<=maxTrim);
        const bool isContained1 = (leftTrim1<=maxTrim) && (rightTrim1<=maxTrim);
        if(isContained0 && !isContained1) {
            // 0 is unambiguously contained in 1.
            return AlignmentType::read0IsContained;
        }
        if(isContained1 && !isContained0) {
            // 1 is unambiguously contained in 0.
            return AlignmentType::read1IsContained;
        }
        if(isContained0 && isContained1) {
            // Near complete overlap.
            return AlignmentType::ambiguous;
        }

        // If getting here, no containment found.
        SHASTA_ASSERT(!isContained0 && !isContained1);

        // Figure out if one of the two reads is backward at both ends.
        const bool read0IsBackward =
            leftTrim0>maxTrim  && rightTrim0<=maxTrim &&
            leftTrim1<=maxTrim && rightTrim1>=maxTrim;
        const bool read1IsBackward =
            leftTrim1>maxTrim  && rightTrim1<=maxTrim &&
            leftTrim0<=maxTrim && rightTrim0>=maxTrim;
        if(read0IsBackward && !read1IsBackward) {
            return AlignmentType::read0IsBackward;
        }
        if(read1IsBackward && !read0IsBackward) {
            return AlignmentType::read1IsBackward;
        }
        return AlignmentType::ambiguous;
    }

    void write(ostream& s) const
    {
        s << "Alignment with " << markerCount << " aligned markers:\n";
        for(int i=0; i<2; i++) {
            const auto& d = data[i];
            s << "    " << i;
            s << ": first " << d.firstOrdinal;
            s << ", last: " << d.lastOrdinal;
            s << ", total: " << d.markerCount << "\n";
        }
        s << "    Offset at center: " << offsetAtCenter() << endl;
    }
};



class shasta::AlignmentData :
    public shasta::OrientedReadPair {
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



// Compute the number of overlapping markers between two
// oriented reads with given number of markers
// and ordinal offset.
inline uint32_t shasta::computeOverlappingMarkerCount(
    uint32_t markerCount0,
    uint32_t markerCount1,
    int32_t ordinalOffset)
{
    const int32_t begin0 = 0;
    const int32_t end0 = int32_t(markerCount0);
    const int32_t begin1 = ordinalOffset;
    const int32_t end1 = begin1 + int32_t(markerCount1);

    const int32_t overlapBegin = max(begin0, begin1);
    const int32_t overlapEnd = min(end0, end1);

    if(overlapEnd < overlapBegin) {
        return 0;
    } else {
        return overlapEnd - overlapBegin;
    }

}

#endif
