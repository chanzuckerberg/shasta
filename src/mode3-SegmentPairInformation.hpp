#ifndef SHASTA_MODE3_SEGMENT_PAIR_INFORMATION_HPP
#define SHASTA_MODE3_SEGMENT_PAIR_INFORMATION_HPP

// Shasta.
#include "invalid.hpp"
#include "SHASTA_ASSERT.hpp"

// Standard library.
#include "algorithm.hpp"
#include "array.hpp"
#include "cstdint.hpp"

namespace shasta {
    namespace mode3 {
        class SegmentPairInformation;
    }
}



// Information for a pair of segments, as computed by
// mode3::AssemblyGraph::analyzeSegmentPair.
class shasta::mode3::SegmentPairInformation {
public:

    // The total number of oriented reads present in each segment.
    array<uint64_t, 2> totalCount = {0, 0};

    // The number of oriented reads present in both segments.
    // If this is zero, the rest of the information is not valid.
    uint64_t commonCount = 0;

    // The offset of segment 1 relative to segment 0, in markers.
    int64_t offset = invalid<int64_t>;

    // The number of oriented reads present in each segment
    // but missing from the other segment,
    // and which should have been present based on the above estimated offset.
    array<uint64_t, 2> unexplainedCount = {0, 0};

    // The number of oriented reads that appear in only one
    // of the two segments, but based on the estimated offset
    // are too short to appear in the other segment.
    array<uint64_t, 2> shortCount = {0, 0};

    // Check that the above counts are consistent.
    void check() const
    {
        for(uint64_t i=0; i<2; i++) {
            SHASTA_ASSERT(commonCount + unexplainedCount[i] + shortCount[i] ==
                totalCount[i]);
        }
    }

    // This computes the fraction of unexplained oriented reads,
    // without counting the short ones.
    double unexplainedFraction(uint64_t i) const
    {
        // return double(unexplainedCount[i]) / double(totalCount[i]);
        return double(unexplainedCount[i]) / double(commonCount + unexplainedCount[i]);
    }
    double maximumUnexplainedFraction() const
    {
        return max(unexplainedFraction(0), unexplainedFraction(1));
    }

    // Jaccard similarity, without counting the short reads.
    double jaccard() const
    {
        return double(commonCount) / double(commonCount + unexplainedCount[0] + unexplainedCount[1]);
    }

    // Raw Jaccard similarity (no special treatment of short reads)
    double rawJaccard() const
    {
        return double(commonCount) / double(totalCount[0] + totalCount[1] - commonCount);
    }
};

#endif
