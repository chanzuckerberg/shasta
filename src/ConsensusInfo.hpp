#ifndef CZI_SHASTA_CONSENSUS_INFO_HPP
#define CZI_SHASTA_CONSENSUS_INFO_HPP



// Class ConsensusInfo is used to summarize information
// at a single position of a multiple sequence alignment.
// It contains read coverage for each base and
// for each repeat count.
// It also stores the base with the best read coverage
// and, for that base, the repeat count with the most coverage.

// Shasta.
#include "Base.hpp"
#include "CZI_ASSERT.hpp"

// Standard library.
#include "algorithm.hpp"
#include "array.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class ConsensusInfo;
    }
}



class ChanZuckerberg::shasta::ConsensusInfo {
public:

    // Coverage for each of ACGT-.
    array<size_t, 5> baseCoverage = {{0, 0, 0, 0, 0}};

    // Character representing the best base of gap. Can be one of "ACGT-".
    char bestBaseCharacter = 'N';

    // Coverage for individual repeat counts for each base (ACGT only, no entry for '-').
    array<vector<size_t>, 4> repeatCountCoverage;
    size_t getRepeatCountCoverage(size_t baseIndex, size_t repeatCount) const;
    void incrementRepeatCountCoverage(size_t baseIndex, size_t repeatCount);
    size_t maxRepeatCount(size_t baseIndex) const;

    // The best repeat count for the best base.
    // Will be 0 if bestBaseCharacter=='-'.
    size_t bestBaseBestRepeatCount = 0;
    void computeBestBaseBestRepeatCount();

    // Get the best base.
    // This asserts if bestBaseCharacter is not a valid base.
    Base bestBase() const;

    // Get base coverage for the best base.
    size_t bestBaseCoverage() const;

};



#endif

