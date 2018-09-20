#ifndef CZI_SHASTA_CONSENSUS_INFO_HPP
#define CZI_SHASTA_CONSENSUS_INFO_HPP



// Class ConsensusInfo is used to summarize coverage information
// at a single position of a multiple sequence alignment.
// It contains read coverage for each base
// and for each repeat count.

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

    // Increment coverage for a given base and repeat count.
    void incrementCoverage(Base, size_t repeatCount);

    // Increment coverage for '-'.
    void incrementGapCoverage();

    // Get coverage for a given base, for all repeat counts.
    // The base can be ACGT or '-'.
    size_t getCoverage(AlignedBase) const;

    // Get coverage for a given base and repeat count.
    // The base cannot be '-'.
    size_t getCoverage(AlignedBase, size_t repeatCount) const;

    // void incrementRepeatCountCoverage(size_t baseIndex, size_t repeatCount);
    size_t maxRepeatCount(size_t baseIndex) const;

    void computeBestBaseBestRepeatCount();

    // Get base coverage for the best base.
    size_t bestBaseCoverage() const;

    // Coverage for each base at this position.
    // Indexed by AlignedBase::value.
    array<size_t, 5> baseCoverage = {{0, 0, 0, 0, 0}};

    // Return the base with the most coverage.
    // This can return ACGT or '-'.
    AlignedBase bestBase() const;

    // Coverage for individual repeat counts for each base.
    // Indexed by Base::value.
    // Note that this includes entries for ACGT only.
    // , no entry for '-').
    array<vector<size_t>, 4> repeatCountCoverage;

    // The best repeat count for the best base.
    // Will be 0 if bestBaseCharacter=='-'.
    size_t bestBaseBestRepeatCount = 0;
};



#endif

