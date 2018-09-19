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

    // Coverage for ACGT-.
    array<size_t, 5> baseCoverage = {{0, 0, 0, 0, 0}};

    // Character representing the best base of gap. Can be one of "ACGT-".
    char bestBaseCharacter = 'N';

    // Coverage for individual repeat counts for each base.
    array<vector<size_t>, 4> repeatCountCoverage;
    size_t getRepeatCountCoverage(size_t baseIndex, size_t repeatCount) const
    {
        CZI_ASSERT(baseIndex < 4);
        const auto& v = repeatCountCoverage[baseIndex];
        if(repeatCount < v.size()) {
            return v[repeatCount];
        } else {
            return 0;
        }
    }
    void incrementRepeatCountCoverage(size_t baseIndex, size_t repeatCount)
    {
        CZI_ASSERT(baseIndex < 4);
        auto& v = repeatCountCoverage[baseIndex];
        if(repeatCount >= v.size()) {
            v.resize(repeatCount+1, 0);
        }
        ++v[repeatCount];
    }
    size_t maxRepeatCount(size_t baseIndex) const {
        return repeatCountCoverage[baseIndex].size() - 1;
    }

    // The best repeat count for the best base.
    // Will be 0 if bestBaseCharacter=='-'.
    size_t bestBaseBestRepeatCount = 0;
    void computeBestBaseBestRepeatCount()
    {
        const auto& v = repeatCountCoverage[bestBase().value];
        bestBaseBestRepeatCount =
            std::max_element(v.begin(), v.end()) - v.begin();
    }

    // Get the best base.
    // This asserts if bestBaseCharacter is not a valid base.
    Base bestBase() const
    {
        return Base::fromCharacter(bestBaseCharacter);
    }

    // Get base coverage for the best base.
    size_t bestBaseCoverage() const
    {
        if(bestBaseCharacter == '-') {
            return baseCoverage[4];
        } else {
            return baseCoverage[bestBase().value];
        }
    }

};



#endif

