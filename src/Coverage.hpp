#ifndef CZI_SHASTA_COVERAGE_HPP
#define CZI_SHASTA_COVERAGE_HPP



// Class ConsensusInfo is used to summarize coverage information
// at a single position of a multiple sequence alignment.
// It stores read coverage for each base (ACGT and '-')
// and for each repeat count.

// Shasta.
#include "Base.hpp"
#include "CZI_ASSERT.hpp"
#include "ReadId.hpp"

// Standard library.
#include "algorithm.hpp"
#include "array.hpp"
#include <set>
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class Coverage;
    }
}



class ChanZuckerberg::shasta::Coverage {
public:

    // Add information about a supporting read.
    // If the AlignedBase is '-',repeatCount must be zero.
    // This is the only public non-const function.
    void addRead(AlignedBase, Strand, size_t repeatCount);



    // Functions to return the best base or best repeat count.
    // Here, best means "with the most coverage".

    // Return the base with the most coverage.
    // This can return ACGT or '-'.
    AlignedBase bestBase() const;

    // Get the repeat count with the most coverage for a given base.
    // The base cannot be '-'.
    size_t bestRepeatCount(Base) const;

    // Get the repeat count with the most coverage for the base
    // with the most coverage.
    // This should only be called if the base with the best coverage
    // is not '-'.
    size_t bestBaseBestRepeatCount() const;



    // Represent a coverage value with a single character.
    static char coverageCharacter(size_t);

    // Get coverage for a given base, for all repeat counts.
    // The base can be ACGT or '-'.
    size_t coverage(AlignedBase) const;
    char coverageCharacter(AlignedBase) const;

    // Get coverage for a given base and repeat count.
    // The base cannot be '-'.
    size_t coverage(Base, size_t repeatCount) const;
    char coverageCharacter(Base, size_t repeatCount) const;

    // Get base coverage for the best base.
    size_t bestBaseCoverage() const;
    char bestBaseCoverageCharacter() const;

    // Get the maximum repeat count for a given base.
    // The base can be ACGT (not '-').
    size_t maxRepeatCount(Base) const;


    // Given a vector of ConsensusInfo objects,
    // find the repeat counts that have non-zero coverage on the best base
    // at any position.
    static std::set<size_t> findRepeatCounts(const vector<Coverage>&);

private:

    // Coverage for each base (ACGT or '-')
    // at the position described by this ConsensusInfo.
    // Indexed by AlignedBase::value.
    array<size_t, 5> baseCoverage = {{0, 0, 0, 0, 0}};

    // Coverage for individual repeat counts for each base.
    // Indexed by Base::value.
    // Note that this includes entries for ACGT only (no entry for '-').
    array<vector<size_t>, 4> repeatCountCoverage;

    // Increment coverage for a given base and repeat count.
    void incrementCoverage(Base, size_t repeatCount);

    // Increment coverage for '-'.
    void incrementGapCoverage();
};



#endif

