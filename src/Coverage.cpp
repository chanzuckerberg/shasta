#include "Coverage.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



// Add information about a supporting read.
// If the AlignedBase is '-',repeatCount must be zero.
void Coverage::addRead(AlignedBase base, Strand strand, size_t repeatCount)
{
    // Store a CoverageData for this read.
    readCoverageData.push_back(CoverageData(base, strand, repeatCount));

    // Increment coverage.
    if(base.isGap()) {
        incrementGapCoverage();
    } else {
        incrementCoverage(Base(base), repeatCount);
    }

}



// Increment coverage for a given base and repeat count.
void Coverage::incrementCoverage(Base base, size_t repeatCount)
{
    // Extract the base value and check it.
    const uint8_t baseValue = base.value;
    CZI_ASSERT(baseValue < 4);

    // Increment total coverage for this base.
    ++baseCoverage[baseValue];

    // Increment coverage for this base and repeat count,
    // extending the vector if necessary.
    auto& v = repeatCountCoverage[baseValue];
    if(v.size() <= repeatCount) {
        v.resize(repeatCount+1, 0);
    }
    ++v[repeatCount];
}



// Increment base coverage for '-'.
void Coverage::incrementGapCoverage()
{
    ++baseCoverage[4];
}



// Return the base with the most coverage.
// This can return ACGT or '-'.
AlignedBase Coverage::bestBase() const
{
    const size_t bestBaseValue =
        std::max_element(baseCoverage.begin(), baseCoverage.end()) - baseCoverage.begin();
    return AlignedBase::fromInteger(bestBaseValue);
}



// Get the repeat count with the most coverage for a given base.
// The base canot be '-'.
size_t Coverage::bestRepeatCount(Base base) const
{
    // Extract the base value and check it.
    const uint8_t baseValue = base.value;
    CZI_ASSERT(baseValue < 4);

    // Access the repeat count coverage vector for this base.
    const auto& v = repeatCountCoverage[baseValue];

    // Return the index with maximum coverage.
    return std::max_element(v.begin(), v.end()) - v.begin();
}



// Get the repeat count with the most coverage for the base
// with the most coverage.
// The should only be called if the base with the best coverage
// is not '-'.
size_t Coverage::bestBaseBestRepeatCount() const
{
    return bestRepeatCount(Base(bestBase()));
}



// Represent a coverage value with a single character.
char Coverage::coverageCharacter(size_t coverage)
{
    if(coverage == 0) {
        return '.';
    } else if(coverage < 10) {
        const string coverageString = to_string(coverage);
        CZI_ASSERT(coverageString.size() == 1);
        return coverageString[0];
    } else {
        return '*';
    }
}



// Get coverage for a given base, for all repeat counts.
// The base can be ACGT or '-'.
size_t Coverage::coverage(AlignedBase base) const
{
    // Extract the base value and check it.
    const uint8_t baseValue = base.value;
    CZI_ASSERT(baseValue < 5);

    // Return total coverage for this base, for all repeat counts.
    return baseCoverage[baseValue];

}
char Coverage::coverageCharacter(AlignedBase base) const
{
    return coverageCharacter(coverage(base));
}



// Get coverage for a given base and repeat count.
// The base cannot be '-'.
size_t Coverage::coverage(Base base, size_t repeatCount) const
{
    // Extract the base value and check it.
    const uint8_t baseValue = base.value;
    CZI_ASSERT(baseValue < 4);

    // Access the coverage vector for this base.
    const auto& v = repeatCountCoverage[baseValue];

    // Return coverage for the given repeat count.
    if(repeatCount < v.size()) {
        return v[repeatCount];
    } else {
        return 0;
    }
}
char Coverage::coverageCharacter(Base base, size_t repeatCount) const
{
    return coverageCharacter(coverage(base, repeatCount));
}



// Get base coverage for the best base.
size_t Coverage::bestBaseCoverage() const
{
    return baseCoverage[bestBase().value];
}
char Coverage::bestBaseCoverageCharacter() const
{
    return coverageCharacter(bestBaseCoverage());
}



// Get the maximum repeat count for a given base.
// The base can be ACGT (not '-').
size_t Coverage::maxRepeatCount(Base base) const
{
    // Extract the base value and check it.
    const uint8_t baseValue = base.value;
    CZI_ASSERT(baseValue < 4);

    // The maximum repeat count is one less than the
    // size of the repeat count coverage vector for this base.
    return repeatCountCoverage[baseValue].size() - 1;
}



// Given a vector of ConsensusInfo objects,
// find the repeat counts that have non-zero coverage on the best base
// at any position.
std::set<size_t> Coverage::findRepeatCounts(const vector<Coverage>& consensusInfos)
{

    std::set<size_t> repeatCounts;
    for(const Coverage& consensusInfo: consensusInfos) {
        const AlignedBase bestBase = consensusInfo.bestBase();
        if(bestBase.isGap()) {
            continue;
        }
        const size_t maxRepeatCount = consensusInfo.maxRepeatCount(Base(bestBase));
        for(size_t repeatCount=0; repeatCount<=maxRepeatCount; repeatCount++) {
            const size_t coverage = consensusInfo.coverage(Base(bestBase), repeatCount);
            if(coverage) {
                repeatCounts.insert(repeatCount);
            }
        }
    }
    return repeatCounts;
}







