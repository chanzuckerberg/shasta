#include "ConsensusInfo.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



// Increment coverage for a given base and repeat count.
void ConsensusInfo::incrementCoverage(Base base, size_t repeatCount)
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
void ConsensusInfo::incrementGapCoverage()
{
    ++baseCoverage[4];
}



// Represent a coverage value with a single character.
char ConsensusInfo::coverageCharacter(size_t coverage)
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
size_t ConsensusInfo::coverage(AlignedBase base) const
{
    // Extract the base value and check it.
    const uint8_t baseValue = base.value;
    CZI_ASSERT(baseValue < 5);

    // Return total coverage for this base, for all repeat counts.
    return baseCoverage[baseValue];

}
char ConsensusInfo::coverageCharacter(AlignedBase base) const
{
    return coverageCharacter(coverage(base));
}



// Get coverage for a given base and repeat count.
// The base cannot be '-'.
size_t ConsensusInfo::coverage(Base base, size_t repeatCount) const
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
char ConsensusInfo::coverageCharacter(Base base, size_t repeatCount) const
{
    return coverageCharacter(coverage(base, repeatCount));
}



// Return the base with the most coverage.
// This can return ACGT or '-'.
AlignedBase ConsensusInfo::bestBase() const
{
    const size_t bestBaseValue =
        std::max_element(baseCoverage.begin(), baseCoverage.end()) - baseCoverage.begin();
    return AlignedBase::fromInteger(bestBaseValue);
}



// Get base coverage for the best base.
size_t ConsensusInfo::bestBaseCoverage() const
{
    return baseCoverage[bestBase().value];
}



size_t ConsensusInfo::maxRepeatCount(size_t baseIndex) const
{
    return repeatCountCoverage[baseIndex].size() - 1;
}



// Get the repeat count with the most coverage for a given base.
// The base canot be '-'.
size_t ConsensusInfo::getBestRepeatCount(Base base) const
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
size_t ConsensusInfo::bestBaseBestRepeatCount() const
{
    return getBestRepeatCount(Base(bestBase()));
}




