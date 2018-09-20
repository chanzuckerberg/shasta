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



// Get coverage for a given base, for all repeat counts.
// The base can be ACGT or '-'.
size_t ConsensusInfo::getCoverage(AlignedBase base) const
{
    // Extract the base value and check it.
    const uint8_t baseValue = base.value;
    CZI_ASSERT(baseValue < 5);

    // Return total coverage for this base, for all repeat counts.
    return baseCoverage[baseValue];

}


// Get coverage for a given base and repeat count.
// The base cannot be '-'.
size_t ConsensusInfo::getCoverage(AlignedBase base, size_t repeatCount) const
{
    CZI_ASSERT(!base.isGap());

    const auto& v = repeatCountCoverage[base.value];
    if(repeatCount < v.size()) {
        return v[repeatCount];
    } else {
        return 0;
    }
}



size_t ConsensusInfo::maxRepeatCount(size_t baseIndex) const
{
    return repeatCountCoverage[baseIndex].size() - 1;
}



void ConsensusInfo::computeBestBaseBestRepeatCount()
{
    const auto& v = repeatCountCoverage[bestBase.value];
    bestBaseBestRepeatCount =
        std::max_element(v.begin(), v.end()) - v.begin();
}



size_t ::ConsensusInfo::bestBaseCoverage() const
{
    return baseCoverage[bestBase.value];
}

