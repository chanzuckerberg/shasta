#include "ConsensusInfo.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



size_t ConsensusInfo::getRepeatCountCoverage(size_t baseIndex, size_t repeatCount) const
{
    CZI_ASSERT(baseIndex < 4);
    const auto& v = repeatCountCoverage[baseIndex];
    if(repeatCount < v.size()) {
        return v[repeatCount];
    } else {
        return 0;
    }
}



void ConsensusInfo::incrementRepeatCountCoverage(size_t baseIndex, size_t repeatCount)
{
    CZI_ASSERT(baseIndex < 4);
    auto& v = repeatCountCoverage[baseIndex];
    if(repeatCount >= v.size()) {
        v.resize(repeatCount+1, 0);
    }
    ++v[repeatCount];
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

