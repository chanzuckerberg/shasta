#include "Alignment.hpp"
using namespace shasta;

uint32_t Alignment::maxSkip() const
{
    uint32_t returnValue = 0;
    for(uint64_t i=1; i<ordinals.size(); i++) {
        const uint32_t skip = max(
            abs(int32_t(ordinals[i][0]) - int32_t(ordinals[i-1][0])),
            abs(int32_t(ordinals[i][1]) - int32_t(ordinals[i-1][1]))
            );
        returnValue = max(returnValue, skip);
    }
    return returnValue;
}



uint32_t Alignment::maxDrift() const
{
    uint32_t returnValue = 0;
    for(uint64_t i=1; i<ordinals.size(); i++) {
        const int32_t offset =  int32_t(ordinals[i][0]) - int32_t(ordinals[i][1]);
        const int32_t previousOffset =  int32_t(ordinals[i-1][0]) - int32_t(ordinals[i-1][1]);
        const uint32_t drift = abs(offset - previousOffset);
        returnValue = max(returnValue, drift);
    }
    return returnValue;
}



void Alignment::swap()
{
    for(auto& o: ordinals) {
        std::swap(o[0], o[1]);
    }
}



void Alignment::reverseComplement(
    uint32_t markerCount0,
    uint32_t markerCount1)
{
    for(auto& o: ordinals) {
        o[0] = markerCount0 - 1 - o[0];
        o[1] = markerCount1 - 1 - o[1];
    }
    reverse(ordinals.begin(), ordinals.end());
}



void Alignment::checkStrictlyIncreasing() const
{
    for(uint64_t i=1; i<ordinals.size(); i++) {
        const auto& previous = ordinals[i-1];
        const auto& current = ordinals[i];
        SHASTA_ASSERT(previous[0] < current[0]);
        SHASTA_ASSERT(previous[1] < current[1]);
    }
}
