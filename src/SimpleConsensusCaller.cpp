#include "SimpleConsensusCaller.hpp"
#include "Coverage.hpp"
using namespace ChanZuckerberg;
using namespace shasta;


pair<AlignedBase, size_t> SimpleConsensusCaller::operator()(
    const Coverage& coverage) const
{
    const AlignedBase base = coverage.mostFrequentBase();
    const size_t repeatCount = coverage.mostFrequentRepeatCount(base);
    return make_pair(base, repeatCount);
}

