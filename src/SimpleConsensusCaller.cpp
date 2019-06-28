#include "SimpleConsensusCaller.hpp"
#include "Coverage.hpp"
using namespace ::shasta;
using namespace ChanZuckerberg::shasta;


Consensus SimpleConsensusCaller::operator()(
    const Coverage& coverage) const
{
    const AlignedBase base = coverage.mostFrequentBase();
    const size_t repeatCount = coverage.mostFrequentRepeatCount(base);
    return Consensus(base, repeatCount);
}

