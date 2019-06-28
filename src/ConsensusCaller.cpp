#include "ConsensusCaller.hpp"
#include "Coverage.hpp"
using namespace ::shasta;
using namespace ChanZuckerberg::shasta;



// Given a vector of Coverage objects,
// find the repeat counts that have non-zero coverage on the called base
// at any position.
std::set<size_t> ConsensusCaller::findRepeatCounts(
    const vector<Coverage>& coverages) const
{

    std::set<size_t> repeatCounts;
    for(const Coverage& coverage: coverages) {
        const AlignedBase base = (*this)(coverage).base;
        if(base.isGap()) {
            continue;
        }
        const size_t repeatCountEnd = coverage.repeatCountEnd(base);
        for(size_t repeatCount=0; repeatCount<repeatCountEnd; repeatCount++) {
            if(coverage.coverage(base, repeatCount)) {
                repeatCounts.insert(repeatCount);
            }
        }
    }
    return repeatCounts;

}
