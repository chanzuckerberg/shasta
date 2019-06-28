#include "BiasedGaussianConsensusCaller.hpp"
#include "Coverage.hpp"
using namespace ::shasta;
using namespace ChanZuckerberg::shasta;

#include <cmath>


Consensus BiasedGaussianConsensusCaller::operator()(
    const Coverage& coverage) const
{
    const AlignedBase base = coverage.mostFrequentBase();

    // Average the repeat counts for this base.
    double sum = 0;
    double count = 0;
    const vector<CoverageData>& coverageData = coverage.getReadCoverageData();
    for(const CoverageData& cd: coverageData) {
        if(cd.base == base) {
            sum += double(cd.repeatCount);
            count += 1.;
        }
    }
    const double average = sum / count;

    const double biasFactor = 1.05;
    const double biasedAverage = biasFactor * average;

    return Consensus(base, std::lround(biasedAverage));
}

