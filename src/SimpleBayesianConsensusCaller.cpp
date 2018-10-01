#include "SimpleBayesianConsensusCaller.hpp"
#include "Coverage.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



// The constructor does not have any parameters.
// All data should be read from a file with fixed name
// in the run directory. We will update the documentation accordingly.
SimpleBayesianConsensusCaller::SimpleBayesianConsensusCaller()
{
    throw runtime_error("SimpleBayesianConsensusCaller constructor: not implemented.");
}



Consensus SimpleBayesianConsensusCaller::operator()(
    const Coverage& coverage) const
{

    // Access the vector containing (base, strand, repeat count)
    // for each read covering the alignment position.
    const vector<CoverageData>& coverageDataVector = coverage.getReadCoverageData();

    // Loop over supporting reads.
    for(const CoverageData& coverageData: coverageDataVector) {
        cout << coverageData.base << " ";
        cout << coverageData.strand << " ";
        cout << coverageData.repeatCount << endl;
    }
    throw runtime_error("SimpleBayesianConsensusCaller::operator(): not implemented.");

    // This returns a A with repeat count 3.
    return Consensus(AlignedBase::fromCharacter('A'), 3);
}
