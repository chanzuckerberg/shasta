#include "SimpleBayesianConsensusCaller.hpp"
#include "Coverage.hpp"
#include "Base.hpp"


using ChanZuckerberg::shasta::SimpleBayesianConsensusCaller;
using ChanZuckerberg::shasta::Coverage;
using ChanZuckerberg::shasta::AlignedBase;
using std::string;
using std::vector;
using std::ifstream;
using std::cout;

int main(){
    SimpleBayesianConsensusCaller classifier;   // Add any required constructor parameters.
    Coverage c;

    c.addRead(AlignedBase::fromInteger((uint8_t)1), 1, 1);    // Arguments are base, strand, repeat count.
    c.addRead(AlignedBase::fromInteger((uint8_t)1), 0, 2);
    c.addRead(AlignedBase::fromInteger((uint8_t)2), 1, 3);

//    int consensus_int;
//    vector<double> loglikelihoods;
//    tie(consensus_int, loglikelihoods) = classifier.predict_runlength(c, AlignedBase::fromInteger((uint8_t)0));

    const Consensus consensus = classifier(c);

    cout << consensus.base << " " << consensus.repeatCount << '\n';

    return 0;
}