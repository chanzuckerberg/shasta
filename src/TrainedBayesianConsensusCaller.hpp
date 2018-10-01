#ifndef CZI_SHASTA_TRAINED_BAYESIAN_CONSENSUS_CALLER_HPP
#define CZI_SHASTA_TRAINED_BAYESIAN_CONSENSUS_CALLER_HPP

/*******************************************************************************

Class TrainedBayesianConsensusCaller uses a Bayesian approach
developed by Jordan Eizenga at UCSC.

*******************************************************************************/

#include "ConsensusCaller.hpp"

#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "tuple.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class TrainedBayesianConsensusCaller;
    }
}


class ChanZuckerberg::shasta::TrainedBayesianConsensusCaller :
    public ChanZuckerberg::shasta::ConsensusCaller {
public:

    // The constructor does not have any parameters.
    // All data should be read from a file with fixed name
    // in the run directory. We will update the documentation accordingly.
    TrainedBayesianConsensusCaller();

    // Function that does the computation for a given alignment position.
    // The Coverage object contains all the necessary information.
    virtual Consensus operator()(const Coverage&) const;

private:
    
    // The maximum length of a homopolymer repeat we allow
    size_t maxRepeatCount = 0;
    
    // Convenience vector that contains all possible homopolymer calls
    vector<Consensus> repeatBases;
    
    // Memo of the log(prob) of a called homopolymer given a true homopolymer
    std::map<tuple<AlignedBase, size_t, AlignedBase, size_t>, double> logConditionalProbabilities;

};

#endif
