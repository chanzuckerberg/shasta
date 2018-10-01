#ifndef CZI_SHASTA_TRAINED_BAYESIAN_CONSENSUS_CALLER_HPP
#define CZI_SHASTA_TRAINED_BAYESIAN_CONSENSUS_CALLER_HPP

/*******************************************************************************

Class TrainedBayesianConsensusCaller uses a Bayesian approach
developed by Jordan Eizenga at UCSC.

*******************************************************************************/

#include "ConsensusCaller.hpp"

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

    // Put here any functions and data needed by the TrainedBayesianConsensusCaller.

};

#endif
