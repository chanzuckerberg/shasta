#ifndef CZI_SHASTA_SIMPLE_BAYESIAN_CONSENSUS_CALLER_HPP
#define CZI_SHASTA_SIMPLE_BAYESIAN_CONSENSUS_CALLER_HPP

/*******************************************************************************

A SimpleBayesianConsensusCaller uses a simple Bayesian approach
to compute the "best" base and repeat count at a position of an alignment.

Based on initial work by Ryan Lorigro at UCSC, the method works as follows.
Here, n is the true repeat count, m is the observed repeat count,
and mi the observed repeat counts in a set of reads.

- Once for a given sequencing technology, estimate conditional probabilities
P(m | n, base read) by mapping reads to portions of a reference
known not to contain variants.

- Use Bayes theorem to estimate

P(n | mi, base) proportional to P(n) times product over i P(m | n, base read)

Where P(n) is the prior probability of a homopolymer run of length n
and can be initially neglected.

Note that in the above, the base read at a given alignment position
must take into account which strand each read is on.

*******************************************************************************/

#include "ConsensusCaller.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class SimpleBayesianConsensusCaller;
    }
}


class ChanZuckerberg::shasta::SimpleBayesianConsensusCaller :
    public ChanZuckerberg::shasta::ConsensusCaller {
public:

    // The constructor does not have any parameters.
    // All data should be read from a file with fixed name
    // in the run directory. We will update the documentation accordingly.
    SimpleBayesianConsensusCaller();

    // Function that does the computation for a given alignment position.
    // The Coverage object contains all the necessary information.
    virtual Consensus operator()(const Coverage&) const;

private:

    // Put here any functions and data needed by the SimpleBayesianConsensusCaller.

};

#endif
