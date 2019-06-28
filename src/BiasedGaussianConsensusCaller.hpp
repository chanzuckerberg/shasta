#ifndef SHASTA_BIASED_GAUSSIAN_CONSENSUS_CALLER_HPP
#define SHASTA_BIASED_GAUSSIAN_CONSENSUS_CALLER_HPP

/*******************************************************************************

A BiasedGaussianConsensusCaller is a ConsensusCaller
that uses the following algorithm to compute
the "best" base and repeat count:

- The best base is the base with the most coverage,
with ties broken in favor of the "earlier" base in
ACGT order.
- The best repeat count computed by averaging the
repeat counts for the best base, and then multiplying
the result by a bias factor.

*******************************************************************************/


#include "ConsensusCaller.hpp"

namespace shasta {
    class BiasedGaussianConsensusCaller;
}


class shasta::BiasedGaussianConsensusCaller :
    public shasta::ConsensusCaller {
public:

    virtual Consensus operator()(const Coverage&) const;

};

#endif
