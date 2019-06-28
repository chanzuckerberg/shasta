#ifndef SHASTA_SIMPLE_CONSENSUS_CALLER_HPP
#define SHASTA_SIMPLE_CONSENSUS_CALLER_HPP

/*******************************************************************************

A SimpleConsensusCaller is a ConsensusCaller
that uses a very simple algorithm to compute
the "best" base and repeat count:

- The best base is the base with the most coverage,
with ties broken in favor of the "earlier" base in
ACGT order.
- The best repeat count is the repeat count with
the most coverage for the best base as defined above.
Ties are broken in favor or larger repeat counts.
Breaking ties in this way, rather than the opposite,
improved sequence identity a bit, but the caller
is still biased in favor of deletions.

*******************************************************************************/


#include "ConsensusCaller.hpp"

        namespace shasta {
            class SimpleConsensusCaller;
        }


class shasta::SimpleConsensusCaller :
    public shasta::ConsensusCaller {
public:

    virtual Consensus operator()(const Coverage&) const;

};

#endif
