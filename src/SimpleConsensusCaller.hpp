#ifndef CZI_SHASTA_SIMPLE_CONSENSUS_CALLER_HPP
#define CZI_SHASTA_SIMPLE_CONSENSUS_CALLER_HPP

/*******************************************************************************

A SimpleConsensusCaller is a ConsensusCaller
that uses a very simple algorithm to compute
the "best" base and repeat count:

- The best base is the base with the most coverage,
with ties broken in favor of the "earlier" base in
ACGT order.
- The best repeat count is the repeat count with
the most coverage for the best base as defined above.
Ties are broken in favor or smaller repeat counts.

*******************************************************************************/


#include "ConsensusCaller.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class SimpleConsensusCaller;
    }
}


class ChanZuckerberg::shasta::SimpleConsensusCaller :
    public ChanZuckerberg::shasta::ConsensusCaller {
public:

    virtual pair<AlignedBase, size_t> operator()(const Coverage&) const;

};

#endif
