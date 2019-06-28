#ifndef SHASTA_CONSENSUS_CALLER_HPP
#define SHASTA_CONSENSUS_CALLER_HPP


/*******************************************************************************

Class ConsensusCaller is an abstract base class.

A concrete ConsensusCaller contains an algorithm for
computing a "best" base and repeat count at a single
position of a multiple sequence alignment.
The Coverage object contains information about read
coverage at that position of the alignment.

Derived classes must implement operator(), which must
run the consensus algorithm and return a pair containing
the "best" base and repeat count.

*******************************************************************************/

// Shasta
#include "Base.hpp"

// Standard libraries.
#include <set>
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    using namespace ChanZuckerberg::shasta;
    class Coverage;
    class ConsensusCaller;
    class Consensus;
}



// Class used to represent the consensus base and repeat count
// at a position of an alignment.
class shasta::Consensus {
public:
    AlignedBase base;
    size_t repeatCount;

    Consensus(AlignedBase base = AlignedBase::gap(), size_t repeatCount = 0) :
        base(base), repeatCount(repeatCount) {}
};



class shasta::ConsensusCaller {
public:

    // Function that, given a Coverage object, returns
    // the "best" base and repeat count, using the
    // algorithm implemented by the derived class.
    // This is the only pure virtual function.
    // It must be implemented by all derived classes.
    virtual Consensus operator()(const Coverage&) const = 0;

    // Virtual destructor, to ensure destruction of derived classes.
    virtual ~ConsensusCaller() {}

    // Given a vector of ConsensusInfo objects,
    // find the repeat counts that have non-zero coverage on the called base
    // at any position.
    std::set<size_t> findRepeatCounts(const vector<Coverage>&) const;
};



#endif
