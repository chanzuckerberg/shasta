#ifndef CZI_SHASTA_CONSENSUS_CALLER_HPP
#define CZI_SHASTA_CONSENSUS_CALLER_HPP


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
#include "utility.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class Coverage;
        class ConsensusCaller;
    }
}



class ChanZuckerberg::shasta::ConsensusCaller {
public:

    // Function that, given a Coverage object, returns
    // the "best" base and repeat count, using the
    // algorithm implemented by the derived class.
    // This is the only pure virtual function.
    // It must be implemented by all derived classes.
    virtual pair<AlignedBase, size_t> operator()(const Coverage&) const = 0;

    // Virtual destructor, to ensure destruction of derived classes.
    virtual ~ConsensusCaller() {}
};



#endif
