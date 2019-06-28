#ifndef SHASTA_COMPUTE_RUN_LENGTH_REPRESENTATION_HPP
#define SHASTA_COMPUTE_RUN_LENGTH_REPRESENTATION_HPP

#include "Base.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {

    // Given the raw representation of a sequence, compute its
    // run-length representation.
    // This returns false if the sequence contains a homopolymer run
    // of more than 255 bases, which cannot be represented
    // with a one-byte repeat count.
    bool computeRunLengthRepresentation(
        const vector<Base>& sequence,
        vector<Base>& runLengthSequence,
        vector<uint8_t>& repeatCount);

    }
}

#endif
