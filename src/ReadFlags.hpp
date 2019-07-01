#ifndef SHASTA_READ_FLAGS_HPP
#define SHASTA_READ_FLAGS_HPP

#include <cstdlib>
#include "cstdint.hpp"

namespace shasta {
    class ReadFlags;
}

class shasta::ReadFlags {
public:

    // This is set for reads that are approximate palindromic,
    // that is, are well aligned with their own reverse complement.
    uint8_t isPalindromic : 1;

    // Set if the read is marked as chimeric.
    uint8_t isChimeric : 1;

    // Set if the read belongs to a small component of the read graph
    // that is not used for assembly.
    // If isChimeric is set, this is also set.
    uint8_t isInSmallComponent : 1;

    // Strand used when assembling this read.
    // If 0, the read is assembled unchanged.
    // If 1, the read is assembled reverse complemented.
    // Not valid if isPalindromic, isChimeric or isInSmallComponent is set.
    uint8_t strand : 1;

    // Unused bits.
    uint8_t bit4 : 1;
    uint8_t bit5 : 1;
    uint8_t bit6 : 1;
    uint8_t bit7 : 1;
    ReadFlags()
    {
        static_assert(sizeof(ReadFlags) == 1, "Unexpected size of ReadFlags.");
        *reinterpret_cast<uint8_t*>(this) = 0;
    }
};

#endif
