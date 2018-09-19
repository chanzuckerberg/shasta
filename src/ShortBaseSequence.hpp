#ifndef CZI_SHASTA_SHORT_BASE_SEQUENCE_HPP
#define CZI_SHASTA_SHORT_BASE_SEQUENCE_HPP

// shasta.
#include "Base.hpp"

// Standard library.
#include "array.hpp"
#include "iostream.hpp"
#include <limits>
#include "stdexcept.hpp"
#include "string.hpp"



namespace ChanZuckerberg {
    namespace shasta {

        // A short sequence of bases.
        // Uses only two integers, so its capacity is limited
        // by the length of the integers used.
        // This class does not keep track of the number of bases
        // actually stored. All unused positions are left set at "A".
        template<class Int> class ShortBaseSequence;
        using ShortBaseSequence8 = ShortBaseSequence<uint8_t>;
        using ShortBaseSequence16 = ShortBaseSequence<uint16_t>;
        using ShortBaseSequence32 = ShortBaseSequence<uint32_t>;
        using ShortBaseSequence64 = ShortBaseSequence<uint64_t>;
        template<class Int> inline ostream& operator<<(ostream&, const ShortBaseSequence<Int>&);

        void testShortBaseSequence();
    }
}



// A short sequence of bases.
// Uses only two integers, so its capacity is limited
// by the length of the integers used.
// Position 0: the LSB bit of the bases (with base 0 corresponding to the MSB bit).
// Position 1: the MSB bit of the bases (with base 0 corresponding to the MSB bit).
// This class does not keep track of the number of bases
// actually stored. All unused positions are left set at "A".
template<class Int> class ChanZuckerberg::shasta::ShortBaseSequence {
public:

    // Sanity check on the Int type.
    static_assert(std::numeric_limits<Int>::is_integer,
        "Int type for ShortBaseSequence must be an integer type.");
    static_assert(!std::numeric_limits<Int>::is_signed,
        "Int type for ShortBaseSequence must be an unsigned integer type.");

    // The number of bases that can be represented equals the number of bits
    // in the Int type.
    static const size_t capacity = std::numeric_limits<Int>::digits;
    static const size_t capacityMinus1 = capacity - 1ULL;

    // The constructor fills the data with 0, which corresponds to all A's.
    ShortBaseSequence()
    {
        std::fill(data.begin(), data.end(), Int(0));
    }

    // Return the base at a given position.
    Base operator[](uint64_t i) const
    {
        const uint64_t bitIndex = capacityMinus1 - (i & capacityMinus1);
        const uint64_t bit0 = (data[0] >> bitIndex) & 1ULL;
        const uint64_t bit1 = (data[1] >> bitIndex) & 1ULL;
        const uint8_t value = uint8_t((bit1 << 1ULL) + bit0);
        return Base::fromInteger(value);
    }

    // Set the base at a given position.
    void set(uint64_t i, Base base) {
        const uint64_t bitIndex = capacityMinus1 - (i & capacityMinus1);
        const Int mask = Int(1ULL << bitIndex);
        const Int maskComplement = Int(~mask);

        const uint64_t bit0 = (base.value) & 1ULL;
        if(bit0 == 0) {
            data[0] &= maskComplement;
        } else {
            data[0] |= mask;
        }

        const uint64_t bit1 = (base.value >> 1ULL) & 1ULL;
        if(bit1 == 0) {
            data[1] &= maskComplement;
        } else {
            data[1] |= mask;
        }

    }

    // Return an integer consisting of the concatenation
    // of the base bits corresponding to the first n bases.
    uint64_t id(uint64_t n) const
    {
        const uint64_t shift = capacity - n;
        const uint64_t lsb = data[0] >> shift;
        const uint64_t msb = data[1] >> shift;
        return (msb << n) | lsb;
    }

    // Opposite of the above: construct the sequence given the id.
    ShortBaseSequence(uint64_t id, uint64_t n)
    {
        const uint64_t mask = (1ULL << n) - 1ULL;
        const uint64_t shift = capacity - n;
        data[0] = Int((id & mask) << shift);
        data[1] = Int(((id >> n) & mask) << shift);
    }

    // Return the reverse complement of the first n bases.
    ShortBaseSequence<Int> reverseComplement(uint64_t n) const
    {
        ShortBaseSequence<Int> reverseComplementedSequence;
        for(size_t i=0; i<n; i++) {
            const Base b = (*this)[i].complement();
            reverseComplementedSequence.set(n-i-1, b);
        }
        return reverseComplementedSequence;
    }

    bool operator==(const ShortBaseSequence<Int>& that) const
    {
        return data == that.data;
    }

    // Write the first n bases.
    ostream& write(ostream& s, uint64_t n) const
    {
        for(uint64_t i=0; i<n; i++) {
            s << (*this)[i];
        }
        return s;
    }

    void shiftLeft() {
        data[0] = Int(data[0] << 1);
        data[1] = Int(data[1] << 1);
    }

    // The data are left public to facilitate low level custom code.
    array<Int, 2> data;
};



template<class Int> inline std::ostream& ChanZuckerberg::shasta::operator<<(
    std::ostream& s,
    const ChanZuckerberg::shasta::ShortBaseSequence<Int>& sequence)
{
    for(size_t i=0; i<sequence.capacity; i++) {
        s << sequence[i];
    }
    return s;
}



#endif
