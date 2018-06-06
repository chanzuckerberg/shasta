#ifndef CZI_NANOPORE2_MARKER_ID_HPP
#define CZI_NANOPORE2_MARKER_ID_HPP

// CZI.
#include "CZI_ASSERT.hpp"

// Standard libraries.
#include "cstdint.hpp"
#include "iostream.hpp"
#include <limits>
#include "string.hpp"



namespace ChanZuckerberg {
    namespace Nanopore2 {

        // Type used to globally identify a marker on a read.
        // This is the global index of the marker in Assembler::markers.
        // For a human assembly with good coverage the total number
        // of markers can be 10 billions or more, so this needs to
        // be uint64_t. There could, however, be situations
        // where uint32_t is sufficient.
        using MarkerId = uint64_t;

        // Class used to globally identify identify a marker on an oriented read,
        // that is a on read, possibly reverse complemented.
        class OrientedMarkerId;

    }
}



// Class used to globally identify a marker on oriented read,
// that is a on a read, possibly reverse complemented.
// The strand stored in the least significant
// bit is 0 if the oriented read is identical
// to the original read and 1 if it is reverse complemented.
class ChanZuckerberg::Nanopore2::OrientedMarkerId {
public:
    OrientedMarkerId() : value(std::numeric_limits<MarkerId>::max()) {}
    OrientedMarkerId(MarkerId readId, uint64_t strand) : value((readId<<1ULL) | strand)
    {
        CZI_ASSERT(strand < 2);
    }
    explicit OrientedMarkerId(MarkerId value) : value(value) {}
    MarkerId getMarkerId() const
    {
        return value >> 1ULL;
    }
    uint64_t getStrand() const
    {
        return value & 1ULL;
    }

    void switchStrand()
    {
        value ^= 1ULL;
    }

    // Return the integer value, which can be used as an index into Assembler::markers.
    MarkerId getValue() const
    {
        return value;
    }

    bool operator==(const OrientedMarkerId& that) const
    {
        return value == that.value;
    }
    bool operator<(const OrientedMarkerId& that) const
    {
        return value < that.value;
    }

private:
    MarkerId value;
};



#endif
