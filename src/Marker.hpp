#ifndef CZI_NANOPORE2_MARKER_HPP
#define CZI_NANOPORE2_MARKER_HPP


/*******************************************************************************

Among all 4^k k-mers of length k, we choose a subset that we call "markers".
The markers are selected at the beginning of an assembly
and never changed, and selected in such a way that,
if (and only if) a k-mer is a marker, its reverse complement
is also a marker.

The k-mer table is a vector of 4^k KmerInfo object,
indexed by k-mer id as computed using Kmer::id(k).
Because of the way markers are selected, the following is
true for all permitted values of i, 0 <= i < 4^k:
kmerTable[i].isMarker == kmerTable[kmerTable[i].reverseComplementKmerId].isMarker

*******************************************************************************/

#include "Kmer.hpp"
#include "Uint.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {

        // The classes with a 0 suffix will be phased out.
        class CompressedMarker0;
        class Marker0;
        class OrderMarkers0ByKmerId;
    }
}



class ChanZuckerberg::Nanopore2::CompressedMarker0 {
public:
    KmerId kmerId;

    // Position difference between the previous marker
    // on the read and this marker.
    // For the first marker in a read, this equals
    // the position of the first base of the marker in the read.
    using Shift = uint16_t;
    Shift shift;
};



class ChanZuckerberg::Nanopore2::Marker0 {
public:
    KmerId kmerId;

    // Position in the read of the first base of the marker.
    uint32_t position;

    // Ordinal of this marker in its read.
    // That is, the leftmost marker has ordinal=0,
    // the next 1, and so on.
    uint32_t ordinal;
};



// Class used to order markers by kmer id.
class ChanZuckerberg::Nanopore2::OrderMarkers0ByKmerId {
public:
    bool operator()(
        const Marker0& x,
        const Marker0& y) const
    {
        return x.kmerId < y.kmerId;
    }
};


#endif
