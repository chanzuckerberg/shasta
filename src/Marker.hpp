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

        // Classes that will be used to represent markers
        // when the restructuring of marker storage is complete.
        class CompressedMarker;
        class Marker;
        class MarkerWithOrdinal;
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



// Markers in shared memory are stored using class CompressedMarker
// which requires only 5 bytes per marker.

// For a run with 120 Gb of coverage and 10% of k-mers
// used as markers, storing all the 24 G markers requires
// 120 GB (we store markers for each read on both strands).
// This compares with 30 GB to store the reads
// (we store reads on one strand only).

// This layout results in unaligned memory accesses.
// This is not a problem as modern processors (beginning with Nehalem)
// have a much lower performance penalty for unaligned memory access
// than older processors did:
// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.455.4198&rep=rep1&type=pdf

class ChanZuckerberg::Nanopore2::CompressedMarker {
public:

    // The id of the k-mer for this marker.
    KmerId kmerId __attribute__ ((packed));

    // The position of this marker in the oriented read.
    // This limits the length of a read to 2^24=16Mib bases.
    Uint24 position;

};
static_assert(sizeof(ChanZuckerberg::Nanopore2::CompressedMarker) == 5,
    "Unexpected size of class CompressedMarker.");



// This stores the same information as CompressedMarker,
// but using built-in, aligned integers.
class ChanZuckerberg::Nanopore2::Marker {
public:

    // The id of the k-mer for this marker.
    KmerId kmerId;

    // The position of this marker in the oriented read.
    uint32_t position;

    // Constructor from a CompressedMarker.
    Marker(const CompressedMarker& compressedMarker) :
        kmerId(compressedMarker.kmerId),
        position(compressedMarker.position)
    {}
};



// This also stores the ordinal, that is the index
// of the marker in the oriented read, when the markers
// are sorted by position in the read.
class ChanZuckerberg::Nanopore2::MarkerWithOrdinal : public Marker {
public:
    uint32_t ordinal;

    // Constructor from a marker and an ordinal.
    MarkerWithOrdinal(const Marker& marker, uint32_t ordinal) :
        Marker(marker),
        ordinal(ordinal)
    {}
};



#endif
