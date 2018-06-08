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

We need to find and store all markers in all reads.
For memory economy, markers on all reads are stored in
a compressed format (class CompressedMarker)
that takes only 4 bytes for each occurrence of a marker in a read
(assuming the Kmer type is ShortBaseSequence8).
This is achieved by storing position offsets
of markers instead of actual positions.

This results in reasonable memory usage in typical cases.
For example, if we randomly choose 10% of all k-mers as markers,
this will require an average of 4 bytes per 10 bases or 0.4 bytes per base.
By comparison, storing the reads requires 0.25 bytes per base.
If the read take 30 GB, the compressed markers take 48 GB.

We also use an uncompressed format for markers (class Marker).
This contains an absolute position in the read rather than
an offset, and it also contains the ordinal of the marker in the read.
This permits storing by KmerId while keeping ordinal and position information,
which is needed when computing alignments.

*******************************************************************************/

#include "Kmer.hpp"
#include "Uint.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {

        class CompressedMarker;
        class Marker;
        class OrderMarkersByKmerId;
    }
}



class ChanZuckerberg::Nanopore2::CompressedMarker {
public:
    KmerId kmerId;

    // Position difference between the previous marker
    // on the read and this marker.
    // For the first marker in a read, this equals
    // the position of the first base of the marker in the read.
    using Shift = uint16_t;
    Shift shift;
};



class ChanZuckerberg::Nanopore2::Marker {
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
class ChanZuckerberg::Nanopore2::OrderMarkersByKmerId {
public:
    bool operator()(
        const Marker& x,
        const Marker& y) const
    {
        return x.kmerId < y.kmerId;
    }
};


#endif
