#ifndef SHASTA_KMER_HPP
#define SHASTA_KMER_HPP

#include "ShortBaseSequence.hpp"
#include <limits>

namespace shasta {

    // Types used to represent a k-mer and a k-mer id.
    // These limit the maximum k-mer length that can be used.
    using Kmer = ShortBaseSequence16;
    using KmerId = uint32_t;

    // Check for consistency of these two types.
    static_assert(
        std::numeric_limits<KmerId>::digits == 2*Kmer::capacity,
        "Kmer and KmerId types are inconsistent.");

    class KmerInfo;
}



class shasta::KmerInfo {
public:

    // Frequency of this k-mer in input reads.
    // Only filled in if selectKmersBasedOnFrequency
    // is used.
    uint64_t frequency = 0;

    KmerId reverseComplementedKmerId;
    bool isMarker;
    bool isRleKmer;

    // Hash function of the KmerId, used for downsampling markers
    // for alignments using method 3.
    uint32_t hash;
};

#endif
