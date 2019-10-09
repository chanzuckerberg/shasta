#ifndef SHASTA_LOW_HASH_NEW_HPP
#define SHASTA_LOW_HASH_NEW_HPP

// Shasta
#include "Kmer.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadFlags.hpp"

namespace shasta {
    class LowHashNew;
    class CompressedMarker;
    class OrientedReadPair;
}


// This class uses the LowHash algorithm to find candidate pairs of aligned reads.
// It uses as features sequences of m consecutive markers.
// This is the new version that also stores alignmentCandidates.featureOrdinals
class shasta::LowHashNew :
    public MultithreadedObject<LowHashNew> {
public:

    // The constructor does all the work.
    LowHashNew(
        size_t m,                       // Number of consecutive markers that define a feature.
        double hashFraction,
        size_t minHashIterationCount,   // Number of minHash iterations.
        size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
        size_t maxBucketSize,           // The maximum size for a bucket to be used.
        size_t minFrequency,            // Minimum number of minHash hits for a pair to be considered a candidate.
        size_t threadCount,
        const MemoryMapped::Vector<KmerInfo>& kmerTable,
        const MemoryMapped::Vector<ReadFlags>& readFlags,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>&,
        MemoryMapped::Vector<OrientedReadPair>& candidates,
        MemoryMapped::VectorOfVectors< array<uint32_t, 2>, uint64_t>& featureOrdinals,
        const string& largeDataFileNamePrefix,
        size_t largeDataPageSize
    );

private:

    // Store some of the arguments passed to the constructor.
    size_t m;                       // Number of consecutive markers that define a feature.
    double hashFraction;
    size_t maxBucketSize;           // The maximum size for a bucket to be used.
    size_t minFrequency;            // Minimum number of minHash hits for a pair to be considered a candidate.
    size_t threadCount;
    const MemoryMapped::Vector<KmerInfo>& kmerTable;
    const MemoryMapped::Vector<ReadFlags>& readFlags;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const string& largeDataFileNamePrefix;
    size_t largeDataPageSize;
};

#endif
