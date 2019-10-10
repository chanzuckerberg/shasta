#ifndef SHASTA_LOW_HASH_NEW_HPP
#define SHASTA_LOW_HASH_NEW_HPP

// Shasta
#include "Kmer.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadFlags.hpp"

// Standard library.
#include "memory.hpp"

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

    // Vectors containing only the k-mer ids of all markers
    // for all oriented reads.
    // Indexed by OrientedReadId.getValue().
    // This is used to speed up the computation of hash functions.
    MemoryMapped::VectorOfVectors<KmerId, uint64_t> kmerIds;
    void createKmerIds();
    void createKmerIds(size_t threadId);

    // The mask used to compute to compute the bucket
    // corresponding to a hash value.
    uint64_t mask;

    // The threshold for a hash value to be considered low.
    uint64_t hashThreshold;

    // The current MinHash iteration.
    // This is used to compute a different MurmurHash function
    // at each iteration.
    size_t iteration;

    // The low hashes of each oriented read and the ordinals at
    // which the corresponding feature occurs.
    // This is recomputed at each iteration.
    // Indexed by OrientedReadId::getValue().
    vector< vector< pair<uint64_t, uint32_t> > > lowHashes;
    void computeLowHashes(size_t threadId);

    // Each bucket entry describes a low hash feature.
    // It consists of an oriented read id and
    // the ordinal where the low hash feature appears.
    class BucketEntry {
    public:
        OrientedReadId orientedReadId;
        uint32_t ordinal;
        BucketEntry(
            OrientedReadId orientedReadId,
            uint32_t ordinal) :
            orientedReadId(orientedReadId),
            ordinal(ordinal) {}
        BucketEntry() {}
    };
    MemoryMapped::VectorOfVectors<BucketEntry, uint64_t> buckets;


    // Compute a histogram of the number of entries in each histogram.
    void computeBucketHistogram();
    void computeBucketHistogramThreadFunction(size_t threadId);
    vector< vector<uint64_t> > threadBucketHistogram;
    ofstream histogramCsv;


    // When two oriented reads appear in the same bucket, we
    // check if that happens by chance or because we found a
    // common feature between the two oriented reads.
    // In the latter case, we store a new CommonFeature
    // containing the two OrientedReadId's and
    // the ordinals where the feature appears.
    // Each thread stores into its own vector of common features.
    class CommonFeature {
    public:
        array<OrientedReadId, 2> orientedReadIds;
        array<uint32_t, 2> ordinals;
        CommonFeature() {}
        CommonFeature(
            OrientedReadId orientedReadId0,
            OrientedReadId orientedReadId1,
            uint32_t ordinal0,
            uint32_t ordinal1
            ) :
            orientedReadIds({orientedReadId0, orientedReadId1}),
            ordinals({ordinal0, ordinal1})
        {}
    };
    vector< shared_ptr<MemoryMapped::Vector<CommonFeature> > > threadCommonFeatures;
    uint64_t countTotalThreadCommonFeatures() const;



    // Thread functions.

    // Pass1: compute the low hashes for each oriented read
    // and count the number of entries in each bucket.
    void pass1ThreadFunction(size_t threadId);

    // Pass 2: fill the buckets.
    void pass2ThreadFunction(size_t threadId);

    // Pass 3: scan the buckets to find common features.
    void pass3ThreadFunction(size_t threadId);
};

#endif
