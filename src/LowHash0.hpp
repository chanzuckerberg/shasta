#ifndef SHASTA_LOW_HASH0_HPP
#define SHASTA_LOW_HASH0_HPP

// Shasta
#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "OrientedReadPair.hpp"
#include "Reads.hpp"

// Standard library.
#include "fstream.hpp"
#include "tuple.hpp"

namespace shasta {
    class LowHash0;
    class Reads;
}



// This class uses the LowHash0 algorithm to find candidate pairs
// of aligned reads. It uses as features
// sequences of m consecutive markers.
class shasta::LowHash0 :
    public MultithreadedObject<LowHash0>{
public:

    // The constructor does all the work.
    LowHash0(
        size_t m,                       // Number of consecutive markers that define a feature.
        double hashFraction,

        // Iteration control. See MinHashOptions for details.
        size_t minHashIterationCount,
        double alignmentCandidatesPerRead,

        size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
        size_t minBucketSize,           // The minimum size for a bucket to be used.
        size_t maxBucketSize,           // The maximum size for a bucket to be used.
        size_t minFrequency,            // Minimum number of minHash hits for a pair to be considered a candidate.
        size_t threadCount,
        const MemoryMapped::Vector<KmerInfo>& kmerTable,
        const Reads& reads,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>&,
        MemoryMapped::Vector<OrientedReadPair>&,
        MemoryMapped::Vector< array<uint64_t, 3> >& readLowHashStatistics,
        const string& largeDataFileNamePrefix,
        size_t largeDataPageSize
);

private:

    // Store some of the arguments passed to the constructor.
    size_t m;                       // Number of consecutive markers that define a feature.
    double hashFraction;
    size_t minBucketSize;           // The minimum size for a bucket to be used.
    size_t maxBucketSize;           // The maximum size for a bucket to be used.
    size_t minFrequency;            // Minimum number of minHash hits for a pair to be considered a candidate.
    size_t threadCount;
    const MemoryMapped::Vector<KmerInfo>& kmerTable;
    const Reads& reads;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    MemoryMapped::Vector< array<uint64_t, 3> > &readLowHashStatistics;
    const string& largeDataFileNamePrefix;
    size_t largeDataPageSize;

    // Vectors containing only the k-mer ids of all markers
    // for all oriented reads.
    // Indexed by OrientedReadId.getValue().
    // This is used to speed up the computation of hash functions.
    MemoryMapped::VectorOfVectors<KmerId, uint64_t> kmerIds;
    void createKmerIds();
    void createKmerIds(size_t threadId);

    // The current MinHash iteration.
    // This is used to compute a different MurmurHash function
    // at each iteration.
    size_t iteration;

    // The low hashes of each oriented read, at the current LowHash0 iteration.
    // Indexed by OrientedReadId::getValue().
    uint64_t hashThreshold;
    vector< vector<uint64_t> > lowHashes;
    void computeLowHashes(size_t threadId);

    // The mask used to compute to compute the bucket
    // corresponding to a hash value.
    uint64_t mask;

    // Each bucket entry consists of an oriented read id and
    // the 32 most significant bits of the hash value, used
    // for collision avoidance.
    class BucketEntry {
    public:
        OrientedReadId orientedReadId;
        uint32_t hashHighBits;
        BucketEntry(
            OrientedReadId orientedReadId,
            uint64_t hash) :
            orientedReadId(orientedReadId),
            hashHighBits(uint32_t(hash >> 32)) {}
        BucketEntry() {}
    };
    MemoryMapped::VectorOfVectors<BucketEntry, uint64_t> buckets;



    // Class used to store candidate pairs.
    class Candidate {
    public:
        ReadId readId1;             // The higher numbered read in the pair, readId1 > readId0.
        uint8_t strand;             // 0=same strand, 1=opposite strands.
        uint16_t frequency;         // Number of times this pair was found during LowHash0.

        // Create a new candidate with frequency 1.
        Candidate(
            ReadId readId1,
            Strand strand) :
            readId1(readId1), strand(uint8_t(strand)), frequency(1) {}

        Candidate() {}

        // Comparison operators.
        bool operator==(const Candidate& that) const
        {
            return (readId1 == that.readId1) && (strand == that.strand);
        }
        bool operator<(const Candidate& that) const
        {
            return tie(readId1, strand) < tie(that.readId1, that.strand);
        }
    };
    static_assert(sizeof(Candidate) == 8, "Unexpected size of LowHash0::Candidate.");



    // Merge two sorted vectors of candidates into a third one.
    // The input vectors can be sorted but can have duplicates.
    // During merging, when two candidates with the same readId1
    // and strand are found, they are combined, adding up their frequency.
    // This is used by pass4ThreadFunction.
    static void merge(
        const vector<Candidate>&,
        const vector<Candidate>&,
        vector<Candidate>&
        );



    // The alignment candidates for each read.
    // Indexed by readId0, the read id of the lower numbered read in the pair.
    // We only store pairs with readId1>readId0.
    // For each readId0, this is kept sorted.
    vector< vector<Candidate> > candidates;



    // Per-iteration statistics for each thread.
    class ThreadStatistics {
    public:
        uint64_t highFrequency;
        uint64_t total;
        uint64_t capacity;
        ThreadStatistics()
        {
            clear();
        }
        void clear()
        {
            highFrequency = 0;
            total = 0;
            capacity = 0;
        }

    };
    vector<ThreadStatistics> threadStatistics;



    // Compute a histogram of the number of entries in each histogram.
    void computeBucketHistogram();
    void computeBucketHistogramThreadFunction(size_t threadId);
    vector< vector<uint64_t> > threadBucketHistogram;
    ofstream histogramCsv;



    // Thread functions.

    // Pass1: compute the low hashes for each oriented read
    // and prepare the buckets for filling.
    void pass1ThreadFunction(size_t threadId);

    // Pass 2: fill the buckets.
    void pass2ThreadFunction(size_t threadId);

    // Pass 3: inspect the buckets to find candidates.
    void pass3ThreadFunction(size_t threadId);

};

#endif
