#ifndef CZI_SHASTA_LOW_HASH_HPP
#define CZI_SHASTA_LOW_HASH_HPP

// Shasta
#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultitreadedObject.hpp"
#include "OrientedReadPair.hpp"
#include "ReadId.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class LowHash;
    }
}



// This class uses the LowHash algorithm to find candidate pairs
// of aligned reads. It uses as features
// sequences of m consecutive markers.
class ChanZuckerberg::shasta::LowHash :
    public MultithreadedObject<LowHash>{
public:

    // The constructor does all the work.
    LowHash(
        size_t m,                       // Number of consecutive markers that define a feature.
        double hashFraction,
        size_t minHashIterationCount,   // Number of minHash iterations.
        size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
        size_t maxBucketSize,           // The maximum size for a bucket to be used.
        size_t minFrequency,            // Minimum number of minHash hits for a pair to be considered a candidate.
        size_t threadCount,
        const MemoryMapped::Vector<KmerInfo>& kmerTable,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>&,
        MemoryMapped::Vector<OrientedReadPair>&,
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

    // The current MinHash iteration.
    // This is used to compute a different MurmurHash function
    // at each iteration.
    size_t iteration;

    // The low hashes of each oriented read, at the current LowHash iteration.
    // Indexed by OrientedReadId::getValue().
    uint64_t hashThreshold;
    vector< vector<uint64_t> > lowHashes;
    void computeLowHashes(size_t threadId);

    // The mask used to compute to compute the bucket
    // corresponding to a hash value.
    uint64_t mask;

    // The buckets containing oriented read ids.
    MemoryMapped::VectorOfVectors<OrientedReadId::Int, uint64_t> buckets;

    // Data structure used to store candidate pairs.
    // Indexed by readId0, the read id of the lower numbered
    // read in the pair.
    // For each readId0, this is sorted by (readId1,isSameStrand) (with readId1 > readId0).
    class Candidate {
    public:
        ReadId readId1;      // The higher numbered read in the pair, readId1 > readId0.
        uint16_t frequency;  // Number of times this pair was found during MinHash.
        bool isSameStrand;   // True if the two reads are on the same strand.
        Candidate(ReadId readId1, uint16_t frequency, bool isSameStrand) :
            readId1(readId1), frequency(frequency), isSameStrand(isSameStrand) {}
    };
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



    // Thread functions.

    // Pass1: compute the low hashes for each oriented read
    // and prepare the buckets for filling.
    void pass1ThreadFunction(size_t threadId);

    // Pass 2: fill the buckets.
    void pass2ThreadFunction(size_t threadId);

    // Pass 3: sort the buckets.
    void pass3ThreadFunction(size_t threadId);

    // Pass 4: inspect the buckets to find candidates.
    void pass4ThreadFunction(size_t threadId);

};

#endif
