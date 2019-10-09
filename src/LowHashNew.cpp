// Shasta.
#include "LowHashNew.hpp"
using namespace shasta;



LowHashNew::LowHashNew(
    size_t m,                       // Number of consecutive markers that define a feature.
    double hashFraction,
    size_t minHashIterationCount,   // Number of minHash iterations.
    size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
    size_t maxBucketSize,           // The maximum size for a bucket to be used.
    size_t minFrequency,            // Minimum number of minHash hits for a pair to be considered a candidate.
    size_t threadCountArgument,
    const MemoryMapped::Vector<KmerInfo>& kmerTable,
    const MemoryMapped::Vector<ReadFlags>& readFlags,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    MemoryMapped::Vector<OrientedReadPair>& candidateAlignments,
    MemoryMapped::VectorOfVectors< array<uint32_t, 2>, uint64_t>& featureOrdinals,
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize
    ) :
    MultithreadedObject(*this),
    m(m),
    hashFraction(hashFraction),
    maxBucketSize(maxBucketSize),
    minFrequency(minFrequency),
    threadCount(threadCountArgument),
    kmerTable(kmerTable),
    readFlags(readFlags),
    markers(markers),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize)

{
    SHASTA_ASSERT(0);
}

