#include "Assembler.hpp"
#include "LowHash.hpp"
using namespace ChanZuckerberg;
using namespace shasta;




// Use the LowHash algorithm to find alignment candidates.
// Use as features sequences of m consecutive special k-mers.
void Assembler::findAlignmentCandidatesLowHash(
    size_t m,                       // Number of consecutive k-mers that define a feature.
    double hashFraction,            // Low hash threshold.
    size_t minHashIterationCount,   // Number of lowHash iterations.
    size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for lowHash.
    size_t maxBucketSize,           // The maximum size for a bucket to be used.
    size_t minFrequency,            // Minimum number of minHash hits for a pair to become a candidate.
    size_t threadCount)
{

    // Check that we have what we need.
    checkKmersAreOpen();
    checkMarkersAreOpen();
    const ReadId readCount = ReadId(markers.size() / 2);
    CZI_ASSERT(readCount > 0);



    // If log2MinHashBucketCount is 0, choose a reasonable value
    // for the current number of reads.
    if(log2MinHashBucketCount == 0) {

        // Compute an approximate base 2 log of the number of reads.
        static_assert(sizeof(readCount) == 4, "Unexpected readCount size.");
        const int leadingZeroBitCount = __builtin_clz(readCount);
        const int log2ReadCount = 32 - leadingZeroBitCount;

        // Make log2MinHashBucketCount reasonably larger
        // than the approximate base 2 log of the number of reads.
        log2MinHashBucketCount = 5 + log2ReadCount;

        cout << "Set log2MinHashBucketCount to " << log2MinHashBucketCount << endl;
    }



    // Check that log2MinHashBucketCount is not unreasonably small.
    if((1ULL << (log2MinHashBucketCount-3ULL)) < readCount) {
        throw runtime_error("log2MinHashBucketCount is unreasonably small.\n"
            "Must at least equal base 2 log of number of reads plus 3.");
    }



    // Create the alignment candidates.
    alignmentCandidates.createNew(largeDataName("AlignmentCandidates"), largeDataPageSize);

    // Run the LowHash computation to find candidate alignments.
    LowHash lowHash(
        m,
        hashFraction,
        minHashIterationCount,
        log2MinHashBucketCount,
        maxBucketSize,
        minFrequency,
        threadCount,
        kmerTable,
        markers,
        alignmentCandidates,
        largeDataFileNamePrefix,
        largeDataPageSize);
}
