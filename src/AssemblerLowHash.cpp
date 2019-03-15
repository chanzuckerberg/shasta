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
        readFlags,
        markers,
        alignmentCandidates,
        largeDataFileNamePrefix,
        largeDataPageSize);
}
