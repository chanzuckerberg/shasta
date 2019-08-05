#include "Assembler.hpp"
#include "MinHash.hpp"
using namespace shasta;




// Use the minHash algorithm to find alignment candidates.
// Use as features sequences of m consecutive special k-mers.
void Assembler::findAlignmentCandidatesMinHash(
    size_t m,                       // Number of consecutive k-mers that define a feature.
    size_t minHashIterationCount,   // Number of minHash iterations.
    size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
    size_t maxBucketSize,           // The maximum size for a bucket to be used.
    size_t minFrequency,            // Minimum number of minHash hits for a pair to become a candidate.
    size_t threadCount
)
{
    checkKmersAreOpen();
    checkMarkersAreOpen();
    const ReadId readCount = ReadId(markers.size() / 2);
    SHASTA_ASSERT(readCount > 0);



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

    // Run the MinHash computation to find candidate alignments.
    MinHash minHash(
        m,
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



void Assembler::accessAlignmentCandidates()
{
    alignmentCandidates.accessExistingReadOnly(largeDataName("AlignmentCandidates"));
}


void Assembler::checkAlignmentCandidatesAreOpen() const
{
    if(!alignmentCandidates.isOpen) {
        throw runtime_error("Alignment candidates are not accessible.");
    }
}



vector<OrientedReadPair> Assembler::getAlignmentCandidates() const
{
    checkAlignmentCandidatesAreOpen();
    vector<OrientedReadPair> v;
    copy(
        alignmentCandidates.begin(),
        alignmentCandidates.end(),
        back_inserter(v));
    return v;
}



// Write the reads that overlap a given read.
void Assembler::writeOverlappingReads(
    ReadId readId0,
    Strand strand0,
    const string& fileName)
{
    // Check that we have what we need.
    checkReadsAreOpen();
    checkAlignmentCandidatesAreOpen();



    // Open the output file and write the oriented read we were given.
    ofstream file(fileName);
    const OrientedReadId orientedReadId0(readId0, strand0);
    writeOrientedRead(orientedReadId0, file);

    const uint64_t length0 = reads[orientedReadId0.getReadId()].baseCount;
    cout << "Reads overlapping " << orientedReadId0 << " length " << length0 << endl;

    // Loop over all overlaps involving this oriented read.
    for(const uint64_t i: alignmentTable[orientedReadId0.getValue()]) {
        const AlignmentData& ad = alignmentData[i];

        // Get the other oriented read involved in this overlap.
        const OrientedReadId orientedReadId1 = ad.getOther(orientedReadId0);

        // Write it out.
        const uint64_t length1 = reads[orientedReadId1.getReadId()].baseCount;
        cout << orientedReadId1 << " length " << length1 << endl;
        writeOrientedRead(orientedReadId1, file);
    }
    cout << "Found " << alignmentTable[orientedReadId0.getValue()].size();
    cout << " overlapping oriented reads." << endl;

}


