#include "Assembler.hpp"
#include "OverlapFinder.hpp"
using namespace ChanZuckerberg;
using namespace shasta;




// Use the minHash algorithm to find pairs of overlapping oriented reads.
// Use as features sequences of m consecutive special k-mers.
void Assembler::findOverlaps(
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

    // Check that log2MinHashBucketCount is not unreasonably small.
    if((1ULL << (log2MinHashBucketCount-3ULL)) < readCount) {
        throw runtime_error("log2MinHashBucketCount is unreasonably small.\n"
            "Must at least equal base 2 log of number of reads plus 3.");
    }

    // Create the overlaps.
    candidateAlignments.createNew(largeDataName("CandidateAlignments"), largeDataPageSize);

    // Call the OverlapFinder to do the MinHash computation.
    OverlapFinder overlapFinder(
        m,
        minHashIterationCount,
        log2MinHashBucketCount,
        maxBucketSize,
        minFrequency,
        threadCount,
        kmerTable,
        markers,
        candidateAlignments,
        largeDataFileNamePrefix,
        largeDataPageSize);
}



void Assembler::accessCandidateAlignments()
{
    candidateAlignments.accessExistingReadOnly(largeDataName("CandidateAlignments"));
}


void Assembler::checkCandidateAlignmentsAreOpen() const
{
    if(!candidateAlignments.isOpen) {
        throw runtime_error("Candidate alignments are not accessible.");
    }
}



// Write the reads that overlap a given read.
void Assembler::writeOverlappingReads(
    ReadId readId0,
    Strand strand0,
    const string& fileName)
{
    // Check that we have what we need.
    checkReadsAreOpen();
    checkCandidateAlignmentsAreOpen();



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


