#include "Assembler.hpp"
#include "LowHash.hpp"
#include "LowHashNew.hpp"
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
    SHASTA_ASSERT(readCount > 0);

    // Create the alignment candidates.
    alignmentCandidates.candidates.createNew(largeDataName("AlignmentCandidates"), largeDataPageSize);

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
        alignmentCandidates.candidates,
        largeDataFileNamePrefix,
        largeDataPageSize);
}



void Assembler::accessAlignmentCandidates()
{
    alignmentCandidates.candidates.accessExistingReadOnly(largeDataName("AlignmentCandidates"));
}


void Assembler::checkAlignmentCandidatesAreOpen() const
{
    if(!alignmentCandidates.candidates.isOpen) {
        throw runtime_error("Alignment candidates are not accessible.");
    }
}



vector<OrientedReadPair> Assembler::getAlignmentCandidates() const
{
    checkAlignmentCandidatesAreOpen();
    vector<OrientedReadPair> v;
    copy(
        alignmentCandidates.candidates.begin(),
        alignmentCandidates.candidates.end(),
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



// New version that also stores alignmentCandidates.featureOrdinals.
// This can be used to filter the alignment candidates.
void Assembler::findAlignmentCandidatesLowHashNew(
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
    SHASTA_ASSERT(readCount > 0);

    // Prepare storage.
    alignmentCandidates.candidates.createNew(
        largeDataName("AlignmentCandidates"), largeDataPageSize);
    alignmentCandidates.featureOrdinals.createNew(
        largeDataName("AlignmentCandidatesFeatureOrdinale"), largeDataPageSize);

    // Do the computation.
    LowHashNew lowHashNew(
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
        alignmentCandidates.candidates,
        alignmentCandidates.featureOrdinals,
        largeDataFileNamePrefix,
        largeDataPageSize);
}
