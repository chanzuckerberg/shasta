#include "Assembler.hpp"
#include "LowHash0.hpp"
#include "LowHash1.hpp"
using namespace shasta;



// Use the LowHash algorithm to find alignment candidates.
// Use as features sequences of m consecutive special k-mers.
void Assembler::findAlignmentCandidatesLowHash0(
    size_t m,                       // Number of consecutive k-mers that define a feature.
    double hashFraction,            // Low hash threshold.

    // Iteration control. See MinHashOptions for details.
    size_t minHashIterationCount,
    double alignmentCandidatesPerRead,

    size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for lowHash.
    size_t minBucketSize,           // The minimum size for a bucket to be used.
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
    readLowHashStatistics.createNew(largeDataName("ReadLowHashStatistics"), largeDataPageSize);

    // Run the LowHash computation to find candidate alignments.
    LowHash0 lowHash(
        m,
        hashFraction,
        minHashIterationCount,
        alignmentCandidatesPerRead,
        log2MinHashBucketCount,
        minBucketSize,
        maxBucketSize,
        minFrequency,
        threadCount,
        kmerTable,
        readFlags,
        markers,
        alignmentCandidates.candidates,
        readLowHashStatistics,
        largeDataFileNamePrefix,
        largeDataPageSize);

    alignmentCandidates.unreserve();
}



void Assembler::accessAlignmentCandidates()
{
    alignmentCandidates.candidates.accessExistingReadOnly(largeDataName("AlignmentCandidates"));
}

void Assembler::accessReadLowHashStatistics()
{
    readLowHashStatistics.accessExistingReadOnly(largeDataName("ReadLowHashStatistics"));
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
void Assembler::findAlignmentCandidatesLowHash1(
    size_t m,                       // Number of consecutive k-mers that define a feature.
    double hashFraction,            // Low hash threshold.
    size_t minHashIterationCount,   // Number of lowHash iterations.
    size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for lowHash.
    size_t minBucketSize,           // The minimum size for a bucket to be used.
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
    LowHash1 lowHash1(
        m,
        hashFraction,
        minHashIterationCount,
        log2MinHashBucketCount,
        minBucketSize,
        maxBucketSize,
        minFrequency,
        threadCount,
        kmerTable,
        readFlags,
        markers,
        alignmentCandidates,
        largeDataFileNamePrefix,
        largeDataPageSize);
    
    alignmentCandidates.unreserve();
}



void Assembler::writeAlignmentCandidates() const
{

    // Sanity checks.
    checkMarkersAreOpen();
    SHASTA_ASSERT(alignmentCandidates.candidates.isOpen);

    if(alignmentCandidates.featureOrdinals.isOpen()) {
        SHASTA_ASSERT(
            alignmentCandidates.candidates.size() ==
            alignmentCandidates.featureOrdinals.size());
    }


    // Write out the candidates.
    ofstream csv("AlignmentCandidates.csv");
    csv << "ReadId0,ReadId1,SameStrand,";
    if(alignmentCandidates.featureOrdinals.isOpen()) {
        csv << "FeatureCount,";
    }
    csv << "\n";

    for(uint64_t i=0; i<alignmentCandidates.candidates.size(); i++) {
        const OrientedReadPair& candidate = alignmentCandidates.candidates[i];
        csv << candidate.readIds[0] << ",";
        csv << candidate.readIds[1] << ",";
        csv << (candidate.isSameStrand ? "Yes" : "No") << ",";
        if(alignmentCandidates.featureOrdinals.isOpen()) {
            csv << alignmentCandidates.featureOrdinals.size(i) << ",";
        }
        csv << "\n";
    }



    // Write out the features.
    if(alignmentCandidates.featureOrdinals.isOpen()) {

        ofstream csv("AlignmentCandidatesFeatures.csv");
        csv << "ReadId0,ReadId1,SameStrand,FeatureCount,Offset,"
            "Ordinal0,Ordinal1,Ordinal0Reversed,Ordinal1Reversed,MarkerCount0,MarkerCount1,"
            "\n";

        for(uint64_t i=0; i<alignmentCandidates.candidates.size(); i++) {
            const OrientedReadPair& candidate = alignmentCandidates.candidates[i];

            // The features for this candidate.
            const auto features = alignmentCandidates.featureOrdinals[i];

            for(const auto& feature: features) {
                const ReadId readId0 = candidate.readIds[0];
                const ReadId readId1 = candidate.readIds[1];
                const uint32_t markerCount0 = uint32_t(markers.size(OrientedReadId(readId0, 0).getValue()));
                const uint32_t markerCount1 = uint32_t(markers.size(OrientedReadId(readId1, 0).getValue()));
                const uint32_t ordinal0 = feature[0];
                const uint32_t ordinal1 = feature[1];

                csv << readId0 << ",";
                csv << readId1 << ",";
                csv << (candidate.isSameStrand ? "Yes" : "No") << ",";
                csv << features.size() << ",";
                csv << int32_t(ordinal1)-int32_t(ordinal0) << ",";
                csv << ordinal0 << ",";
                csv << ordinal1 << ",";
                csv << markerCount0 - 1 - ordinal0 << ",";
                csv << markerCount1 - 1 - ordinal1 << ",";
                csv << markerCount0 << ",";
                csv << markerCount1 << ",";
                csv << "\n";
            }
        }
    }
}



// This marks all pairs as alignment candidates.
void Assembler::markAlignmentCandidatesAllPairs()
{
    // Create the alignment candidates.
    alignmentCandidates.candidates.createNew(largeDataName("AlignmentCandidates"), largeDataPageSize);

    // Add all pairs on both orientations.
    const ReadId n = readCount();
    for(ReadId r0=0; r0<n-1; r0++) {
        for(ReadId r1=r0+1; r1<n; r1++) {
            alignmentCandidates.candidates.push_back(OrientedReadPair(r0, r1, true));
            alignmentCandidates.candidates.push_back(OrientedReadPair(r0, r1, false));
        }
    }

    alignmentCandidates.unreserve();
    cout << "Marked all pairs of reads as alignment candidates on both orientations." << endl;
}
