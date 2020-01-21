// Shasta.
#include "LowHash1.hpp"
#include "AlignmentCandidates.hpp"
#include "Marker.hpp"
using namespace shasta;

// Standad library.
#include "algorithm.hpp"
#include "chrono.hpp"



LowHash1::LowHash1(
    size_t m,                       // Number of consecutive markers that define a feature.
    double hashFraction,
    size_t minHashIterationCount,   // Number of minHash iterations.
    size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
    size_t minBucketSize,           // The minimum size for a bucket to be used.
    size_t maxBucketSize,           // The maximum size for a bucket to be used.
    size_t minFrequency,            // Minimum number of minHash hits for a pair to be considered a candidate.
    size_t threadCountArgument,
    const MemoryMapped::Vector<KmerInfo>& kmerTable,
    const MemoryMapped::Vector<ReadFlags>& readFlags,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    AlignmentCandidates& candidates,
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize
    ) :
    MultithreadedObject(*this),
    m(m),
    hashFraction(hashFraction),
    minBucketSize(minBucketSize),
    maxBucketSize(maxBucketSize),
    minFrequency(minFrequency),
    threadCount(threadCountArgument),
    kmerTable(kmerTable),
    readFlags(readFlags),
    markers(markers),
    candidates(candidates),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize),
    histogramCsv("LowHashBucketHistogram.csv")

{
    cout << timestamp << "LowHash1 begins." << endl;
    const auto tBegin = steady_clock::now();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Estimate the total number of low hashes and its base 2 log.
    // Except for very short reads, each marker generates a feature,
    // and each feature generates a low hash with probability hashFraction.
    // So an estimate of the total number of hashes is:
    const uint64_t totalLowHashCountEstimate =
        uint64_t(hashFraction * double(markers.totalSize()));
    const uint64_t leadingZeroBitCount = __builtin_clzl(totalLowHashCountEstimate);
    const uint64_t log2TotalLowHashCountEstimate = 64 - leadingZeroBitCount;

    // If log2MinHashBucketCount is 0, choose a reasonable value
    // for the current number of reads.
    // Otherwise, check that log2MinHashBucketCount is not unreasonably small.
    if(log2MinHashBucketCount == 0) {
        log2MinHashBucketCount = 5 + log2TotalLowHashCountEstimate;
    } else {
        if(log2MinHashBucketCount < log2TotalLowHashCountEstimate) {
            throw runtime_error("LowHash1: log2MinHashBucketCount is unreasonably small.");
        }
    }

    // Set the number of buckets and the corresponding mask.
    const uint64_t bucketCount = 1ULL << log2MinHashBucketCount;
    mask = bucketCount - 1;
    cout << "LowHash1 algorithm will use 2^" << log2MinHashBucketCount;
    cout << " = " << bucketCount << " buckets. "<< endl;
    cout << "Estimated number of low hashes per iteration " << totalLowHashCountEstimate << endl;
    cout << "Estimated load factor " << double(totalLowHashCountEstimate)/double(bucketCount) << endl;

    // Create vectors containing only the k-mer ids of all markers.
    // This is used to speed up the computation of hash functions.
    cout << timestamp << "Creating kmer ids for oriented reads." << endl;
    createKmerIds();

    // Compute the threshold for a hash value to be considered low.
    hashThreshold = uint64_t(hashFraction * double(std::numeric_limits<uint64_t>::max()));

    // The number of oriented reads, each with its own vector of markers.
    const OrientedReadId::Int orientedReadCount = OrientedReadId::Int(markers.size());
    const ReadId readCount = orientedReadCount / 2;
    SHASTA_ASSERT(orientedReadCount == 2*readCount);

    // Set up work areas.
    buckets.createNew(
            largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "tmp-LowHash-Buckets"),
            largeDataPageSize);
    lowHashes.resize(orientedReadCount);
    threadCommonFeatures.resize(threadCount);
    for(size_t threadId=0; threadId!=threadCount; threadId++) {
        threadCommonFeatures[threadId] = make_shared<MemoryMapped::Vector<CommonFeature> >();
        threadCommonFeatures[threadId]->createNew(
            largeDataFileNamePrefix.empty() ? "" :
            (largeDataFileNamePrefix + "tmp-LowHash-ThreadCommonFeatures-" + to_string(threadId)),
            largeDataPageSize);
    }

    // Write the header of the histogram file.
    histogramCsv << "Iteration,BucketSize,BucketCount,FeatureCount\n";

    // LowHash iteration loop.
    for(iteration=0; iteration<minHashIterationCount; iteration++) {
        cout << timestamp << "LowHash iteration " << iteration << " begins." << endl;

        // Compute the low hashes for each oriented read
        // and count the number of low hash features in each bucket.
        buckets.clear();
        buckets.beginPass1(bucketCount);
        size_t batchSize = 10000;
        setupLoadBalancing(readCount, batchSize);
        runThreads(&LowHash1::computeHashesThreadFunction, threadCount);

        // Fill the buckets.
        buckets.beginPass2();
        setupLoadBalancing(readCount, batchSize);
        runThreads(&LowHash1::fillBucketsThreadFunction, threadCount);
        buckets.endPass2(false, false);
        cout << "Load factor at this iteration " <<
            double(buckets.totalSize()) / double(buckets.size()) << endl;
        computeBucketHistogram();

        // Scan the buckets to find common features.
        // Each thread stores the common features it finds in its own vector.
        const uint64_t oldCommonFeatureCount = countTotalThreadCommonFeatures();
        batchSize = 10000;
        setupLoadBalancing(bucketCount, batchSize);
        runThreads(&LowHash1::scanBucketsThreadFunction, threadCount);
        const uint64_t newCommonFeatureCount = countTotalThreadCommonFeatures();
        cout << "Stored " << newCommonFeatureCount-oldCommonFeatureCount <<
            " common features at this iteration." << endl;
    }

    // Gather together all the common features found by all threads.
    cout << timestamp << "Gathering common features found by all threads." << endl;
    gatherCommonFeatures();
    cout << timestamp << "Total number of common features including duplicates is " <<
        commonFeatures.totalSize() << endl;

    // We no longer need the common features by thread.
    for(size_t threadId=0; threadId!=threadCount; threadId++) {
        threadCommonFeatures[threadId]->remove();
        threadCommonFeatures[threadId] = 0;
    }
    threadCommonFeatures.clear();

    // Process the common features.
    // For each orientedReadId0, we look at all the CommonFeatureInfo we have
    // and sort them by orientedReadId1, then by ordinals, and remove duplicates.
    // We then find groups of at least minFrequency common features involving the
    // same pair(orientedReadId0, orientedReadId1)
    cout << timestamp << "Processing the common features we found." << endl;
    processCommonFeatures();

    // Clean up.
    buckets.remove();
    kmerIds.remove();
    lowHashes.clear();
    commonFeatures.remove();

    // Done.
    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "LowHash1 completed in " << tTotal << " s." << endl;
}



void LowHash1::createKmerIds()
{
    kmerIds.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "tmp-LowHash-Markers"),
        largeDataPageSize);
    const ReadId orientedReadCount = ReadId(markers.size());
    const ReadId readCount = orientedReadCount / 2;
    kmerIds.beginPass1(orientedReadCount);
    for(ReadId readId=0; readId!=readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto markerCount = markers.size(orientedReadId.getValue());
            kmerIds.incrementCount(orientedReadId.getValue(), markerCount);
        }
    }
    kmerIds.beginPass2();
    kmerIds.endPass2(false);
    const size_t batchSize = 10000;
    setupLoadBalancing(readCount, batchSize);
    runThreads(&LowHash1::createKmerIds, threadCount);
}



// Thread function for createKmerIds.
void LowHash1::createKmerIds(size_t threadId)
{

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];

                SHASTA_ASSERT(kmerIds.size(orientedReadId.getValue()) == orientedReadMarkers.size());

                auto pointer = kmerIds.begin(orientedReadId.getValue());
                for(const CompressedMarker& marker: orientedReadMarkers) {
                    *pointer++ = marker.kmerId;
                }
            }
        }
    }
}



// Thread function to compute the low hashes for each oriented read
// and count the number of entries in each bucket.
void LowHash1::computeHashesThreadFunction(size_t threadId)
{
    const int featureByteCount = int(m * sizeof(KmerId));
    const uint64_t seed = iteration * 37;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {
            if(readFlags[readId].isPalindromic) {
                continue;
            }
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);

                vector< pair<uint64_t, uint32_t> >& orientedReadLowHashes = lowHashes[orientedReadId.getValue()];
                orientedReadLowHashes.clear();
                const size_t markerCount = kmerIds.size(orientedReadId.getValue());

                // Handle the pathological case where there are fewer than m markers.
                // This oriented read ends up in no bucket.
                if(markerCount < m) {
                    continue;
                }

                // Get the markers for this oriented read.
                KmerId* kmerIdsPointer = kmerIds.begin(orientedReadId.getValue());
                const size_t featureCount = markerCount - m + 1;

                // Loop over features of this oriented read.
                // Features are sequences of m consecutive markers.
                for(size_t j=0; j<featureCount; j++, kmerIdsPointer++) {
                    const uint64_t hash = MurmurHash64A(kmerIdsPointer, featureByteCount, seed);
                    if(hash < hashThreshold) {
                        orientedReadLowHashes.push_back(make_pair(hash, j));
                        const uint64_t bucketId = hash & mask;
                        buckets.incrementCountMultithreaded(bucketId);
                    }
                }
            }
        }
    }

}



// Thread function to fill the buckets.
void LowHash1::fillBucketsThreadFunction(size_t threadId)
{

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {
            if(readFlags[readId].isPalindromic) {
                continue;
            }
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);
                const vector< pair<uint64_t, uint32_t> > & orientedReadLowHashes = lowHashes[orientedReadId.getValue()];

                for(const auto& p: orientedReadLowHashes) {
                    const uint64_t hash = p.first;
                    const uint64_t bucketId = hash & mask;
                    const uint32_t ordinal = p.second;
                    buckets.storeMultithreaded(bucketId, BucketEntry(orientedReadId, ordinal));
                }
            }
        }
    }
}



void LowHash1::computeBucketHistogram()
{
    threadBucketHistogram.clear();
    threadBucketHistogram.resize(threadCount);
    const uint64_t batchSize = 10000;
    setupLoadBalancing(buckets.size(), batchSize);
    runThreads(&LowHash1::computeBucketHistogramThreadFunction, threadCount);

    // Combine the histograms found by each thread.
    uint64_t largestBucketSize = 0;
    for(const vector<uint64_t>& histogram: threadBucketHistogram) {
        largestBucketSize = max(largestBucketSize, uint64_t(histogram.size()));
    }
    vector<uint64_t> bucketHistogram(largestBucketSize, 0);
    for(const vector<uint64_t>& histogram: threadBucketHistogram) {
        for(uint64_t bucketSize=0; bucketSize<histogram.size(); bucketSize++) {
            bucketHistogram[bucketSize] += histogram[bucketSize];
        }
    }

    for(uint64_t bucketSize=0; bucketSize<bucketHistogram.size(); bucketSize++) {
        const uint64_t frequency = bucketHistogram[bucketSize];
        if(frequency) {
            histogramCsv <<
                iteration << "," <<
                bucketSize << "," <<
                frequency << "," <<
                bucketSize*frequency << "\n";
        }
    }


}
void LowHash1::computeBucketHistogramThreadFunction(size_t threadId)
{
    vector<uint64_t>& histogram = threadBucketHistogram[threadId];
    histogram.clear();
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t bucketId=begin; bucketId!=end; bucketId++) {
            const uint64_t bucketSize = buckets.size(bucketId);
            if(bucketSize >= histogram.size()) {
                histogram.resize(bucketSize + 1, 0);
            }
            ++histogram[bucketSize];
        }
    }
}



// Thread function to scan the buckets to find common features.
void LowHash1::scanBucketsThreadFunction(size_t threadId)
{
    // Access the vector where this thread will store
    // the common features it finds.
    MemoryMapped::Vector<CommonFeature>& commonFeatures = *threadCommonFeatures[threadId];

    const uint64_t mLocal = uint64_t(m);

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over buckets in this batch.
        for(uint64_t bucketId=begin; bucketId!=end; bucketId++) {

            // Access this bucket.
            const span<BucketEntry> bucket = buckets[bucketId];
            if(bucket.size() < max(size_t(2), minBucketSize)) {
                continue;
            }
            if(bucket.size() > maxBucketSize) {
                continue;
            }

            // Loop over pairs of bucket entries.
            for(const BucketEntry& feature0: bucket) {
                const OrientedReadId orientedReadId0 = feature0.orientedReadId;
                const ReadId readId0 = orientedReadId0.getReadId();
                const Strand strand0 = orientedReadId0.getStrand();
                const uint32_t ordinal0 = feature0.ordinal;
                const auto allKmerIds0 = kmerIds[orientedReadId0.getValue()];
                const auto featureKmerIds0 = allKmerIds0.begin() + ordinal0;
                const uint32_t markerCount0 = uint32_t(allKmerIds0.size());

                for(const BucketEntry& feature1: bucket) {
                    const OrientedReadId orientedReadId1 = feature1.orientedReadId;
                    const ReadId readId1 = orientedReadId1.getReadId();

                    // Only consider the ones where readId0 < readId1.
                    if(readId0 >= readId1) {
                        continue;
                    }

                    const Strand strand1 = orientedReadId1.getStrand();
                    const uint32_t ordinal1 = feature1.ordinal;
                    const auto allKmerIds1 = kmerIds[orientedReadId1.getValue()];
                    const auto featureKmerIds1 = allKmerIds1.begin() + ordinal1;
                    const uint32_t markerCount1 = uint32_t(allKmerIds1.size());

                    // If the k-mers are not the same, this is a collision. Discard.
                    if(not std::equal(featureKmerIds0, featureKmerIds0+mLocal, featureKmerIds1)) {
                        continue;
                    }

                    // We found a common feature. Store it.
                    // If read0 is on strand 1, we have to reverse the ordinals.
                    if(strand0 == 0) {
                        commonFeatures.push_back(CommonFeature(
                            readId0,
                            readId1,
                            strand0==strand1,
                            ordinal0,
                            ordinal1));
                    } else {
                        commonFeatures.push_back(CommonFeature(
                            readId0,
                            readId1,
                            strand0==strand1,
                            markerCount0-1-ordinal0,
                            markerCount1-1-ordinal1));
                    }
                }
            }
        }
    }
}


// Add up the number of common feature found by all threads.
uint64_t LowHash1::countTotalThreadCommonFeatures() const
{
    uint64_t n = 0;
    for(const auto& v: threadCommonFeatures) {
        n += v->size();
    }
    return n;
}



void LowHash1::gatherCommonFeatures()
{
    commonFeatures.createNew(
            largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "tmp-CommonFeatures"),
            largeDataPageSize);
    commonFeatures.beginPass1(kmerIds.size()/2);
    runThreads(&LowHash1::gatherCommonFeaturesPass1, threadCount);
    commonFeatures.beginPass2();
    runThreads(&LowHash1::gatherCommonFeaturesPass2, threadCount);
    commonFeatures.endPass2(false);
}
void LowHash1::gatherCommonFeaturesPass1(size_t threadId)
{
    const MemoryMapped::Vector<CommonFeature>& v = *threadCommonFeatures[threadId];
    for(const CommonFeature& commonFeature: v) {
        commonFeatures.incrementCountMultithreaded(commonFeature.orientedReadPair.readIds[0]);
    }
}
void LowHash1::gatherCommonFeaturesPass2(size_t threadId)
{
    const MemoryMapped::Vector<CommonFeature>& v = *threadCommonFeatures[threadId];
    for(const CommonFeature& commonFeature: v) {
        commonFeatures.storeMultithreaded(
            commonFeature.orientedReadPair.readIds[0],
            CommonFeatureInfo(commonFeature));
    }
}



// Process the common features.
// For each readId0, we look at all the CommonFeatureInfo we have
// and sort them by readId1, then by ordinals, and remove duplicates.
// We then find groups of at least minFrequency common features involving the
// same pair(orientedReadId0, orientedReadId1)
// Each group generates an alignment candidate and the
// corresponding common features.
// Each thread stores the alignment candidates it finds in its own vector.
void LowHash1::processCommonFeatures()
{
    const uint64_t readCount = kmerIds.size() / 2;
    const uint64_t batchSize = 1000;

    // Prepare areas where each thread will store what it finds.
    threadCandidateTable.resize(readCount);
    threadAlignmentCandidates.resize(threadCount);
    threadCandidateHistogram.resize(threadCount);

    // Extract the candidates and features.
    setupLoadBalancing(readCount, batchSize);
    runThreads(&LowHash1::processCommonFeaturesThreadFunction, threadCount);



    // Gather the candidates and the features.
    for(ReadId readId0=0; readId0<readCount; readId0++) {

        // Figure out where the candidates are stored.
        const auto& info = threadCandidateTable[readId0];
        const uint64_t threadId = info[0];
        const uint64_t begin = info[1];
        const uint64_t end = info[2];

        // Loop over all these candidates.
        for(uint64_t i=begin; i!=end; ++i) {
            const OrientedReadPair& orientedReadPair =
                threadAlignmentCandidates[threadId]->candidates[i];
            SHASTA_ASSERT(orientedReadPair.readIds[0] == readId0);
            candidates.candidates.push_back(orientedReadPair);
            const auto features = threadAlignmentCandidates[threadId]->featureOrdinals[i];
            candidates.featureOrdinals.appendVector(features.begin(), features.end());
        }
    }
    SHASTA_ASSERT(candidates.candidates.size() == candidates.featureOrdinals.size());
    cout << timestamp << "Found " << candidates.candidates.size() <<
        " alignment candidates with a total " <<
        candidates.featureOrdinals.totalSize() <<
        " features." << endl;



    // Combine the histograms found by each thread.
    for(size_t threadId=0; threadId!=threadCount; threadId++) {
        const vector<uint64_t>& v = threadCandidateHistogram[threadId];
        for(uint64_t i=0; i<v.size(); i++){
            const uint64_t n = v[i];
            if(n > 0) {
                if(candidateHistogram.size() <= n){
                    candidateHistogram.resize(n+1, 0);
                }
                candidateHistogram[i] += n;
            }
        }
    }
    ofstream csv("LowHashCandidateHistogram.csv");
    csv << "CommonFeatureCount,Frequency\n";
    for(uint64_t i=0; i<candidateHistogram.size(); i++) {
        const uint64_t n = candidateHistogram[i];
        if(n > 0) {
            csv << i << "," << n << "\n";
        }
    }



    // Clean up.
    threadCandidateTable.clear();
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        threadAlignmentCandidates[threadId]->candidates.remove();
        threadAlignmentCandidates[threadId]->featureOrdinals.remove();
    }
    threadAlignmentCandidates.clear();
}



void LowHash1::processCommonFeaturesThreadFunction(size_t threadId)
{
    // Access the vector where this thread will store
    // the alignment candidates it finds.
    threadAlignmentCandidates[threadId] = make_shared<AlignmentCandidates>();
    AlignmentCandidates& alignmentCandidates = *threadAlignmentCandidates[threadId];
    alignmentCandidates.candidates.createNew(
        largeDataFileNamePrefix.empty() ? "" :
        (largeDataFileNamePrefix + "tmp-ThreadAlignmentCandidates-" + to_string(threadId)),
        largeDataPageSize);
    alignmentCandidates.featureOrdinals.createNew(
        largeDataFileNamePrefix.empty() ? "" :
        (largeDataFileNamePrefix + "tmp-ThreadAlignmentCandidatesOrdinals-" + to_string(threadId)),
        largeDataPageSize);
    vector<uint64_t>& histogram = threadCandidateHistogram[threadId];

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over ReadId's in this batch.
        for(ReadId readId0=ReadId(begin); readId0!=ReadId(end); readId0++) {
            // std::lock_guard<std::mutex> lock(mutex); // ************************** TAKE OUT!
            // cout << "Working on readId0 " << readId0 << endl;
            const span<CommonFeatureInfo> features = commonFeatures[readId0];
            threadCandidateTable[readId0][0] = uint64_t(threadId);
            threadCandidateTable[readId0][1] = alignmentCandidates.candidates.size();;

            /*
            cout << features.size() << " features before deduplication:" << endl;
            for(auto it=features.begin(); it!=features.end(); ++it) {
                const CommonFeatureInfo& feature = *it;
                cout <<
                    feature.readId1 << " " <<
                    (feature.isSameStrand ? "same strand " : " opposite strands ") <<
                    feature.ordinals[0] << " " <<
                    feature.ordinals[1] << " " <<
                    int32_t(feature.ordinals[1]) - int32_t(feature.ordinals[0]) << "\n";
            }
            */

            // Deduplicate.
            const auto uniqueBegin = features.begin();
            auto uniqueEnd = features.end();
            sort(uniqueBegin, uniqueEnd);
            uniqueEnd = unique(uniqueBegin, uniqueEnd);

            /*
            cout << uniqueEnd-uniqueBegin << " features after deduplication:" << endl;
            for(auto it=uniqueBegin; it!=uniqueEnd; ++it) {
                const CommonFeatureInfo& feature = *it;
                cout <<
                    feature.readId1 << " " <<
                    (feature.isSameStrand ? "same strand " : " opposite strands ") <<
                    feature.ordinals[0] << " " <<
                    feature.ordinals[1] << " " <<
                    int32_t(feature.ordinals[1]) - int32_t(feature.ordinals[0]) << "\n";
            }
            */

            // Loop over streaks of features with the same readId1 and isSameStrand.
            for(auto it=uniqueBegin; it!=uniqueEnd;) {
                auto streakBegin = it;
                auto streakEnd = streakBegin;
                const ReadId readId1 = streakBegin->readId1;
                const bool isSameStrand = streakBegin->isSameStrand;
                while(streakEnd!=uniqueEnd and streakEnd->readId1==readId1 and streakEnd->isSameStrand==isSameStrand) {
                    ++streakEnd;
                }

                // Increment the histogram.
                const int64_t streakLength = streakEnd - streakBegin;
                if(histogram.size() <= uint64_t(streakLength)) {
                    histogram.resize(streakLength + 1, 0);
                }
                ++histogram[streakLength];

                // If too few, skip.
                if(streakLength < int64_t(minFrequency)) {
                    it = streakEnd;
                    continue;
                }

                /*
                cout << "Common features of reads " <<
                    readId0 << " " <<
                    readId1 << (isSameStrand ? " same strand" : " opposite strands") << ":\n";
                for(auto it=streakBegin; it!=streakEnd; ++it) {
                    const CommonFeatureInfo& feature = *it;
                    cout <<
                        feature.ordinals[0] << " " <<
                        feature.ordinals[1] << " " <<
                        int32_t(feature.ordinals[1]) - int32_t(feature.ordinals[0]) << "\n";
                }
                cout << "Marker count " <<
                    kmerIds[OrientedReadId(readId0, 0).getValue()].size() << " " <<
                    kmerIds[OrientedReadId(readId1, 0).getValue()].size() << ":\n";
                */

                // This streak generates an alignment candidate
                // and the corresponding common features.
                alignmentCandidates.candidates.push_back(OrientedReadPair(readId0, readId1, isSameStrand));
                alignmentCandidates.featureOrdinals.appendVector();
                for(auto it=streakBegin; it!=streakEnd; ++it) {
                    const CommonFeatureInfo& feature = *it;
                    alignmentCandidates.featureOrdinals.append(feature.ordinals);
                }

                // Prepare for the next streak.
                it = streakEnd;
            }
            threadCandidateTable[readId0][2] = alignmentCandidates.candidates.size();;
        }
    }
}
