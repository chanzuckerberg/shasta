// Shasta.
#include "LowHashNew.hpp"
#include "Marker.hpp"
using namespace shasta;

// Standad library.
#include "algorithm.hpp"
#include "chrono.hpp"



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
    largeDataPageSize(largeDataPageSize),
    histogramCsv("LowHashBucketHistogram.csv")

{
    cout << timestamp << "LowHashNew begins." << endl;
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
            throw runtime_error("LowHashNew: log2MinHashBucketCount is unreasonably small.");
        }
    }

    // Set the number of buckets and the corresponding mask.
    const uint64_t bucketCount = 1 << log2MinHashBucketCount;
    mask = bucketCount - 1;
    cout << "LowHashNew algorithm will use 2^" << log2MinHashBucketCount;
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
        runThreads(&LowHashNew::computeHashesThreadFunction, threadCount);

        // Fill the buckets.
        buckets.beginPass2();
        setupLoadBalancing(readCount, batchSize);
        runThreads(&LowHashNew::fillBucketsThreadFunction, threadCount);
        buckets.endPass2(false, false);
        cout << "Load factor at this iteration " <<
            double(buckets.totalSize()) / double(buckets.size()) << endl;
        computeBucketHistogram();

        // Scan the buckets to find common features.
        // Each thread stores the common features it finds in its own vector.
        const uint64_t oldCommonFeatureCount = countTotalThreadCommonFeatures();
        batchSize = 10000;
        setupLoadBalancing(bucketCount, batchSize);
        runThreads(&LowHashNew::scanBucketsThreadFunction, threadCount);
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
    cout << timestamp << "LowHashNew completed in " << tTotal << " s." << endl;

    SHASTA_ASSERT(0);
}



void LowHashNew::createKmerIds()
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
    runThreads(&LowHashNew::createKmerIds, threadCount);
}



// Thread function for createKmerIds.
void LowHashNew::createKmerIds(size_t threadId)
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
void LowHashNew::computeHashesThreadFunction(size_t threadId)
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
void LowHashNew::fillBucketsThreadFunction(size_t threadId)
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



void LowHashNew::computeBucketHistogram()
{
    threadBucketHistogram.clear();
    threadBucketHistogram.resize(threadCount);
    const uint64_t batchSize = 10000;
    setupLoadBalancing(buckets.size(), batchSize);
    runThreads(&LowHashNew::computeBucketHistogramThreadFunction, threadCount);

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
void LowHashNew::computeBucketHistogramThreadFunction(size_t threadId)
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
void LowHashNew::scanBucketsThreadFunction(size_t threadId)
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
            const MemoryAsContainer<BucketEntry> bucket = buckets[bucketId];
            if(bucket.size() < 2) {
                continue;
            }
            if(bucket.size() > maxBucketSize) {
                continue;
            }

            // Loop over pairs of bucket entries.
            for(const BucketEntry& feature0: bucket) {
                const OrientedReadId orientedReadId0 = feature0.orientedReadId;
                const uint32_t ordinal0 = feature0.ordinal;
                const auto kmerIds0 = kmerIds[orientedReadId0.getValue()].begin() + ordinal0;

                for(const BucketEntry& feature1: bucket) {
                    const OrientedReadId orientedReadId1 = feature1.orientedReadId;

                    // Only consider the ones where readId0 < readId1.
                    if(orientedReadId0.getReadId() >= orientedReadId1.getReadId()) {
                        continue;
                    }

                    const uint32_t ordinal1 = feature1.ordinal;
                    const auto kmerIds1 = kmerIds[orientedReadId1.getValue()].begin() + ordinal1;

                    // If the k-mers are not the same, this is a collision. Discard.
                    if(not std::equal(kmerIds0, kmerIds0+mLocal, kmerIds1)) {
                        continue;
                    }

                    // We found a common feature. Store it.
                    commonFeatures.push_back(CommonFeature(
                        orientedReadId0,
                        orientedReadId1,
                        ordinal0,
                        ordinal1));
                }
            }
        }
    }
}


// Add up the number of common feature found by all threads.
uint64_t LowHashNew::countTotalThreadCommonFeatures() const
{
    uint64_t n = 0;
    for(const auto& v: threadCommonFeatures) {
        n += v->size();
    }
    return n;
}



void LowHashNew::gatherCommonFeatures()
{
    commonFeatures.createNew(
            largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "tmp-CommonFeatures"),
            largeDataPageSize);
    commonFeatures.beginPass1(kmerIds.size());
    runThreads(&LowHashNew::gatherCommonFeaturesPass1, threadCount);
    commonFeatures.beginPass2();
    runThreads(&LowHashNew::gatherCommonFeaturesPass2, threadCount);
    commonFeatures.endPass2();
}
void LowHashNew::gatherCommonFeaturesPass1(size_t threadId)
{
    const MemoryMapped::Vector<CommonFeature>& v = *threadCommonFeatures[threadId];
    for(const CommonFeature& commonFeature: v) {
        commonFeatures.incrementCountMultithreaded(commonFeature.orientedReadIds[0].getValue());
    }
}
void LowHashNew::gatherCommonFeaturesPass2(size_t threadId)
{
    const MemoryMapped::Vector<CommonFeature>& v = *threadCommonFeatures[threadId];
    for(const CommonFeature& commonFeature: v) {
        commonFeatures.storeMultithreaded(
            commonFeature.orientedReadIds[0].getValue(),
            CommonFeatureInfo(commonFeature));
    }
}



// Process the common features.
// For each orientedReadId0, we look at all the CommonFeatureInfo we have
// and sort them by orientedReadId1, then by ordinals, and remove duplicates.
// We then find groups of at least minFrequency common features involving the
// same pair(orientedReadId0, orientedReadId1)
void LowHashNew::processCommonFeatures()
{
    const uint64_t readCount = kmerIds.size() / 2;
    const uint64_t batchSize = 1000;
    setupLoadBalancing(readCount, batchSize);
    runThreads(&LowHashNew::processCommonFeaturesThreadFunction, threadCount);
}
void LowHashNew::processCommonFeaturesThreadFunction(size_t threadId)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over ReadId's in this batch.
        for(ReadId readId0=ReadId(begin); readId0!=ReadId(end); readId0++) {
            for(Strand strand0=0; strand0<2; strand0++) {
                const OrientedReadId orientedReadId0(readId0, strand0);
                std::lock_guard<std::mutex> lock(mutex); // ************************** TAKE OUT!
                cout << " Working on orientedReadId0 " << orientedReadId0 << endl;
                const MemoryAsContainer<CommonFeatureInfo> features = commonFeatures[orientedReadId0.getValue()];

                // Deduplicate.
                const auto uniqueBegin = features.begin();
                auto uniqueEnd = features.end();
                sort(uniqueBegin, uniqueEnd);
                uniqueEnd = unique(uniqueBegin, uniqueEnd);

                /*
                for(auto it=uniqueBegin; it!=uniqueEnd; ++it) {
                    const CommonFeatureInfo& feature = *it;
                    cout <<
                        feature.orientedReadId1 << " " <<
                        feature.ordinals[0] << " " <<
                        feature.ordinals[1] << " " <<
                        int32_t(feature.ordinals[1]) - int32_t(feature.ordinals[0]) << "\n";
                }
                */

                // Loop over streaks of features with the same orientedReadId1.
                for(auto it=uniqueBegin; it!=uniqueEnd;) {
                    auto streakBegin = it;
                    auto streakEnd = streakBegin;
                    const OrientedReadId orientedReadId1 = streakBegin->orientedReadId1;
                    while(streakEnd!=uniqueEnd and streakEnd->orientedReadId1==orientedReadId1) {
                        ++streakEnd;
                    }

                    // If too few, skip.
                    if(streakEnd - streakBegin < int64_t(minFrequency)) {
                        it = streakEnd;
                        continue;
                    }

                    cout << "Common features of oriented reads " <<
                        orientedReadId0 << " " <<
                        orientedReadId1 << ":\n";
                    for(auto it=streakBegin; it!=streakEnd; ++it) {
                        const CommonFeatureInfo& feature = *it;
                        cout <<
                            feature.ordinals[0] << " " <<
                            feature.ordinals[1] << " " <<
                            int32_t(feature.ordinals[1]) - int32_t(feature.ordinals[0]) << "\n";
                    }
                    cout << "Marker count " <<
                        kmerIds[orientedReadId0.getValue()].size() << " " <<
                        kmerIds[orientedReadId1.getValue()].size() << ":\n";
                    it = streakEnd;
                }
            }
        }
    }

}
