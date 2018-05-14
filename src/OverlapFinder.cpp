// Nanopore2.
#include "OverlapFinder.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Standard library.
#include <chrono>



// Class OverlapFinder uses the MinHash algorithm to find pairs
// of overlapping oriented reads. It uses as features
// sequences of m consecutive markers.



OverlapFinder::OverlapFinder(
    size_t m,                       // Number of consecutive markers that define a feature.
    size_t minHashIterationCount,   // Number of minHash iterations.
    size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
    size_t maxBucketSize,           // The maximum size for a bucket to be used.
    size_t minFrequency,            // Minimum number of minHash hits for a pair to be considered an overlap.
    size_t threadCount,
    const MemoryMapped::Vector<KmerInfo>& kmerTable,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& compressedMarkers,
    MemoryMapped::Vector<Overlap>& overlaps,
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t>& overlapTable,
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize
    ) :
    MultithreadedObject(*this),
    m(m),
    maxBucketSize(maxBucketSize),
    minFrequency(minFrequency),
    kmerTable(kmerTable),
    compressedMarkers(compressedMarkers),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize)

{
    cout << timestamp << "MinHash begins." << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Create vectors containing only the k-mer ids of all markers.
    // This is used to speed up the computation of hash functions.
    cout << timestamp << "Creating markers for oriented reads." << endl;
    createMarkers();

    // The number of oriented reads, each with its own vector of markers.
    const ReadId readCount = ReadId(compressedMarkers.size());
    const size_t n = markers.size();
    cout << "There are " << readCount << " reads, " << n << " oriented reads." << endl;
    CZI_ASSERT(n == 2*readCount);
    orientedReadBucket.resize(n);
    overlapCandidates.resize(n);

    // Compute the number of buckets and the corresponding mask
    // used to convert a hash value to a bucket index.
    if(log2MinHashBucketCount > 31) {
        throw runtime_error("log2MinHashBucketCount can be up to 31.");
    }
    const uint32_t bucketCount = 1 << log2MinHashBucketCount;
    mask = bucketCount - 1;

    // Make space to count the number of overlaps found so far.
    totalOverlapCountByThread.resize(threadCount);

    // MinHash iteration loop.
    for(iteration=0; iteration<minHashIterationCount; iteration++) {
        cout << timestamp << "MinHash iteration " << iteration << " begins." << endl;
        const auto t0 = std::chrono::steady_clock::now();

        // Compute the bucket that each read belongs to.
        const size_t batchSize = 10000;
        setupLoadBalancing(n, batchSize);
        runThreads(&OverlapFinder::computeBuckets, threadCount);
        const auto t1 = std::chrono::steady_clock::now();

        // Construct the buckets.
        cout << timestamp << "Initializing the buckets." << endl;
        buckets.clear();
        buckets.resize(bucketCount);
        /*
        // This does not reallocate. Faster, But total bucket capacity increases a lot.
        for(auto& bucket: buckets) {
            bucket.clear();
        }
        */
        cout << timestamp << "Filling in the buckets." << endl;
        for(ReadId i=0; i<n; i++) {
            buckets[orientedReadBucket[i]].push_back(i);
        }

        // Inspect the buckets to find overlap candidates.
        cout << timestamp << "Inspecting the buckets." << endl;
        const auto t2 = std::chrono::steady_clock::now();
        setupLoadBalancing(n, batchSize);
        runThreads(&OverlapFinder::inspectBuckets, threadCount, "threadLogs/inspectBuckets");
        const auto t3 = std::chrono::steady_clock::now();

        const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
        const double t12 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count());
        const double t23 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2)).count());
        cout << "Times for this iteration: " << t01 << " " << t12 << " " << t23 << endl;
        const size_t totalOverlapCount =
            accumulate(totalOverlapCountByThread.begin(), totalOverlapCountByThread.end(), 0);
        cout << "Found " << totalOverlapCount;
        cout << " overlaps with frequency at least " << minFrequency << " so far." << endl;

        // Write out total capacity for buckets and overlap candidates.
        size_t totalBucketCapacity = 0;
        for(const auto& bucket: buckets) {
            totalBucketCapacity += bucket.capacity();
        }
        size_t totalCandidateCapacity = 0;
        size_t totalCandidateSize = 0;
        for(const auto& v: overlapCandidates) {
            totalCandidateCapacity += v.capacity();
            totalCandidateSize += v.size();
        }
        cout << "Total bucket capacity: " << totalBucketCapacity << endl;
        cout << "Overlap candidates: total size " << totalCandidateSize;
        cout << ", total capacity  " << totalCandidateCapacity << endl;

    }



    // Create the overlaps.
    cout << timestamp << "Storing overlaps." << endl;
    CZI_ASSERT(n == 2*readCount);
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        for(Strand strand0=0; strand0<2; strand0++) {
            const OrientedReadId orientedReadId0(readId0, strand0);
            const auto& candidates = overlapCandidates[orientedReadId0.getValue()];
            for(const auto& candidate: candidates) {
                if(candidate.second >= minFrequency) {
                    const OrientedReadId orientedReadId1(candidate.first);
                    CZI_ASSERT(orientedReadId1.getValue() == candidate.first);
                    CZI_ASSERT(orientedReadId1 < orientedReadId0);
                    overlaps.push_back(Overlap(orientedReadId1, orientedReadId0, candidate.second));
                }
            }

        }
    }
    cout << "Found " << overlaps.size() << " overlaps."<< endl;
    cout << "Average number of overlaps per oriented read is ";
    cout << (2.* double(overlaps.size())) / double(n)  << endl;


    // Create the overlap table.
    cout << timestamp << "Creating overlap table." << endl;
    overlapTable.beginPass1(n);
    for(const Overlap& overlap: overlaps) {
        CZI_ASSERT(overlap.orientedReadIds[0].getValue() < n);
        CZI_ASSERT(overlap.orientedReadIds[1].getValue() < n);
        overlapTable.incrementCount(overlap.orientedReadIds[0].getValue());
        overlapTable.incrementCount(overlap.orientedReadIds[1].getValue());
    }
    overlapTable.beginPass2();
    for(size_t i=overlaps.size()-1; ; i--) {
        const Overlap& overlap = overlaps[i];
        overlapTable.store(overlap.orientedReadIds[0].getValue(), i);
        overlapTable.store(overlap.orientedReadIds[1].getValue(), i);
        if(i==0) {
            break;
        }
    }
    overlapTable.endPass2();
    CZI_ASSERT(overlapTable.totalSize() == 2*overlaps.size());



    // Cleanup.
    markers.remove();
    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    cout << timestamp << "MinHash completed in " << tTotal << " s." << endl;

}



void OverlapFinder::createMarkers()
{
    markers.createNew(
        largeDataFileNamePrefix + "tmp-OverlapFinder-Markers",
        largeDataPageSize);

    vector<KmerId> kmerIds;

    const ReadId readCount = ReadId(compressedMarkers.size());
    for(ReadId readId=0; readId!=readCount; readId++) {
        const auto readCompressedMarkers = compressedMarkers[readId];
        kmerIds.clear();

        // Add the markers of the read without reverse complementing.
        markers.appendVector();
        for(const CompressedMarker& compressedMarker: readCompressedMarkers) {
            const KmerId kmerId = compressedMarker.kmerId;
            markers.append(kmerId);
            kmerIds.push_back(kmerId);
        }

        // Add the markers of the read with reverse complementing.
        markers.appendVector();
        for(auto it=kmerIds.rbegin(); it!=kmerIds.rend(); ++it) {
            markers.append(kmerTable[*it].reverseComplementedKmerId);
        }
    }
    CZI_ASSERT(markers.size() == 2*compressedMarkers.size());
}



// Thread function used to compute the bucket that each
// oriented read belongs to.
void OverlapFinder::computeBuckets(size_t threadId)
{
    const int featureByteCount = int(m * sizeof(KmerId));
    const uint32_t seed = uint32_t(iteration*37);

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads assigned to this batch.
        for(size_t i=begin; i!=end; i++) {
            const size_t markerCount = markers.size(i);
            uint32_t minHash = std::numeric_limits<uint32_t>::max();

            // If there are not enough markers, the minHash
            // stays to this initial value and the oriented read
            // ends up in the final bucket.
            // We will ignore all reads in the final bucket.
            if(markerCount >= m) {

                // Get the markers for this oriented read.
                KmerId* kmerIds = markers.begin(i);
                const size_t featureCount = markers.size(i) - m + 1;

                // Loop over features of this oriented read.
                // Features are sequences of m consecutive markers.
                for(size_t j=0; j<featureCount; j++, kmerIds++) {
                    const uint32_t hash = MurmurHash2(kmerIds, featureByteCount, seed);
                    minHash = min(hash, minHash);
                }
            }
            orientedReadBucket[i] = minHash & mask;
        }
    }

}



// Thread function used to inspect the buckets to find overlap candidates.
void OverlapFinder::inspectBuckets(size_t threadId)
{
    ostream& out = getLog(threadId);
    size_t& totalCountForThread = totalOverlapCountByThread[threadId];
    totalCountForThread = 0;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << "Working on block " << begin << " " << end << endl;

        // Loop over oriented reads assigned to this batch.
        for(size_t i=begin; i!=end; i++) {

            // Locate the candidate overlaps for this oriented read.
            auto& candidates = overlapCandidates[i];

            const uint32_t bucketId = orientedReadBucket[i];

            // Ignore the last bucket.
            // It contains (mostly) reads with less than m markers.
            if(bucketId != mask) {

                // Locate the bucket that contains this oriented read.
                const auto& bucket = buckets[bucketId];


                // Loop over bucket entries, but do nothing if the bucket is too big.
                if(bucket.size() <= maxBucketSize) {
                    for(const auto j: bucket) {
                        if(j >= i) {
                            break;
                        }

                        bool found = false;
                        for(auto& candidate: candidates) {
                            if(candidate.first == j) {
                                // We found it. Increment the frequency.
                                ++candidate.second;
                                found = true;
                                break;
                            }
                        }
                        if(!found) {
                            candidates.push_back(make_pair(j, 1));
                        }
                    }
                }
            }

            // Increment the total number of overlaps
            // (used for statistics only).
            for(const auto& candidate: candidates) {
                if(candidate.second >= minFrequency) {
                    ++totalCountForThread;
                }
            }
        }
    }

}

