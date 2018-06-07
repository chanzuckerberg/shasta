// Nanopore2.
#include "OverlapFinder.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Standard library.
#include "chrono.hpp"



// Class OverlapFinder uses the MinHash algorithm to find pairs
// of overlapping oriented reads. It uses as features
// sequences of m consecutive markers.



OverlapFinder::OverlapFinder(
    size_t m,                       // Number of consecutive markers that define a feature.
    size_t minHashIterationCount,   // Number of minHash iterations.
    size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
    size_t maxBucketSize,           // The maximum size for a bucket to be used.
    size_t minFrequency,            // Minimum number of minHash hits for a pair to be considered an overlap.
    size_t threadCountArgument,
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
    threadCount(threadCountArgument),
    kmerTable(kmerTable),
    compressedMarkers(compressedMarkers),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize)

{
    cout << timestamp << "MinHash begins." << endl;
    const auto tBegin = steady_clock::now();

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
    const OrientedReadId::Int orientedReadCount = OrientedReadId::Int(markers.size());
    cout << "There are " << readCount << " reads, " << orientedReadCount << " oriented reads." << endl;
    CZI_ASSERT(orientedReadCount == 2*readCount);

    // Compute the number of buckets and the corresponding mask
    // used to convert a hash value to a bucket index.
    if(log2MinHashBucketCount > 31) {
        throw runtime_error("log2MinHashBucketCount can be up to 31.");
    }
    const uint32_t bucketCount = 1 << log2MinHashBucketCount;
    mask = bucketCount - 1;

    // Prepare work areas.
    buckets.createNew(
            largeDataFileNamePrefix + "tmp-OverlapFinder-Buckets",
            largeDataPageSize);
    minHash.resize(orientedReadCount);
    overlapCandidates.resize(readCount);
    totalOverlapCountByThread.resize(threadCount);

    // MinHash iteration loop.
    for(iteration=0; iteration<minHashIterationCount; iteration++) {
        cout << timestamp << "MinHash iteration " << iteration << " begins." << endl;
        const auto t0 = steady_clock::now();

        // Compute the min hash of each oriented read.
        const size_t batchSize = 10000;
        setupLoadBalancing(orientedReadCount, batchSize);
        runThreads(&OverlapFinder::computeMinHash, threadCount);
        const auto t1 = steady_clock::now();

        // Construct the buckets.
        const auto tb0 = steady_clock::now();
        buckets.clear();
        const auto tb1 = steady_clock::now();
        buckets.beginPass1(bucketCount);
        const auto tb2 = steady_clock::now();
        for(OrientedReadId::Int orientedReadId=0; orientedReadId<orientedReadCount; orientedReadId++) {
            const uint64_t bucketId = minHash[orientedReadId] & mask;
            buckets.incrementCount(bucketId);
        }
        const auto tb3 = steady_clock::now();
        buckets.beginPass2();
        const auto tb4 = steady_clock::now();
        for(OrientedReadId::Int orientedReadId=0; orientedReadId<orientedReadCount; orientedReadId++) {
            const uint64_t bucketId = minHash[orientedReadId] & mask;
            buckets.store(bucketId, orientedReadId);
        }
        const auto tb5 = steady_clock::now();
        buckets.endPass2(false, false);
        const auto tb6 = steady_clock::now();
        cout << "Bucket filling times: ";
        cout << seconds(tb1-tb0) << " ";
        cout << seconds(tb2-tb1) << " ";
        cout << seconds(tb3-tb2) << " ";
        cout << seconds(tb4-tb3) << " ";
        cout << seconds(tb5-tb4) << " ";
        cout << seconds(tb6-tb5) << " ";
        cout << endl;

        // Inspect the buckets to find overlap candidates.
        const auto t2 = steady_clock::now();
        setupLoadBalancing(readCount, batchSize);
        runThreads(&OverlapFinder::inspectBuckets, threadCount, "threadLogs/inspectBuckets");
        const auto t3 = steady_clock::now();

        const double t01 = seconds(t1 - t0);
        const double t12 = seconds(t2 - t1);
        const double t23 = seconds(t3 - t2);
        cout << "Times for this iteration: hash " << t01 << ", fill " << t12 << ", inspect " << t23 << endl;
        const size_t totalOverlapCount =
            accumulate(totalOverlapCountByThread.begin(), totalOverlapCountByThread.end(), 0ULL);
        cout << "Found " << totalOverlapCount;
        cout << " overlaps with frequency at least " << minFrequency << " so far." << endl;

        // Write out total capacity for overlap candidates.
        size_t totalCandidateCapacity = 0;
        size_t totalCandidateSize = 0;
        for(const auto& v: overlapCandidates) {
            totalCandidateCapacity += v.capacity();
            totalCandidateSize += v.size();
        }
        cout << "Overlap candidates: total size " << totalCandidateSize;
        cout << ", total capacity  " << totalCandidateCapacity << endl;

    }



    // Remove the buckets and the markers. They are no longer needed.
    buckets.remove();
    markers.remove();



    // Create the overlaps.
    cout << timestamp << "Storing overlaps." << endl;
    CZI_ASSERT(orientedReadCount == 2*readCount);
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        const auto& candidates = overlapCandidates[readId0];
        for(const OverlapCandidate& candidate: candidates) {
            if(candidate.frequency >= minFrequency) {
                const ReadId readId1 = candidate.readId1;
                CZI_ASSERT(readId0 < readId1);
                    overlaps.push_back(Overlap(readId0, readId1, candidate.isSameStrand, candidate.frequency));
            }
        }
    }
    cout << "Found " << overlaps.size() << " overlaps."<< endl;
    cout << "Average number of overlaps per oriented read is ";
    cout << (2.* double(overlaps.size())) / double(orientedReadCount)  << endl;



    // Create the overlap table.
    // It contains indexes of the overlaps that each OrientedReadId
    // is involved in.
    cout << timestamp << "Creating overlap table." << endl;
    overlapTable.beginPass1(orientedReadCount);
    for(const Overlap& overlap: overlaps) {
        OrientedReadId orientedReadId0(overlap.readIds[0], 0);
        OrientedReadId orientedReadId1(overlap.readIds[1], overlap.isSameStrand ? 0 : 1);
        overlapTable.incrementCount(orientedReadId0.getValue());
        overlapTable.incrementCount(orientedReadId1.getValue());
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        overlapTable.incrementCount(orientedReadId0.getValue());
        overlapTable.incrementCount(orientedReadId1.getValue());
    }
    overlapTable.beginPass2();
    if(overlaps.size() > 0) {
        for(size_t i=overlaps.size()-1; ; i--) {
            const Overlap& overlap = overlaps[i];
            OrientedReadId orientedReadId0(overlap.readIds[0], 0);
            OrientedReadId orientedReadId1(overlap.readIds[1], overlap.isSameStrand ? 0 : 1);
            overlapTable.store(orientedReadId0.getValue(), i);
            overlapTable.store(orientedReadId1.getValue(), i);
            orientedReadId0.flipStrand();
            orientedReadId1.flipStrand();
            overlapTable.store(orientedReadId0.getValue(), i);
            overlapTable.store(orientedReadId1.getValue(), i);
            if(i==0) {
                break;
            }
        }
    }
    overlapTable.endPass2();
    CZI_ASSERT(overlapTable.totalSize() == 4*overlaps.size());



    // Done.
    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "MinHash completed in " << tTotal << " s." << endl;

}


#if 0
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
#else



void OverlapFinder::createMarkers()
{
    markers.createNew(
        largeDataFileNamePrefix + "tmp-OverlapFinder-Markers",
        largeDataPageSize);
    const ReadId readCount = ReadId(compressedMarkers.size());
    markers.beginPass1(2*readCount);
    for(ReadId readId=0; readId!=readCount; readId++) {
        const auto markerCount = compressedMarkers.size(readId);
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            markers.incrementCount(orientedReadId.getValue(), markerCount);
        }
    }
    markers.beginPass2();
    markers.endPass2(false);
    const size_t batchSize = 10000;
    setupLoadBalancing(readCount, batchSize);
    runThreads(&OverlapFinder::createMarkers, threadCount, "threadLogs/createMarkers");
}



// Thread function for createMarkers.
void OverlapFinder::createMarkers(size_t threadId)
{
    ostream& out = getLog(threadId);

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        out << begin << " " << end << endl;

        // Loop over reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {
            const auto readCompressedMarkers = compressedMarkers[readId];
            const OrientedReadId orientedReadIdStrand0(readId, 0);
            const OrientedReadId orientedReadIdStrand1(readId, 1);

            CZI_ASSERT(markers.size(orientedReadIdStrand0.getValue()) == readCompressedMarkers.size());
            CZI_ASSERT(markers.size(orientedReadIdStrand1.getValue()) == readCompressedMarkers.size());

            auto pointer = markers.begin(orientedReadIdStrand0.getValue());
            for(const CompressedMarker& compressedMarker: readCompressedMarkers) {
                *pointer++ = compressedMarker.kmerId;
            }

            pointer = markers.begin(orientedReadIdStrand1.getValue());
            for(auto it=readCompressedMarkers.end()-1; it>=readCompressedMarkers.begin(); --it) {
                *pointer++ = kmerTable[it->kmerId].reverseComplementedKmerId;
            }

        }
    }

}
#endif


// Thread function used to compute the min hash of each oriented read.
void OverlapFinder::computeMinHash(size_t threadId)
{
    const int featureByteCount = int(m * sizeof(KmerId));
    const uint64_t seed = iteration * 37;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads assigned to this batch.
        for(size_t i=begin; i!=end; i++) {
            const size_t markerCount = markers.size(i);
            uint64_t minHashValue = std::numeric_limits<uint64_t>::max();

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
                    const uint64_t hash = MurmurHash64A(kmerIds, featureByteCount, seed);
                    minHashValue = min(hash, minHashValue);
                }
            }
            minHash[i] = minHashValue;
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

        // Loop over reads assigned to this batch.
        for(ReadId readId0=ReadId(begin); readId0!=ReadId(end); readId0++) {

            // Locate the candidate overlaps for this oriented read.
            auto& candidates = overlapCandidates[readId0];

            // Loop over two strands.
            for(Strand strand0=0; strand0<2; strand0++) {
                const OrientedReadId orientedReadId0(readId0, strand0);
                const OrientedReadId::Int orientedReadIdInt0 = orientedReadId0.getValue();
                const uint64_t minHash0 = minHash[orientedReadIdInt0];

                const uint64_t bucketId = minHash0 & mask;

                // Ignore the last bucket.
                // It contains (mostly) reads with less than m markers.
                if(bucketId != mask) {

                    // Locate the bucket that contains this oriented read.
                    const auto& bucket = buckets[bucketId];


                    // Loop over bucket entries, but do nothing if the bucket is too big.
                    if(bucket.size() <= maxBucketSize) {
                        for(const OrientedReadId::Int orientedReadIdInt1: bucket) {
                            const OrientedReadId orientedReadId1 = OrientedReadId(orientedReadIdInt1);
                            const ReadId readId1 = orientedReadId1.getReadId();
                            const Strand strand1 = orientedReadId1.getStrand();
                            const bool isSameStrand = (strand1 == strand0);

                            // Only consider pairs with readId0 < readId1.
                            // Ignore overlaps with self, on either strand.
                            if(readId1 <= readId0) {
                                continue;
                            }

                            // To avoid most collision, only consider it
                            // if not the min hash values are the same.
                            if(minHash[orientedReadIdInt1] == minHash[orientedReadIdInt0]) {

                                // Look for this candidate.
                                // If found, increase its frequency.
                                // Otherwise, add it with frequency 1.
                                bool found = false;
                                for(OverlapCandidate& candidate: candidates) {
                                    if(candidate.readId1==readId1 && candidate.isSameStrand==isSameStrand) {
                                        // We found it. Increment the frequency.
                                        ++candidate.frequency;
                                        found = true;
                                        break;
                                    }
                                }
                                if(!found) {
                                    candidates.push_back(OverlapCandidate(readId1, 1, isSameStrand));
                                }
                            }
                        }
                    }
                }
            }

            // Increment the total number of overlaps
            // (used for statistics only).
            for(const auto& candidate: candidates) {
                if(candidate.frequency >= minFrequency) {
                    ++totalCountForThread;
                }
            }
        }
    }

}

