// Shasta.
#include "MinHash.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Standard library.
#include "chrono.hpp"
#include <numeric>



// Class MinHash uses the MinHash algorithm to find candidate pairs
// of aligned reads. It uses as features
// sequences of m consecutive markers.



MinHash::MinHash(
    size_t m,                       // Number of consecutive markers that define a feature.
    size_t minHashIterationCount,   // Number of minHash iterations.
    size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
    size_t maxBucketSize,           // The maximum size for a bucket to be used.
    size_t minFrequency,            // Minimum number of minHash hits for a pair to be considered a candidate.
    size_t threadCountArgument,
    const MemoryMapped::Vector<KmerInfo>& kmerTable,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    MemoryMapped::Vector<OrientedReadPair>& candidateAlignments,
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize
    ) :
    MultithreadedObject(*this),
    m(m),
    maxBucketSize(maxBucketSize),
    minFrequency(minFrequency),
    threadCount(threadCountArgument),
    kmerTable(kmerTable),
    markers(markers),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize)

{
    cout << timestamp << "MinHash begins." << endl;
    const auto tBegin = steady_clock::now();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Create vectors containing only the k-mer ids of all markers.
    // This is used to speed up the computation of hash functions.
    cout << timestamp << "Creating kmer ids for oriented reads." << endl;
    createKmerIds();

    // The number of oriented reads, each with its own vector of markers.
    const OrientedReadId::Int orientedReadCount = OrientedReadId::Int(markers.size());
    const ReadId readCount = orientedReadCount / 2;
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
            largeDataFileNamePrefix + "tmp-MinHash-Buckets",
            largeDataPageSize);
    minHash.resize(orientedReadCount);
    candidates.resize(readCount);
    totalCandidateCountByThread.resize(threadCount);

    // MinHash iteration loop.
    for(iteration=0; iteration<minHashIterationCount; iteration++) {
        cout << timestamp << "MinHash iteration " << iteration << " begins." << endl;
        const auto t0 = steady_clock::now();

        // Compute the min hash of each oriented read.
        const size_t batchSize = 10000;
        setupLoadBalancing(orientedReadCount, batchSize);
        runThreads(&MinHash::computeMinHash, threadCount);
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

        // Inspect the buckets to find candidates.
        const auto t2 = steady_clock::now();
        setupLoadBalancing(readCount, batchSize);
        runThreads(&MinHash::inspectBuckets, threadCount);
        const auto t3 = steady_clock::now();

        // Write a summary fo rthis iteration.
        const double t01 = seconds(t1 - t0);
        const double t12 = seconds(t2 - t1);
        const double t23 = seconds(t3 - t2);
        cout << "Times for this iteration: hash " << t01 << ", fill " << t12 << ", inspect " << t23 << endl;
        const size_t totalCandidateCount =
            std::accumulate(totalCandidateCountByThread.begin(), totalCandidateCountByThread.end(), 0ULL);
        cout << "Found " << totalCandidateCount;
        cout << " candidates with frequency at least " << minFrequency << " so far." << endl;

        // Write out total capacity for candidates.
        size_t totalCandidateCapacity = 0;
        size_t totalCandidateSize = 0;
        for(const auto& v: candidates) {
            totalCandidateCapacity += v.capacity();
            totalCandidateSize += v.size();
        }
        cout << "Candidates: total size " << totalCandidateSize;
        cout << ", total capacity  " << totalCandidateCapacity << endl;

    }



    // Remove the buckets and the markers. They are no longer needed.
    buckets.remove();
    kmerIds.remove();



    // Create the candidate alignments.
    cout << timestamp << "Storing candidate alignments." << endl;
    CZI_ASSERT(orientedReadCount == 2*readCount);
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        const auto& candidates0 = candidates[readId0];
        for(const Candidate& candidate: candidates0) {
            if(candidate.frequency >= minFrequency) {
                const ReadId readId1 = candidate.readId1;
                CZI_ASSERT(readId0 < readId1);
                candidateAlignments.push_back(
                    OrientedReadPair(readId0, readId1, candidate.isSameStrand));
            }
        }
    }
    cout << "Found " << candidateAlignments.size() << " candidates."<< endl;
    cout << "Average number of candidates per oriented read is ";
    cout << (2.* double(candidateAlignments.size())) / double(orientedReadCount)  << endl;



    // Done.
    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "MinHash completed in " << tTotal << " s." << endl;

}



void MinHash::createKmerIds()
{
    kmerIds.createNew(
        largeDataFileNamePrefix + "tmp-MinHash-Markers",
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
    runThreads(&MinHash::createKmerIds, threadCount);
}



// Thread function for createKmerIds.
void MinHash::createKmerIds(size_t threadId)
{

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];

                CZI_ASSERT(kmerIds.size(orientedReadId.getValue()) == orientedReadMarkers.size());

                auto pointer = kmerIds.begin(orientedReadId.getValue());
                for(const CompressedMarker& marker: orientedReadMarkers) {
                    *pointer++ = marker.kmerId;
                }
            }
        }
    }
}



// Thread function used to compute the min hash of each oriented read.
void MinHash::computeMinHash(size_t threadId)
{
    const int featureByteCount = int(m * sizeof(KmerId));
    const uint64_t seed = iteration * 37;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads assigned to this batch.
        for(size_t i=begin; i!=end; i++) {
            const size_t markerCount = kmerIds.size(i);
            uint64_t minHashValue = std::numeric_limits<uint64_t>::max();

            // If there are not enough markers, the minHash
            // stays to this initial value and the oriented read
            // ends up in the final bucket.
            // We will ignore all reads in the final bucket.
            if(markerCount >= m) {

                // Get the markers for this oriented read.
                KmerId* kmerIdsPointer = kmerIds.begin(i);
                const size_t featureCount = kmerIds.size(i) - m + 1;

                // Loop over features of this oriented read.
                // Features are sequences of m consecutive markers.
                for(size_t j=0; j<featureCount; j++, kmerIdsPointer++) {
                    const uint64_t hash = MurmurHash64A(kmerIdsPointer, featureByteCount, seed);
                    minHashValue = min(hash, minHashValue);
                }
            }
            minHash[i] = minHashValue;
        }
    }

}



// Thread function used to inspect the buckets to find candidates.
void MinHash::inspectBuckets(size_t threadId)
{
    size_t& totalCountForThread = totalCandidateCountByThread[threadId];
    totalCountForThread = 0;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads assigned to this batch.
        for(ReadId readId0=ReadId(begin); readId0!=ReadId(end); readId0++) {

            // Locate the candidates for this oriented read.
            auto& candidates0 = candidates[readId0];

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
                            // Don't generate candidates with self, on either strand.
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
                                for(Candidate& candidate: candidates0) {
                                    if(candidate.readId1==readId1 && candidate.isSameStrand==isSameStrand) {
                                        // We found it. Increment the frequency.
                                        ++candidate.frequency;
                                        found = true;
                                        break;
                                    }
                                }
                                if(!found) {
                                    candidates0.push_back(Candidate(readId1, 1, isSameStrand));
                                }
                            }
                        }
                    }
                }
            }

            // Increment the total number of candidates
            // (used for statistics only).
            for(const auto& candidate: candidates0) {
                if(candidate.frequency >= minFrequency) {
                    ++totalCountForThread;
                }
            }
        }
    }

}

