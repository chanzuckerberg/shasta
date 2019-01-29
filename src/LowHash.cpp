// Shasta.
#include "LowHash.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Standard library.
#include "chrono.hpp"
#include <numeric>



// Class LowHash uses the LowHash algorithm to find candidate pairs
// of aligned reads. It uses as features
// sequences of m consecutive markers.



LowHash::LowHash(
    size_t m,                       // Number of consecutive markers that define a feature.
    double hashFraction,
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
    hashFraction(hashFraction),
    maxBucketSize(maxBucketSize),
    minFrequency(minFrequency),
    threadCount(threadCountArgument),
    kmerTable(kmerTable),
    markers(markers),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize)

{
    cout << timestamp << "LowHash begins." << endl;
    const auto tBegin = steady_clock::now();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Create vectors containing only the k-mer ids of all markers.
    // This is used to speed up the computation of hash functions.
    cout << timestamp << "Creating kmer ids for oriented reads." << endl;
    createKmerIds();

    // Compute the threshold for a hash value to be considered low.
    hashThreshold = uint64_t(double(hashFraction) * double(std::numeric_limits<uint64_t>::max()));

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
    log2MinHashBucketCount = 26;   // ************************** FOR NOW
    const uint32_t bucketCount = 1 << log2MinHashBucketCount;
    mask = bucketCount - 1;



    // Set up work areas.
    buckets.createNew(
            largeDataFileNamePrefix + "tmp-LowHash-Buckets",
            largeDataPageSize);
    lowHashes.resize(orientedReadCount);
    candidates.resize(readCount);
    threadStatistics.resize(threadCount);



    // LowHash iteration loop.
    for(iteration=0; iteration<minHashIterationCount; iteration++) {
        cout << timestamp << "LowHash iteration " << iteration << " begins." << endl;

        // Pass1: compute the low hashes for each oriented read
        // and prepare the buckets for filling.
        buckets.clear();
        buckets.beginPass1(bucketCount);
        size_t batchSize = 10000;
        setupLoadBalancing(orientedReadCount, batchSize);
        runThreads(&LowHash::pass1ThreadFunction, threadCount);

        // Pass 2: fill the buckets.
        buckets.beginPass2();
        batchSize = 10000;
        setupLoadBalancing(orientedReadCount, batchSize);
        runThreads(&LowHash::pass2ThreadFunction, threadCount);
        buckets.endPass2(false, false);

        // Pass 3: sort the buckets.
        batchSize = 1000;
        setupLoadBalancing(bucketCount, batchSize);
        runThreads(&LowHash::pass3ThreadFunction, threadCount);

        // Pass 4: inspect the buckets to find candidates.
        batchSize = 10000;
        setupLoadBalancing(readCount, batchSize);
        runThreads(&LowHash::pass4ThreadFunction, threadCount,
            "threadLogs/LowHash-pass4-iteration" + to_string(iteration));

        // Write a summary for this iteration.
        uint64_t highFrequency = 0;
        uint64_t total = 0;
        uint64_t capacity = 0;
        for(const auto& s: threadStatistics) {
            highFrequency += s.highFrequency;
            total += s.total;
            capacity += s.capacity;
        }
        cout << "Alignment candidates after iteration " << iteration;
        cout << ": high frequency " << highFrequency;
        cout << ", total " << total;
        cout << ", capacity " << capacity << "." << endl;
    }



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
    cout << "Found " << candidateAlignments.size() << " alignment candidates."<< endl;
    cout << "Average number of alignment candidates per oriented read is ";
    cout << (2.* double(candidateAlignments.size())) / double(orientedReadCount)  << "." << endl;



    // Clean up work areas.
    buckets.remove();
    kmerIds.remove();



    // Done.
    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "LowHash completed in " << tTotal << " s." << endl;
}



void LowHash::createKmerIds()
{
    kmerIds.createNew(
        largeDataFileNamePrefix + "tmp-LowHash-Markers",
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
    runThreads(&LowHash::createKmerIds, threadCount, "threadLogs/createKmerIds");
}



// Thread function for createKmerIds.
void LowHash::createKmerIds(size_t threadId)
{
    ostream& out = getLog(threadId);

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        out << begin << " " << end << endl;

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



// Pass1: compute the low hashes for each oriented read
// and prepare the buckets for filling.
void LowHash::pass1ThreadFunction(size_t threadId)
{
    const int featureByteCount = int(m * sizeof(KmerId));
    const uint64_t seed = iteration * 37;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads assigned to this batch.
        for(size_t i=begin; i!=end; i++) {
            vector<uint64_t>& orientedReadLowHashes = lowHashes[i];
            orientedReadLowHashes.clear();
            const size_t markerCount = kmerIds.size(i);

            // Handle the pathological case where there are fewer than m markers.
            // This oriented read ends up in no bucket.
            if(markerCount < m) {
                continue;
            }


            // Get the markers for this oriented read.
            KmerId* kmerIdsPointer = kmerIds.begin(i);
            const size_t featureCount = markerCount - m + 1;

            // Loop over features of this oriented read.
            // Features are sequences of m consecutive markers.
            for(size_t j=0; j<featureCount; j++, kmerIdsPointer++) {
                const uint64_t hash = MurmurHash64A(kmerIdsPointer, featureByteCount, seed);
                if(hash < hashThreshold) {
                    orientedReadLowHashes.push_back(hash);
                    const uint64_t bucketId = hash & mask;
                    buckets.incrementCountMultithreaded(bucketId);
                }
            }
        }
    }

}



// Pass 2: fill the buckets.
void LowHash::pass2ThreadFunction(size_t threadId)
{

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over oriented reads assigned to this batch.
        for(size_t i=begin; i!=end; i++) {
            const vector<uint64_t>& orientedReadLowHashes = lowHashes[i];
            for(const uint64_t hash: orientedReadLowHashes) {
                const uint64_t bucketId = hash & mask;
                buckets.storeMultithreaded(bucketId, OrientedReadId::Int(i));
            }
        }
    }
}



// Pass 3: sort the buckets.
void LowHash::pass3ThreadFunction(size_t threadId)
{
    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over buckets assigned to this batch.
        for(size_t bucketId=begin; bucketId!=end; bucketId++) {
            const MemoryAsContainer<OrientedReadId::Int> bucket = buckets[bucketId];
            sort(bucket.begin(), bucket.end());
        }
    }

}



// Pass 4: inspect the buckets to find candidates.
void LowHash::pass4ThreadFunction(size_t threadId)
{
    ostream& out = getLog(threadId);
    const bool debug = true;

    // Work area used to store the new candidates for a single readId0,
    // for both strands. The second member of the pair is isSameStrand.
    vector< pair<ReadId, bool> > workArea;

    // Work area to store the merged candidates for s single readId0.
    vector<Candidate> mergedCandidates0;

    ThreadStatistics& thisThreadStatistics = threadStatistics[threadId];
    thisThreadStatistics.clear();

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads assigned to this batch.
        for(ReadId readId0=ReadId(begin); readId0!=ReadId(end); readId0++) {
            workArea.clear();

            // Loop over two strands.
            for(Strand strand0=0; strand0<2; strand0++) {
                const OrientedReadId orientedReadId0(readId0, strand0);
                const OrientedReadId::Int orientedReadIdInt0 = orientedReadId0.getValue();

                // Loop over the low hashes for this oriented read.
                const vector<uint64_t>& orientedReadLowHashes = lowHashes[orientedReadIdInt0];
                for(const uint64_t hash: orientedReadLowHashes) {

                    // Loop over oriented read ids in the bucket corresponding to this hash.
                    const uint64_t bucketId = hash & mask;
                    const MemoryAsContainer<OrientedReadId::Int> bucket = buckets[bucketId];
                    for(OrientedReadId::Int orientedReadIdInt1: bucket) {
                        const OrientedReadId orientedReadId1 = OrientedReadId(orientedReadIdInt1);
                        const ReadId readId1 = orientedReadId1.getReadId();

                        // Only consider it if readId1 > readId0.
                        if(readId1 <= readId0) {
                            continue;
                        }

                        // Add it to our work area.
                        const bool isSameStrand = orientedReadId1.getStrand() == strand0;
                        workArea.push_back(make_pair(readId1, isSameStrand));
                    }
                }
            }

            // Sort  and deduplicate the work area.
            sort(workArea.begin(), workArea.end());
            workArea.resize(unique(workArea.begin(), workArea.end()) - workArea.begin());






            // Merge the contents of the work area
            // with the candidates previously stored.
            vector<Candidate>& candidates0 = candidates[readId0];
            if(debug) {
                out << "Read id " << readId0 << endl;
                out << "workArea:" << endl;
                for(const auto& p: workArea) {
                    out << p.first << " " << p.second << endl;
                }
                out << "candidates0:" << endl;
                for(const Candidate& candidate: candidates0) {
                    out << candidate.readId1 << " " << candidate.isSameStrand << " " << candidate.frequency << endl;
                }
            }
            mergedCandidates0.clear();
            vector<Candidate>::iterator begin0 = candidates0.begin();
            vector<Candidate>::iterator end0 = candidates0.end();
            vector< pair<ReadId, bool> >::iterator beginW = workArea.begin();
            vector< pair<ReadId, bool> >::iterator endW = workArea.end();
            vector<Candidate>::iterator it0 = begin0;
            vector< pair<ReadId, bool> >::iterator itW = beginW;
            while(true) {

                // If we reached the end of the work area, just copy the
                // remaining entries in the candidate list.
                if(itW == endW) {
                    copy(it0, end0, back_inserter(mergedCandidates0));
                    break;
                }

                // If we reached the end of the existing candidate list,
                // copy the remaining work area entries, with frequency 1.
                if(it0 == end0) {
                    for(; itW != endW; ++itW) {
                        mergedCandidates0.push_back(Candidate(itW->first, 1, itW->second));
                    }
                    break;
                }

                if(*itW < make_pair(it0->readId1, it0->isSameStrand)) {
                    // The work area has an entry that is not in the existing
                    // candidate list for this read.
                    // Create a new entry with frequency 1.
                    mergedCandidates0.push_back(Candidate(itW->first, 1, itW->second));
                    ++itW;
                }

                else if(make_pair(it0->readId1, it0->isSameStrand) < *itW) {
                    // The existing candidate list for this read has an entry
                    // that is not in the work area.
                    // Just copy it.
                    mergedCandidates0.push_back(*it0);
                    ++it0;
                }

                else {
                    // The work area has an entry that is already in the existing
                    // candidate list for this read.
                    // Merged then by incrementing its frequency.
                    mergedCandidates0.push_back(*it0);
                    mergedCandidates0.back().frequency++;
                    ++itW;
                    ++it0;
                }

            }
            // out << readId0 << " " << candidates0.size() << " " << workArea.size();
            // out << " " << mergedCandidates0.size() << endl;

            if(debug) {
                out << "mergedCandidates0:" << endl;
                for(const Candidate& candidate: mergedCandidates0) {
                    out << candidate.readId1 << " " << candidate.isSameStrand << " " << candidate.frequency << endl;
                }
            }


            // Store the merged candidates in place of the old ones.
            candidates0.swap(mergedCandidates0);

            // Update thread statistics.
            thisThreadStatistics.total += candidates0.size();
            thisThreadStatistics.capacity += candidates0.capacity();
            for(const Candidate& candidate: candidates0) {
                if(candidate.frequency >= minFrequency) {
                    ++thisThreadStatistics.highFrequency;
                }
            }
        }
    }
}

