// Shasta.
#include "LowHash0.hpp"
#include "ReadFlags.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Standard library.
#include "chrono.hpp"
#include <numeric>



// Class LowHash0 uses the LowHash0 algorithm to find candidate pairs
// of aligned reads. It uses as features
// sequences of m consecutive markers.



LowHash0::LowHash0(
    size_t m,                       // Number of consecutive markers that define a feature.
    double hashFraction,

    // Iteration control. See MinHashOptions for details.
    size_t minHashIterationCount,
    double alignmentCandidatesPerRead,

    size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
    size_t minBucketSize,           // The minimum size for a bucket to be used.
    size_t maxBucketSize,           // The maximum size for a bucket to be used.
    size_t minFrequency,            // Minimum number of minHash hits for a pair to be considered a candidate.
    size_t threadCountArgument,
    const MemoryMapped::Vector<KmerInfo>& kmerTable,
    const MemoryMapped::Vector<ReadFlags>& readFlags,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    MemoryMapped::Vector<OrientedReadPair>& candidateAlignments,
    MemoryMapped::Vector< array<uint64_t, 3> >& readLowHashStatistics,
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
    readLowHashStatistics(readLowHashStatistics),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize),
    histogramCsv("LowHashBucketHistogram.csv")

{
    cout << timestamp << "LowHash0 begins." << endl;
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
    const uint32_t leadingZeroBitCount = uint32_t(__builtin_clzl(totalLowHashCountEstimate));
    const uint32_t log2TotalLowHashCountEstimate = 64 - leadingZeroBitCount;


    // If log2MinHashBucketCount is 0, choose a reasonable value
    // for the current number of reads.
    // Otherwise, check that log2MinHashBucketCount is not unreasonably small.
    if(log2MinHashBucketCount == 0) {
        log2MinHashBucketCount = 5 + log2TotalLowHashCountEstimate;
    } else {
        if(log2MinHashBucketCount < log2TotalLowHashCountEstimate) {
            throw runtime_error("log2MinHashBucketCount is unreasonably small.");
        }
    }
    if(log2MinHashBucketCount > 31) {

        cout << "log2MinHashBucketCount reduced from " << log2MinHashBucketCount <<
            " to maximum allowed value 31."  << endl;
        log2MinHashBucketCount = 31;
    }
    const uint32_t bucketCount = 1 << log2MinHashBucketCount;
    mask = bucketCount - 1;
    cout << "LowHash0 algorithm will use 2^" << log2MinHashBucketCount;
    cout << " = " << bucketCount << " buckets. "<< endl;




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
    

    // Set up work areas.
    buckets.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "tmp-LowHash0-Buckets"),
        largeDataPageSize);
    lowHashes.resize(orientedReadCount);
    candidates.resize(readCount);
    threadStatistics.resize(threadCount);
    readLowHashStatistics.resize(readCount);
    fill(readLowHashStatistics.begin(), readLowHashStatistics.end(),
        array<uint64_t, 3>({0, 0, 0}));

    // Write the header of the histogram file.
    histogramCsv << "Iteration,BucketSize,BucketCount,FeatureCount\n";



    // LowHash0 iteration loop.
    uint64_t highFrequency = 0;
    for(iteration=0; ; iteration++) {

        // Decide if it is time to stop the iteration.
        if(minHashIterationCount==0) {

            // minHashIterationCount is zero, so alignmentCandidatesPerRead
            // controls the iteration.
            const double currentAlignmentCandidatesPerRead =
                2. * double(highFrequency) / double(readCount);
            if(iteration != 0) {
                cout << "Average number of alignment candidates that each read is involved in is " <<
                    currentAlignmentCandidatesPerRead << endl;
            }
            if(currentAlignmentCandidatesPerRead >= alignmentCandidatesPerRead) {
                break;
            }

        } else {

            // minHashIterationCount is not zero, so it controls the number of iterations.
            if(iteration==minHashIterationCount) {
                break;
            }
        }



        cout << timestamp << "LowHash0 iteration " << iteration << " begins." << endl;

        // Pass1: compute the low hashes for each oriented read
        // and prepare the buckets for filling.
        buckets.clear();
        buckets.beginPass1(bucketCount);
        size_t batchSize = 10000;
        setupLoadBalancing(readCount, batchSize);
        runThreads(&LowHash0::pass1ThreadFunction, threadCount);

        // Pass 2: fill the buckets.
        buckets.beginPass2();
        batchSize = 10000;
        setupLoadBalancing(readCount, batchSize);
        runThreads(&LowHash0::pass2ThreadFunction, threadCount);
        buckets.endPass2(false, false);
        computeBucketHistogram();

        // Pass 3: inspect the buckets to find candidates.
        batchSize = 10000;
        setupLoadBalancing(readCount, batchSize);
        runThreads(&LowHash0::pass3ThreadFunction, threadCount);

        // Write a summary for this iteration.
        highFrequency = 0;
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
    SHASTA_ASSERT(orientedReadCount == 2*readCount);
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        const auto& candidates0 = candidates[readId0];
        for(const Candidate& candidate: candidates0) {
            if(candidate.frequency >= minFrequency) {
                const ReadId readId1 = candidate.readId1;
                SHASTA_ASSERT(readId0 < readId1);
                candidateAlignments.push_back(
                    OrientedReadPair(readId0, readId1, candidate.strand==0));
            }
        }
    }
    cout << "Found " << candidateAlignments.size() << " alignment candidates."<< endl;
    cout << "Average number of alignment candidates per oriented read is ";
    cout << (2.* double(candidateAlignments.size())) / double(orientedReadCount)  << "." << endl;

    // Write read bucket statistics.
    ofstream csv("ReadLowHashStatistics.csv");
    csv << "ReadId,Palindromic,Features,Sparse,Good,Crowded,Total,FeatureSampling,"
        "SparseFraction,GoodFraction,CrowdedFraction\n";
    for(ReadId readId=0; readId<readCount; readId++) {
        const array<uint64_t, 3>& counters = readLowHashStatistics[readId];
        const uint64_t total = std::accumulate(counters.begin(), counters.end(), 0);
        const uint64_t featureCount = markers.size(OrientedReadId(readId, 0).getValue()) - (m-1);
        const double featureSampling = double(total) / double(featureCount);
        csv << readId << ",";
        csv << (readFlags[readId].isPalindromic ? "Yes," : "No,");
        csv << featureCount << ",";
        csv << counters[0] << ",";
        csv << counters[1] << ",";
        csv << counters[2] << ",";
        csv << total << ",";
        csv << featureSampling << ",";
        if(total == 0) {
            csv << ",,\n";
        } else {
            csv << double(counters[0]) / double(total) << ",";
            csv << double(counters[1]) / double(total) << ",";
            csv << double(counters[2]) / double(total) << "\n";
        }
    }



    // Clean up work areas.
    buckets.remove();
    kmerIds.remove();



    // Done.
    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "LowHash0 completed in " << tTotal << " s." << endl;
}



void LowHash0::createKmerIds()
{
    kmerIds.createNew(
    	largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "tmp-LowHash0-Markers"),
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
    runThreads(&LowHash0::createKmerIds, threadCount);
}



// Thread function for createKmerIds.
void LowHash0::createKmerIds(size_t threadId)
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



// Pass1: compute the low hashes for each oriented read
// and prepare the buckets for filling.
void LowHash0::pass1ThreadFunction(size_t threadId)
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

                vector<uint64_t>& orientedReadLowHashes = lowHashes[orientedReadId.getValue()];
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
                        orientedReadLowHashes.push_back(hash);
                        const uint64_t bucketId = hash & mask;
                        buckets.incrementCountMultithreaded(bucketId);
                    }
                }
            }
        }
    }

}



// Pass 2: fill the buckets.
void LowHash0::pass2ThreadFunction(size_t threadId)
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
                const vector<uint64_t>& orientedReadLowHashes = lowHashes[orientedReadId.getValue()];

                for(const uint64_t hash: orientedReadLowHashes) {
                    const uint64_t bucketId = hash & mask;
                    buckets.storeMultithreaded(bucketId, BucketEntry(orientedReadId, hash));

                    // Update statistics for this read.
                    const uint64_t bucketSize = buckets.size(bucketId);
                    if(bucketSize < minBucketSize) {
                        ++readLowHashStatistics[readId][0];
                    } else if(bucketSize > maxBucketSize) {
                        ++readLowHashStatistics[readId][2];
                    } else {
                        ++readLowHashStatistics[readId][1];
                    }
                }
            }
        }
    }
}



// Pass 3: inspect the buckets to find candidates.
void LowHash0::pass3ThreadFunction(size_t threadId)
{

    // The alignment candidates found at this iteration for a single read.
    vector<Candidate> newCandidates;

    // The merged candidates for a single readId0.
    vector<Candidate> mergedCandidates;

    ThreadStatistics& thisThreadStatistics = threadStatistics[threadId];
    thisThreadStatistics.clear();

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads assigned to this batch.
        for(ReadId readId0=ReadId(begin); readId0!=ReadId(end); readId0++) {
            newCandidates.clear();

            // Loop over two strands.
            for(Strand strand0=0; strand0<2; strand0++) {
                const OrientedReadId orientedReadId0(readId0, strand0);
                const OrientedReadId::Int orientedReadIdInt0 = orientedReadId0.getValue();

                // Loop over the low hashes for this oriented read.
                const vector<uint64_t>& orientedReadLowHashes = lowHashes[orientedReadIdInt0];
                for(const uint64_t hash: orientedReadLowHashes) {
                    const uint32_t hashHighBits = uint32_t(hash >> 32);

                    // Loop over oriented read ids in the bucket corresponding to this hash.
                    const uint64_t bucketId = hash & mask;
                    const span<BucketEntry> bucket = buckets[bucketId];
                    if(bucket.size() < max(size_t(2), minBucketSize)) {
                        continue;
                    }
                    if(bucket.size() > maxBucketSize) {
                        continue;   // The bucket is too big. Skip it.
                    }
                    for(const BucketEntry& bucketEntry: bucket) {
                        if(bucketEntry.hashHighBits != hashHighBits) {
                            continue;   // Collision.
                        }
                        const OrientedReadId orientedReadId1 = bucketEntry.orientedReadId;
                        const ReadId readId1 = orientedReadId1.getReadId();

                        // Only consider it if readId1 > readId0.
                        if(readId1 <= readId0) {
                            continue;
                        }

                        // Add it to our work area.
                        const bool isSameStrand = orientedReadId1.getStrand() == strand0;
                        newCandidates.push_back(Candidate(readId1, isSameStrand? 0 : 1));
                    }
                }
            }

            // Sort the candidates found during this iteration.
            sort(newCandidates.begin(), newCandidates.end());

            // Merge the contents of the work area
            // with the candidates previously stored.
            vector<Candidate>& storedCandidates = candidates[readId0];
            mergedCandidates.clear();
            merge(storedCandidates, newCandidates, mergedCandidates);

            // Store the merged candidates in place of the old ones.
            storedCandidates.resize(mergedCandidates.size());
            copy(mergedCandidates.begin(), mergedCandidates.end(), storedCandidates.begin());

            // Update thread statistics.
            thisThreadStatistics.total += storedCandidates.size();
            thisThreadStatistics.capacity += storedCandidates.capacity();
            for(const Candidate& candidate: storedCandidates) {
                if(candidate.frequency >= minFrequency) {
                    ++thisThreadStatistics.highFrequency;
                }
            }
        }
    }
}



// Merge two sorted vectors of candidates into a third one.
// The input vectors can be sorted but can have duplicates.
// During merging, when two candidates with the same readId1
// and strand are found, they are combined, adding up their frequency.
// This is used by pass4ThreadFunction.
void LowHash0::merge(
    const vector<LowHash0::Candidate>& x0,
    const vector<LowHash0::Candidate>& x1,
    vector<LowHash0::Candidate>& y
    )
{
    using Iterator = vector<LowHash0::Candidate>::const_iterator;
    Iterator begin0 = x0.begin();
    Iterator begin1 = x1.begin();
    Iterator end0 = x0.end();
    Iterator end1 = x1.end();



    // Merge loop..
    // At each step, we find the lowest of the two candidates
    // currently pointed by the two iterators (ties allowed and are ok).
    // If the merged vector is not empty and ends with an entry that compares equal,
    // increment its frequency. Otherwise, just create a new entry in the merged vector.
    Iterator it0 = begin0;
    Iterator it1 = begin1;
    while(true) {

        // If we reached the end of x0, process the remaining portion of x1,
        // then exit the merge loop.
        if(it0 == end0) {
            for(; it1!=end1; ++it1) {
                if(!y.empty() && y.back() == *it1) {
                    y.back().frequency = uint16_t(y.back().frequency + it1->frequency);
                } else {
                    y.push_back(*it1);
                }
            }
            break;
        }

        // If we reached the end of x1, process the remaining portion of x0,
        // then exit the merge loop.
        if(it1 == end1) {
            for(; it0!=end0; ++it0) {
                if(!y.empty() && y.back() == *it0) {
                    y.back().frequency = uint16_t(y.back().frequency + it0->frequency);
                } else {
                    y.push_back(*it0);
                }
            }
            break;
        }

        // If the current x0 entry is less than the current x1 entry,
        // process the x0 entry.
        if(*it0 < *it1) {
            if(!y.empty() && y.back() == *it0) {
                y.back().frequency  = uint16_t(y.back().frequency + it0->frequency);
            } else {
                y.push_back(*it0);
            }
            ++it0;
        } else {
        // If the current x0 entry is not less than the current x1 entry,
        // process the current x1 entry.
            if(!y.empty() && y.back() == *it1) {
                y.back().frequency = uint16_t(y.back().frequency + it1->frequency);
            } else {
                y.push_back(*it1);
            }
            ++it1;
        }
    }
}



void LowHash0::computeBucketHistogram()
{
    threadBucketHistogram.clear();
    threadBucketHistogram.resize(threadCount);
    const uint64_t batchSize = 10000;
    setupLoadBalancing(buckets.size(), batchSize);
    runThreads(&LowHash0::computeBucketHistogramThreadFunction, threadCount);

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
void LowHash0::computeBucketHistogramThreadFunction(size_t threadId)
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
