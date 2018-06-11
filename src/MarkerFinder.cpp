// Nanopore2.
#include "MarkerFinder.hpp"
#include "LongBaseSequence.hpp"
#include "ReadId.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Standard library.
#include <chrono>
#include <limits>


MarkerFinder::MarkerFinder(
    size_t k,
    const MemoryMapped::Vector<KmerInfo>& kmerTable,
    LongBaseSequences& reads,
    MemoryMapped::VectorOfVectors<CompressedMarker0, uint64_t>& markers,
    size_t threadCountArgument) :
    MultithreadedObject(*this),
    k(k),
    kmerTable(kmerTable),
    reads(reads),
    markers(markers),
    threadCount(threadCountArgument)
{
    // Initial message.
    cout << timestamp << "Finding markers in " << reads.size() << " reads." << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    const size_t batchSize = 100000;
    markers.beginPass1(reads.size());
    setupLoadBalancing(reads.size(), batchSize);
    pass = 1;
    runThreads(&MarkerFinder::threadFunction, threadCount, "threadLogs/MarkerFinder-Pass1");
    markers.beginPass2();
    markers.endPass2(false);
    setupLoadBalancing(reads.size(), batchSize);
    pass = 2;
    runThreads(&MarkerFinder::threadFunction, threadCount, "threadLogs/MarkerFinder-Pass2");

    // Final message.
    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    cout << timestamp << "Finding markers completed in " << tTotal << " s." << endl;
}



void MarkerFinder::threadFunction(size_t threadId)
{
    ostream& out = getLog(threadId);

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << "Working on reads ";
        out << begin << " through " << end << " of " << reads.size() << endl;

        // Loop over reads of this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            const LongBaseSequenceView read = reads[readId];
            size_t markerCount = 0; // For this read.
            CompressedMarker0* markerPointer = 0;
            if(pass == 2) {
                markerPointer = markers.begin(readId);
            }

            if(read.baseCount >= k) {   // Avoid pathological case.

                // Loop over k-mers of this read.
                Kmer kmer;
                for(size_t position=0; position<k; position++) {
                    kmer.set(position, read[position]);
                }
                uint32_t oldPosition = 0;
                for(uint32_t position=0; /*The check is done later */; position++) {
                    const KmerId kmerId = KmerId(kmer.id(k));
                    if(kmerTable[kmerId].isMarker) {
                        // This k-mer is a marker.

                        if(pass == 1) {
                            ++markerCount;
                        } else {
                            markerPointer->kmerId = kmerId;
                            const uint32_t shift = position - oldPosition;
                            CZI_ASSERT(shift <= std::numeric_limits<CompressedMarker0::Shift>::max());   // We may have to handle this.
                            markerPointer->shift = CompressedMarker0::Shift(shift);
                            ++markerPointer;
                            oldPosition = position;
                        }
                    }

                    if(position+k == read.baseCount) {
                        break;
                    }

                    // Update the k-mer.
                    kmer.shiftLeft();
                    kmer.set(k-1, read[position+k]);
                }
            }

            if(pass == 1) {
                markers.incrementCount(readId, markerCount);
            } else {
                CZI_ASSERT(markerPointer == markers.end(readId));
            }
        }
    }

}
