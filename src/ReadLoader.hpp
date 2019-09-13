#ifndef SHASTA_READ_LOADER_HPP
#define SHASTA_READ_LOADER_HPP

// shasta
#include "LongBaseSequence.hpp"
#include "MemoryMappedObject.hpp"
#include "MultitreadedObject.hpp"

// Standard library.
#include "memory.hpp"
#include "string.hpp"

namespace shasta {
    class ReadLoader;
}



// Class used to load reads from a fasta file.
class shasta::ReadLoader :
    public MultithreadedObject<ReadLoader>{
public:

    // The constructor does all the work.
    ReadLoader(
        const string& fileName,
        size_t minReadLength,
        size_t threadCountForReading,
        size_t threadCountForProcessing,
        const string& dataNamePrefix,
        size_t pageSize,
        LongBaseSequences& reads,
        MemoryMapped::VectorOfVectors<char, uint64_t>& readNames,
        MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& readRepeatCounts);

    // The number of reads and raw bases discarded because the read length
    // was less than minReadLength.
    uint64_t discardedShortReadReadCount;
    uint64_t discardedShortReadBaseCount;

    // The number of reads and raw bases discarded because the read
    // contained repeat counts greater than 255.
    uint64_t discardedBadRepeatCountReadCount;
    uint64_t discardedBadRepeatCountBaseCount;

private:

    // The file descriptor for the input file.
    int fileDescriptor = -1;

    // The size, in bytes, of the input file.
    size_t fileSize;
    void getFileSize();

    // The minimum read length. Shorter reads are not stored.
    size_t minReadLength;

    // Buffer to keep the the input file being processed.
    vector<char> buffer;

    // The number of threads to be used for read and for processing.
    size_t threadCountForReading;
    size_t threadCountForProcessing;

};



#endif
