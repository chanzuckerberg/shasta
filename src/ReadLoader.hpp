#ifndef SHASTA_READ_LOADER_HPP
#define SHASTA_READ_LOADER_HPP

// shasta
#include "LongBaseSequence.hpp"
#include "MemoryMappedObject.hpp"
#include "MultitreadedObject.hpp"

// Standard library.
#include "memory.hpp"
#include "string.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class ReadLoader;
    }
}



// Class used to load reads from a fasta file.
class ChanZuckerberg::shasta::ReadLoader :
    public MultithreadedObject<ReadLoader>{
public:

    // The constructor does all the work.
    ReadLoader(
        const string& fileName,
        size_t minReadLength,
        size_t blockSize,
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

    // The block size we are using.
    size_t blockSize;

    // The begin/end offset of the block being read.
    size_t blockBegin;
    size_t blockEnd;

    // Buffer to keep the the input file being processed.
    vector<char> buffer;

    // Characters left over from the previous block.
    vector<char> leftOver;

    // The number of threads to be used for read and for processing.
    size_t threadCountForReading;
    size_t threadCountForProcessing;

    // Read one block into the above buffer.
    // This copies the leftOver data to buffer, then
    // reads into the rest of the buffer the portion of the input file
    // at offset in [blockBegin, blockEnd).
    // It finally moves to the leftOver data the
    // final, possibly partial, read in the buffer
    // (this is not done for the final block).
    void readBlock(size_t threadCount);
    void readBlockSequential();
    void readBlockParallel(size_t threadCount);

    // Functions called by each thread.
    void readThreadFunction(size_t threadId);
    void processThreadFunction(size_t threadId);

    // Return true if a read begins at this position in the buffer.
    bool readBeginsHere(size_t bufferIndex) const;

    // Vectors where each thread stores the reads it found.
    // Indexed by threadId.
    vector< shared_ptr<MemoryMapped::VectorOfVectors<char, uint64_t> > > threadReadNames;
    vector< shared_ptr<LongBaseSequences> > threadReads;
    vector< shared_ptr<MemoryMapped::VectorOfVectors<uint8_t, uint64_t> > > threadReadRepeatCounts;

    // Given the raw representation of a read, compute its
    // run-length representation.
    // This returns false if the read contains a homopolymer run
    // of more than 255 bases, which cannot be represented
    // with a one-byte repeat count.
    static bool computeRunLengthRead(
        const vector<Base>& read,
        vector<Base>& runLengthRead,
        vector<uint8_t>& readRepeatCount);

    // Create the name to be used for a MemoryMapped
    // object to be used by a thread.
    static string threadDataName(
        const string& dataNamePrefix,
        size_t threadId,
        const string& dataName);

};



#endif
