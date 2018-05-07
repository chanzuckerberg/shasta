#ifndef CZI_NANOPORE2_READ_LOADER_HPP
#define CZI_NANOPORE2_READ_LOADER_HPP

// Nanopore2
#include "LongBaseSequence.hpp"
#include "MemoryMappedObject.hpp"
#include "MultitreadedObject.hpp"

// Standard library.
#include "string.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {
        class ReadLoader;
        class AssemblerInfo;
    }
}



// Class used to load reads from a fasta file.
class ChanZuckerberg::Nanopore2::ReadLoader :
    public MultithreadedObject<ReadLoader>{
public:

    // The constructor does all the work.
    ReadLoader(
        const string& fileName,
        size_t blockSize,
        size_t threadCountForReading,
        size_t threadCountForProcessing,
        LongBaseSequences& reads,
        MemoryMapped::VectorOfVectors<char, uint64_t>& readNames);
private:

    // The file descriptor for the input file.
    int fileDescriptor = -1;

    // The size, in bytes, of the input file.
    size_t fileSize;
    void getFileSize();

    // The block size we are using.
    size_t blockSize;

    // The begin/end offset of the block being processed.
    size_t blockBegin;
    size_t blockEnd;

    // Buffer to keep a block of the input file.
    vector<char> buffer;

    // Characters left over from the previous block.
    vector<char> leftOver;

    // The number of threads to be used for processing.
    size_t threadCountForProcessing;

    // Read one block into the above buffer.
    // This reads into the buffer the portion of the input file
    // at offset in [blockBegin, blockEnd).
    // It then reads the rest of the last read in the block
    // and updates blockEnd accordingly.
    void readBlock(size_t threadCount);
    void readBlockSequential();
    void readBlockParallel(size_t threadCount);

    // Function called by each thread.
    void threadFunction(size_t threadId);

    // Return true if a read begins at this position in the buffer.
    bool readBeginsHere(size_t bufferIndex) const;

};



#endif
