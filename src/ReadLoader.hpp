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
class CompressedRunnieReader;



// Class used to load reads from a fasta file.
class shasta::ReadLoader :
    public MultithreadedObject<ReadLoader>{
public:

    // The constructor does all the work.
    ReadLoader(
        const string& fileName,
        size_t minReadLength,
        size_t threadCount,
        const string& dataNamePrefix,
        size_t pageSize,
        LongBaseSequences& reads,
        MemoryMapped::VectorOfVectors<char, uint64_t>& readNames,
        MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& readRepeatCounts);

    // The number of reads and raw bases discarded because the read
    // contained invalid bases.
    uint64_t discardedInvalidBaseReadCount = 0;
    uint64_t discardedInvalidBaseBaseCount = 0; // Only counts the valid bases in those reads.

    // The number of reads and raw bases discarded because the read length
    // was less than minReadLength.
    uint64_t discardedShortReadReadCount = 0;
    uint64_t discardedShortReadBaseCount = 0;

    // The number of reads and raw bases discarded because the read
    // contained repeat counts greater than 255.
    uint64_t discardedBadRepeatCountReadCount = 0;
    uint64_t discardedBadRepeatCountBaseCount = 0;

private:

    // The name of the file we are processing.
    const string& fileName;

    // The minimum read length. Shorter reads are not stored.
    const size_t minReadLength;

    // The number of threads to be used for processing.
    // Reading is done single-threaded as there is usually no benefit
    // frm multithreaded reading.
    size_t threadCount;
    void adjustThreadCount();

    // Information that we can use to create temporary
    // memory mapped binary data structures.
    const string& dataNamePrefix;
    const size_t pageSize;

    // The data structure that the reads will be added to.
    LongBaseSequences& reads;
    MemoryMapped::VectorOfVectors<char, uint64_t>& readNames;
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& readRepeatCounts;

    // Create the name to be used for a MemoryMapped
    // object to be used by a thread.
    string threadDataName(
        size_t threadId,
        const string& dataName) const;

    // Read an entire file into a buffer,
    // using threadCountForReading threads.
    vector<char> buffer;
    void readFile();

    // Vectors where each thread stores the reads it found.
    // Indexed by threadId.
    vector< shared_ptr<MemoryMapped::VectorOfVectors<char, uint64_t> > > threadReadNames;
    vector< shared_ptr<LongBaseSequences> > threadReads;
    vector< shared_ptr<MemoryMapped::VectorOfVectors<uint8_t, uint64_t> > > threadReadRepeatCounts;
    void allocatePerThreadDataStructures();
    void allocatePerThreadDataStructures(size_t threadId);

    // Store the reads computed by each thread and free
    // the per-thread data structures.
    void storeReads();

    // Functions used for fasta files.
    void processFastaFile();
    void processFastaFileThreadFunction(size_t threadId);
    // Function that returns true if a read begins
    // at this position in Fasta format.
    bool fastaReadBeginsHere(uint64_t offset) const;



    // Functions used for compressed runnie files.
    void processCompressedRunnieFile();


};



#endif
