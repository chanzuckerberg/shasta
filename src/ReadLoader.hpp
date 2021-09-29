#ifndef SHASTA_READ_LOADER_HPP
#define SHASTA_READ_LOADER_HPP

// shasta
#include "PalindromeQuality.hpp"
#include "LongBaseSequence.hpp"
#include "MemoryMappedObject.hpp"
#include "MultithreadedObject.hpp"
#include "Reads.hpp"

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
        uint64_t minReadLength,
        bool noCache,
        size_t threadCount,
        const string& dataNamePrefix,
        size_t pageSize,
        bool detectPalindromesOnFastqLoad,
        double qScoreRelativeMeanDifference,
        double qScoreMinimumMean,
        double qScoreMinimumVariance,
        bool writePalindromesToCsv,
        Reads& reads);

    ~ReadLoader();
    
    // The number of reads and raw bases discarded because the read
    // contained invalid bases.
    uint64_t discardedInvalidBaseReadCount = 0;
    uint64_t discardedInvalidBaseBaseCount = 0; // Only counts the valid bases in those reads.

    // The number of reads and raw bases discarded because the read length
    // was less than minReadLength.
    uint64_t discardedShortReadReadCount = 0;
    uint64_t discardedShortReadBaseCount = 0;

    // The number of reads and raw bases discarded because the read had
    // a q score distribution that was indicative of a palindrome.
    uint64_t discardedPalindromicReadCount = 0;
    uint64_t discardedPalindromicBaseCount = 0;

    // The number of reads and raw bases discarded because the read
    // contained repeat counts greater than 255.
    uint64_t discardedBadRepeatCountReadCount = 0;
    uint64_t discardedBadRepeatCountBaseCount = 0;

private:

    // The name of the file we are processing.
    const string& fileName;

    // The minimum read length. Shorter reads are not stored.
    const uint64_t minReadLength;

    // If set, use the O_DIRECT flag when opening input files (Linux only).
    bool noCache;

    // The number of threads to be used for processing.
    // Reading is done single-threaded as there is usually no benefit
    // frm multithreaded reading.
    size_t threadCount;
    void adjustThreadCount();

    // Information that we can use to create temporary
    // memory mapped binary data structures.
    const string& dataNamePrefix;
    const size_t pageSize;

    // Boolean switch to use quality scores to skip reads that have an indication of palindromic sequence
    bool detectPalindromesOnFastqLoad;

    // Each of the 3 thresholds necessary for calling isPalindromic()
    double qScoreRelativeMeanDifference;
    double qScoreMinimumMean;
    double qScoreMinimumVariance;

    // This is true if shasta was run with command "filterReads"
    bool writePalindromesToCsv;

    // The data structure that the reads will be added to.
    Reads& reads;

    // Create the name to be used for a MemoryMapped object.
    string dataName(
        const string& dataName) const;
    string threadDataName(
        size_t threadId,
        const string& dataName) const;

    // Read an entire file into a buffer,
    // using threadCountForReading threads.
    int64_t fileSize;
    MemoryMapped::Vector<char> buffer;
    void allocateBuffer();
    void readFile();
    void allocateBufferAndReadFile();

    // Vectors where each thread stores the reads it found.
    // Indexed by threadId.
    vector< unique_ptr<MemoryMapped::VectorOfVectors<char, uint64_t> > > threadReadNames;
    vector< unique_ptr<MemoryMapped::VectorOfVectors<char, uint64_t> > > threadReadMetaData;
    vector< unique_ptr<LongBaseSequences> > threadReads;
    vector< unique_ptr<MemoryMapped::VectorOfVectors<uint8_t, uint64_t> > > threadReadRepeatCounts;
    vector< vector<string> > threadPalindromicReadNames;
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

    // Functions used for fastq files.
    void processFastqFile();
    void processFastqFileThreadFunction(size_t threadId);

    // Find all line ends in the file.
    void findLineEnds();
    void findLineEndsThreadFunction(size_t threadId);
    vector< vector<uint64_t> > threadLineEnds;
    vector<uint64_t> lineEnds;

#ifdef __linux__
    int tryDirectIO(const string& fileName);
#endif
};



#endif
