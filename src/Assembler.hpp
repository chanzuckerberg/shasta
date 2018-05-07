#ifndef CZI_NANOPORE2_ASSEMBLER_HPP
#define CZI_NANOPORE2_ASSEMBLER_HPP

// Nanopore2
#include "LongBaseSequence.hpp"
#include "MemoryMappedObject.hpp"
#include "MultitreadedObject.hpp"

// Standard library.
#include "string.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {
        class Assembler;
        class AssemblerInfo;
    }
}



// Class used to store various pieces of assembler information in shared memory.
class ChanZuckerberg::Nanopore2::AssemblerInfo {
public:

    // The length of k-mers used to define markers.
    size_t k;
};



class ChanZuckerberg::Nanopore2::Assembler :
    public MultithreadedObject<Assembler> {
public:



    // The constructor specifies the file name prefixes for binary data files.
    // There are two prefixes, one used for small data and one used for large data.
    // If these are directory names, they must include the final "/".
    // It also specifies the page size for small and large binary data files.
    // Typically, small binary data files will reside in a regular
    // directory on disk or on /dev/shm mapped backed by 4K pages,
    // while large binary data file will reside in a huge page
    // file system backed by 2MB pages.
    // 1GB huge pages are also supported.
    // The page sizes specified here must be equal to, or be an exact multiple of,
    // the actual size of the pages backing the data.
    // If the system has no large pages and it is not possible to change that,
    // use 4096 for both page sizes, and for performance place both
    // the small and large binary data under /dev/shm (in-memory filesystem).
    Assembler(
        const string& smallDataFileNamePrefix,
        const string& largeDataFileNamePrefix,
        size_t smallDataPageSize,
        size_t largeDataPageSize
        );


    // Add reads from a fasta file.
    // The reads are added to those already previously present.
    void addReadsFromFasta(
        const string& fileName,
        size_t blockSize,
        size_t threadCountForReading,
        size_t threadCountForProcessing);

private:

    // Data filled in by the constructor.
    string smallDataFileNamePrefix;
    string largeDataFileNamePrefix;
    size_t smallDataPageSize;
    size_t largeDataPageSize;

    // Functions to construct names for small and large binary objects.
    string smallDataName(const string& name) const
    {
        return smallDataFileNamePrefix + name;
    }
    string largeDataName(const string& name) const
    {
        return largeDataFileNamePrefix + name;
    }

    // Various pieces of assembler information stored in shared memory.
    // See class AssemblerInfo for more information.
    MemoryMapped::Object<AssemblerInfo> assemblerInfo;

    // The reads used for this assembly.
    // Indexed by ReadId.
    LongBaseSequences reads;

    // The names of the reads from the input fasta or fastq files.
    // Indexed by ReadId.
    // Note that we don't enforce uniqueness of read names.
    // We don't use read names to identify reads.
    // These names are only used as an aid in tracing each read
    // back to its origin.
    MemoryMapped::VectorOfVectors<char, uint64_t> readNames;

};

#endif
