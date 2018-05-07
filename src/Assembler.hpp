#ifndef CZI_NANOPORE2_ASSEMBLER_HPP
#define CZI_NANOPORE2_ASSEMBLER_HPP

// Nanopore2
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
private:
    string smallDataFileNamePrefix;
    string largeDataFileNamePrefix;
    string smallDataName(const string& name) const
    {
        return smallDataFileNamePrefix + name;
    }
    string largeDataName(const string& name) const
    {
        return largeDataFileNamePrefix + name;
    }
    size_t smallDataPageSize;
    size_t largeDataPageSize;

};

#endif
