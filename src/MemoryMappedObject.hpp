// Class to describe a single class object stored in a file mapped to memory.

#ifndef SHASTA_MEMORY_MAPPED_OBJECT_HPP
#define SHASTA_MEMORY_MAPPED_OBJECT_HPP

// Shasta.
#include "SHASTA_ASSERT.hpp"
#include "filesystem.hpp"
#include "touchMemory.hpp"

// Standard libraries.
#include "array.hpp"
#include <cstring>
#include "cstddef.hpp"
#include "iostream.hpp"
#include "stdexcept.hpp"
#include "string.hpp"

// Linux.
#include <fcntl.h>
#include <sys/mman.h>
#ifdef __linux__
#include <linux/mman.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>



namespace shasta {
    namespace MemoryMapped {
        template<class T> class Object;
    }
}



template<class T> class shasta::MemoryMapped::Object {
public:

    // Constructor and destructor.
    Object();
    ~Object();

    // Disallow copy and assignment.
    Object(const Object&) = delete;
    Object& operator=(const Object&) = delete;

    // Create a mapped object.
    void createNew(const string& name, size_t pageSize);

    // Open a previously created vector with read-only or read-write access.
    void accessExisting(const string& name, bool readWriteAccess);
    void accessExistingReadOnly(const string& name);
    void accessExistingReadWrite(const string& name);

    // Sync the mapped memory to disk.
    // This guarantees that the data on disk reflect all the latest changes in memory.
    // This is automatically called by close, and therefore also by the destructor.
    void syncToDisk();

    // Sync the mapped memory to disk, then unmap it.
    // This is automatically called by the destructor.
    void close();

    // Close and remove the supporting file.
    void remove();

    // Return a pointer to the stored object.
    T* operator->();
    const T* operator->() const;

private:

    // Compute the number of pages needed to hold n bytes.
    static size_t computePageCount(size_t n, size_t pageSize)
    {
        return (n - 1ULL ) / pageSize  + 1ULL;
    }

    // The header begins at the beginning of the mapped file.
    class Header {
        public:

        // The size of the header in bytes, including padding.
        size_t headerSize;

        // The size of the object, in bytes.
        size_t objectSize;

        // The number of objects. It is always 1.
        size_t objectCount;

        // The mapped file is always allocated with size equal to
        // a multiple of this page size, specified to createNew.
        // If the mapped file is backed by Linux huge pages,
        // this value must equal the Linux huge page size used
        // or a multiple of it.
        size_t pageSize;

        // The number of pages in the mapped file.
        // this equals exactly fileSize/pageSize.
        size_t pageCount;

        // The total number of allocated bytes in the mapped file.
        // This equals headerSize + dataSize, rounded up to the next
        // multiple of a page size.
        size_t fileSize;

        // The capacity is always 1.
        size_t capacity;

        // Magic number used for sanity check.
        static const size_t constantMagicNumber =  0xb7756f4515d8bc94ULL;
        size_t magicNumber;

        // Pad to 256 bytes to make sure the object is aligned with cache lines.
        array<size_t, 24> padding;

        // Constructor.
        Header(size_t pageSizeArgument)
        {
            clear();
            headerSize = sizeof(Header);
            objectSize = sizeof(T);
            objectCount = 1;
            pageSize = pageSizeArgument;
            pageCount = computePageCount(headerSize + objectSize * objectCount, pageSize);
            fileSize = pageCount * pageSize;
            capacity = 1;
            magicNumber = constantMagicNumber;
        }



        // Set the header to all zero bytes.
        void clear()
        {
            std::memset(this, 0, sizeof(Header));
        }
    };
    static_assert(sizeof(Header) == 256, "Unexpected header size for MemoryMapped::Object.");
    Header* header;

    // The data immediately follow the header.
    T* data;

    // Flags that indicate if the mapped file is open, and if so,
    // whether it is open for read-only or read-write.
public:
    bool isOpen;
    bool isOpenWithWriteAccess;
private:

    // The file name. If not open, this is an empty string.
    string fileName;

    // Unmap the memory.
    void unmap();

    // Some private utility functions;

    // Open the given file name as new (create if not existing, truncate if existing)
    // and with write access.
    // Return the file descriptor.
    static int openNew(const string& name);

    // Open the given existing file.
    // Return the file descriptor.
    static int openExisting(const string& name, bool readWriteAccess);

    // Truncate the given file descriptor to the specified size.
    static void truncate(int fileDescriptor, size_t fileSize);

    // Map to memory the given file descriptor for the specified size.
    static void* map(int fileDescriptor, size_t fileSize, bool writeAccess);

    // Find the size of the file corresponding to an open file descriptor.
    size_t getFileSize(int fileDescriptor);

    void createNewAnonymous(size_t pageSize);
};



// Default constructor.
template<class T> inline shasta::MemoryMapped::Object<T>::Object() :
    header(0),
    data(0),
    isOpen(false),
    isOpenWithWriteAccess(false)
{
}



// Destructor.
template<class T> inline shasta::MemoryMapped::Object<T>::~Object()
{
    if(isOpen) {
        close();
    }
}

// Open the given file name as new (create if not existing, truncate if existing)
// and with write access.
// Return the file descriptor.
template<class T> inline int shasta::MemoryMapped::Object<T>::openNew(const string& name)
{
    const int fileDescriptor = ::open(
            name.c_str(),
            O_CREAT | O_TRUNC | O_RDWR,
            S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if(fileDescriptor == -1) {
        throw runtime_error("Error " + to_string(errno)
            + " opening MemoryMapped::Object " + name + ": " + strerror(errno));
    }
    return fileDescriptor;
}

// Open the given existing file.
// Return the file descriptor.
template<class T> inline int shasta::MemoryMapped::Object<T>::openExisting(const string& name, bool readWriteAccess)
{
    const int fileDescriptor = ::open(
        name.c_str(),
        readWriteAccess ? O_RDWR : O_RDONLY);
    if(fileDescriptor == -1) {
        throw runtime_error("Error " + to_string(errno)
            + " opening MemoryMapped::Object " + name + ": " + strerror(errno));
    }
    return fileDescriptor;
}

// Truncate the given file descriptor to the specified size.
template<class T> inline void shasta::MemoryMapped::Object<T>::truncate(int fileDescriptor, size_t fileSize)
{
    const int ftruncateReturnCode = ::ftruncate(fileDescriptor, fileSize);
    if(ftruncateReturnCode == -1) {
        ::close(fileDescriptor);
        throw runtime_error(
        "The following error occurred during ftruncate:\n"
        "Error " + to_string(errno) + ": " +  string(strerror(errno)) +
        "\nThe most likely cause for this error is insufficient memory. "
        "Run on a larger machine.");
    }
}

// Map to memory the given file descriptor for the specified size.
template<class T> inline void* shasta::MemoryMapped::Object<T>::map(int fileDescriptor, size_t fileSize, bool writeAccess)
{
    void* pointer = ::mmap(0, fileSize, PROT_READ | (writeAccess ? PROT_WRITE : 0), MAP_SHARED, fileDescriptor, 0);
    if(pointer == reinterpret_cast<void*>(-1LL)) {
        ::close(fileDescriptor);
        throw runtime_error("Error during mmap.");
    }
    return pointer;
}

// Find the size of the file corresponding to an open file descriptor.
template<class T> inline size_t shasta::MemoryMapped::Object<T>::getFileSize(int fileDescriptor)
{
    struct stat fileInformation;
    const int fstatReturnCode = ::fstat(fileDescriptor, &fileInformation);
    if(fstatReturnCode == -1) {
        ::close(fileDescriptor);
        throw runtime_error("Error during fstat.");
    }
    return fileInformation.st_size;
}



// Create a new mapped object.
template<class T> inline void shasta::MemoryMapped::Object<T>::createNew(
    const string& name,
    size_t pageSize)
{
    SHASTA_ASSERT(pageSize==4096 || pageSize==2*1024*1024);

    if(name.empty()) {
        createNewAnonymous(pageSize);
        return;
    }

    try {
        // If already open, should have called close first.
        SHASTA_ASSERT(!isOpen);

        // Create the header.
        const Header headerOnStack(pageSize);
        const size_t fileSize = headerOnStack.fileSize;

        // Create the file.
        const int fileDescriptor = openNew(name);

        // Make it the size we want.
        truncate(fileDescriptor, fileSize);

        // Map it in memory.
        void* pointer = map(fileDescriptor, fileSize, true);

        // There is no need to keep the file descriptor open.
        // Closing the file descriptor as early as possible will make it possible to use large
        // numbers of Vector objects all at the same time without having to increase
        // the limit on the number of concurrently open descriptors.
        ::close(fileDescriptor);

        // Figure out where the data and the header go.
        header = static_cast<Header*>(pointer);
        data = reinterpret_cast<T*>(header+1);

        // Store the header.
        *header = headerOnStack;

        // Call the default constructor on the data.
        new(data) T();

        // Indicate that the mapped vector is open with write access.
        isOpen = true;
        isOpenWithWriteAccess = true;
        fileName = name;

    } catch(std::exception& e) {
        cout << e.what() << endl;
        throw runtime_error("Error " + to_string(errno)
            + " creating MemoryMapped::Object " + name + ": " + strerror(errno));
    }

}



// Create a new mapped object.
template<class T> inline void shasta::MemoryMapped::Object<T>::createNewAnonymous(
    size_t pageSize)
{
    try {
        // If already open, should have called close first.
        SHASTA_ASSERT(!isOpen);

        // Create the header.
        const Header headerOnStack(pageSize);
        const size_t fileSize = headerOnStack.fileSize;

        // Map it in memory.
        int flags = MAP_PRIVATE | MAP_ANONYMOUS;
#ifdef __linux__
        if(pageSize == 2*1024*1024) {
            flags |= MAP_HUGETLB | MAP_HUGE_2MB;
        }
#endif
        void* pointer = ::mmap(0, fileSize,
            PROT_READ | PROT_WRITE, flags,
            -1, 0);
        if(pointer == reinterpret_cast<void*>(-1LL)) {
            throw runtime_error("Error " + to_string(errno)
                + " during mmap call for MemoryMapped::Vector: " + string(strerror(errno)));
        }

        // Figure out where the data and the header go.
        header = static_cast<Header*>(pointer);
        data = reinterpret_cast<T*>(header+1);

        // Store the header.
        *header = headerOnStack;

        // Call the default constructor on the data.
        new(data) T();

        // Indicate that the mapped vector is open with write access.
        isOpen = true;
        isOpenWithWriteAccess = true;
        fileName = "";

    } catch(std::exception& e) {
        cout << e.what() << endl;
        throw runtime_error("Error " + to_string(errno)
            + " creating MemoryMapped::Object: " + strerror(errno));
    }

}



// Open a previously created object with read-only or read-write access.
template<class T> inline void shasta::MemoryMapped::Object<T>::accessExisting(const string& name, bool readWriteAccess)
{
    try {
        // If already open, should have called close first.
        SHASTA_ASSERT(!isOpen);

        // Create the file.
        const int fileDescriptor = openExisting(name, readWriteAccess);

        // Find the size of the file.
        const size_t fileSize = getFileSize(fileDescriptor);

        // Now map it in memory.
        void* pointer = map(fileDescriptor, fileSize, readWriteAccess);

        // There is no need to keep the file descriptor open.
        // Closing the file descriptor as early as possible will make it possible to use large
        // numbers of Vector objects all at the same time without having to increase
        // the limit on the number of concurrently open descriptors.
        ::close(fileDescriptor);

        // Figure out where the data and the header are.
        header = static_cast<Header*>(pointer);
        data = reinterpret_cast<T*>(header+1);

        // Sanity checks.
        if(header->magicNumber != Header::constantMagicNumber) {
            throw runtime_error("Error accessing " + name +
                ": unexpected magic number in header. " +
                "The binary format of this file is not recognized. " +
                "Perhaps a file mixup?"
                );
        }
        if(header->fileSize != fileSize) {
            throw runtime_error("Error accessing " + name +
                ": file size not consistent with file header. " +
                "Perhaps a file mixup?"
                );
        }
        if(header->objectSize != sizeof(T)) {
            throw runtime_error("Error accessing " + name +
                ": unexpected object size. Expected " + to_string(sizeof(T)) +
                ", found " + to_string(header->objectSize) +
                ". You may be attempting to access an assembly created by a different version of Shasta."
                );
        }

        // Indicate that the mapped vector is open.
        isOpen = true;
        isOpenWithWriteAccess = readWriteAccess;
        fileName = name;

    } catch(std::exception& e) {
        throw runtime_error("Error accessing " + name + ": " + e.what());
    }

}
template<class T> inline void shasta::MemoryMapped::Object<T>::accessExistingReadOnly(const string& name)
{
    accessExisting(name, false);
}
template<class T> inline void shasta::MemoryMapped::Object<T>::accessExistingReadWrite(const string& name)
{
    accessExisting(name, true);
}



// Sync the mapped memory to disk.
template<class T> inline void shasta::MemoryMapped::Object<T>::syncToDisk()
{
    SHASTA_ASSERT(isOpen);
    const int msyncReturnCode = ::msync(header, header->fileSize, MS_SYNC);
    if(msyncReturnCode == -1) {
        throw runtime_error("Error during msync for " + fileName);
    }
}

// Unmap the memory.
template<class T> inline void shasta::MemoryMapped::Object<T>::unmap()
{
    SHASTA_ASSERT(isOpen);

    const int munmapReturnCode = ::munmap(header, header->fileSize);
    if(munmapReturnCode == -1) {
        throw runtime_error("Error unmapping " + fileName);
    }

    // Mark it as not open.
    isOpen = false;
    isOpenWithWriteAccess = false;
    header = 0;
    data = 0;
    fileName = "";

}



// Sync the mapped memory to disk, then unmap it.
template<class T> inline void shasta::MemoryMapped::Object<T>::close()
{
    SHASTA_ASSERT(isOpen);

    if(!fileName.empty()) {
        syncToDisk();
    }
    unmap();
}

// Close it and remove the supporting file.
template<class T> inline void shasta::MemoryMapped::Object<T>::remove()
{
    const string savedFileName = fileName;
    close();    // This forgets the fileName.
    filesystem::remove(savedFileName);
}


// Return a pointer to the stored object.
template<class T> inline T* shasta::MemoryMapped::Object<T>::operator->()
{
    SHASTA_ASSERT(isOpen);
    SHASTA_ASSERT(data);
    return data;
}
template<class T> inline const T* shasta::MemoryMapped::Object<T>::operator->() const
{
    SHASTA_ASSERT(isOpen);
    SHASTA_ASSERT(data);
    return data;
}



#endif
