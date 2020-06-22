// Class to describe a vector stored in a file mapped to memory.

#ifndef SHASTA_MEMORY_MAPPED_VECTOR_HPP
#define SHASTA_MEMORY_MAPPED_VECTOR_HPP

// Shasta.
#include "SHASTA_ASSERT.hpp"
#include "filesystem.hpp"
#include "MarkerInterval.hpp"
#include "MurmurHash2.hpp"
#include "touchMemory.hpp"

// Boost libraries.
#include <boost/lexical_cast.hpp>

// Standard libraries.
#include <cstring>
#include "algorithm"
#include "cstddef.hpp"
#include "iostream.hpp"
#include "stdexcept.hpp"
#include "string.hpp"
#include "vector.hpp"

// Linux.
#include <fcntl.h>
#include <sys/mman.h>
#ifdef __linux__
#include <linux/mman.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>
#include "array.hpp"


// Forward declarations.
namespace shasta {
    namespace MemoryMapped {
        template<class T> class Vector;
    }
    void testMemoryMappedVector();
}



template<class T> class shasta::MemoryMapped::Vector {
public:

    // The access functions work as in std::vector.
    size_t size() const;
    bool empty() const;
    size_t capacity() const;
    T& operator[](size_t);
    const T& operator[](size_t) const;
    T& front();
    const T& front() const;
    T& back();
    const T& back() const;
    T* begin();
    const T* begin() const;
    T* end();
    const T* end() const;
    void push_back(const T&);

    // Constructor and destructor.
    Vector();
    ~Vector();

    // Disallow C++ copy and assignment
    // (see below for how to make a copy).
    Vector(const Vector&) = delete;
    Vector& operator=(const Vector&) = delete;

    // Create a new mapped vector with n objects.
    // The last argument specifies the required capacity.
    // Actual capacity will be a bit larger due to rounding up to the next page boundary.
    // The vector is stored in a memory mapped file with the specified name.
    void createNew(const string& name, size_t pageSize, size_t n=0, size_t requiredCapacity=0);

    // Open a previously created vector with read-only or read-write access.
    // If accessExistingReadWrite is called with allowReadOnly=true,
    // it attempts to open with read-write access, but if that fails falls back to
    // read-only access.
    void accessExisting(const string& name, bool readWriteAccess);
    void accessExistingReadOnly(const string& name);
    void accessExistingReadWrite(const string& name);

    // Attempt to open a previously created name with read write access.
    // If this fails, create the vector, empty, and access it with read write access.
    void accessExistingReadWriteOrCreateNew(const string& name, size_t pageSize);

    // Sync the mapped memory to disk.
    // This guarantees that the data on disk reflect all the latest changes in memory.
    // This is automatically called by close, and therefore also by the destructor.
    void syncToDisk();

    // Sync the mapped memory to disk, then unmap it.
    // This is automatically called by the destructor.
    void close();


    // Close and remove the supporting file.
    void remove();

    // Resize works as for std::vector;
    void resize(size_t);
    void clear()
    {
        resize(0);
    }

    // Touch a range of memory in order to cause the
    // supporting pages of virtual memory to be loaded in real memory.
    // The return value can be ignored.
    size_t touchMemory() const
    {
        return shasta::touchMemory(begin(), end());
    }


    void reserve();
    void reserve(size_t capacity);

    void unreserve();

    // Use this instead of resize when it is known that the size
    // will not further increase. This results in reduce
    // memory requirement, because resize increases capacity to 1.5
    // times the new size.
    void reserveAndResize(size_t n)
    {
        reserve(n);
        resize(n);
    }

    // Make a copy of the Vector.
    void makeCopy(Vector<T>& copy, const string& newName) const;

    // Comparison operator.
    bool operator==(const Vector<T>& that) const
    {
        return
            size() == that.size() &&
            std::equal(begin(), end(), that.begin());
    }

    // Return a hash function of the stored data.
    // Can be used to check integrity.
    uint64_t hash() const;

    // Rename the supporting memory mapped file, if any.
    void rename(const string& newFileName);

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

        // The size of each object stored in the vector, in bytes.
        size_t objectSize;

        // The number of objects currently stored in the vector.
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

        // The current capacity of the vector (number of objects that can be stored
        // in the currently allocated memory).
        size_t capacity;

        // Magic number used for sanity check.
        static const size_t constantMagicNumber =  0xa3756fd4b5d8bcc1ULL;
        size_t magicNumber;

        // Pad to 4096 bytes to make sure the data are page aligned.
        array<size_t, 4096/sizeof(size_t) - 8 > padding;



        // Constructor with a given size and capacity.
        // Actual capacity will a bit larger, rounded up to the next oage boundary.
        Header(size_t n, size_t requestedCapacity, size_t pageSizeArgument)
        {
            SHASTA_ASSERT(requestedCapacity >= n);
            clear();
            headerSize = sizeof(Header);
            objectSize = sizeof(T);
            objectCount = n;
            pageSize = pageSizeArgument;
            pageCount = computePageCount(headerSize + objectSize * requestedCapacity, pageSize);
            fileSize = pageCount * pageSize;
            capacity = (fileSize - headerSize) / objectSize;
            magicNumber = constantMagicNumber;
        }



        // Set the header to all zero bytes.
        void clear()
        {
            std::memset(this, 0, sizeof(Header));
        }

    };
    static_assert(sizeof(Header) == 4096, "Unexpected header size for MemoryMapped::Vector.");
    Header* header;

    // The data immediately follow the header.
    T* data;

public:

    // Flags that indicate if the mapped file is open, and if so,
    // whether it is open for read-only or read-write.
    bool isOpen;
    bool isOpenWithWriteAccess;

    // The file name. If not open, this is an empty string.
    string fileName;

private:
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


    void createNewAnonymous(size_t pageSize, size_t n=0, size_t requiredCapacity=0);
    void resizeAnonymous(size_t newSize);
    void reserveAnonymous(size_t newSize);
    void unmapAnonymous();
};



// Access functions for class Vector.
// As a result, invalid uses will result in segmentation faults.
// Note that in the non-const access functions we assert for isOpen, not isOpenWithWriteAccess.
// This is necessary to allow legimitate patterns, such as having a non-const reference
// to a Vector that is open read-only.
template<class T> inline size_t shasta::MemoryMapped::Vector<T>::size() const
{
    return isOpen ? header->objectCount : 0ULL;
}
template<class T> inline bool shasta::MemoryMapped::Vector<T>::empty() const
{
    return isOpen ? (size()==0) : 0ULL;
}
template<class T> inline size_t shasta::MemoryMapped::Vector<T>::capacity() const
{
    return isOpen ? header->capacity : 0ULL;
}

template<class T> inline T& shasta::MemoryMapped::Vector<T>::operator[](size_t i)
{
    // SHASTA_ASSERT(isOpen);
    return data[i];
}
template<class T> inline const T& shasta::MemoryMapped::Vector<T>::operator[](size_t i) const
{
    // SHASTA_ASSERT(isOpen);
    return data[i];
}

template<class T> inline T& shasta::MemoryMapped::Vector<T>::front()
{
    SHASTA_ASSERT(isOpen);
    return *data;
}
template<class T> inline const T& shasta::MemoryMapped::Vector<T>::front() const
{
    SHASTA_ASSERT(isOpen);
    SHASTA_ASSERT(size() > 0);
    return *data;
}

template<class T> inline T& shasta::MemoryMapped::Vector<T>::back()
{
    SHASTA_ASSERT(isOpen);
    return data[size() - 1ULL];
}
template<class T> inline const T& shasta::MemoryMapped::Vector<T>::back() const
{
    SHASTA_ASSERT(isOpen);
    return data[size() - 1ULL];
}
template<class T> inline T* shasta::MemoryMapped::Vector<T>::begin()
{
    SHASTA_ASSERT(isOpen);
    return data;
}
template<class T> inline const T* shasta::MemoryMapped::Vector<T>::begin() const
{
    SHASTA_ASSERT(isOpen);
    return data;
}

template<class T> inline T* shasta::MemoryMapped::Vector<T>::end()
{
    SHASTA_ASSERT(isOpen);
    return data + size();
}

template<class T> inline const T* shasta::MemoryMapped::Vector<T>::end() const
{
    SHASTA_ASSERT(isOpen);
    return data + size();
}
template<class T> inline void shasta::MemoryMapped::Vector<T>::push_back(const T& t)
{
    SHASTA_ASSERT(isOpen);
    resize(size()+1ULL);
    back() = t;

}


// Destructor.
template<class T> inline shasta::MemoryMapped::Vector<T>::~Vector()
{
    if(isOpen) {

        if(fileName.empty()) {

            unmapAnonymous();

        } else {

            if(isOpenWithWriteAccess) {
                unreserve();
            }
            close();
        }
    }
}


// Default constructor.
template<class T> inline shasta::MemoryMapped::Vector<T>::Vector() :
    header(0),
    data(0),
    isOpen(false),
    isOpenWithWriteAccess(false)
{
}



// Open the given file name as new (create if not existing, truncate if existing)
// and with write access.
// Return the file descriptor.
template<class T> inline int shasta::MemoryMapped::Vector<T>::openNew(const string& name)
{

    // The specified name is not a directory.
    // Open or create a file with this name.
    const int fileDescriptor = ::open(
            name.c_str(),
            O_CREAT | O_TRUNC | O_RDWR,
            S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if(fileDescriptor == -1) {
        throw runtime_error("Error opening " + name);
    }

    return fileDescriptor;
}

// Open the given existing file.
// Return the file descriptor.
template<class T> inline int shasta::MemoryMapped::Vector<T>::openExisting(const string& name, bool readWriteAccess)
{
    const int fileDescriptor = ::open(
        name.c_str(),
        readWriteAccess ? O_RDWR : O_RDONLY);
    if(fileDescriptor == -1) {
        throw runtime_error("Error " + boost::lexical_cast<string>(errno)
            + " opening MemoryMapped::Vector " + name + ": " + string(strerror(errno)));
    }
    return fileDescriptor;
}

// Truncate the given file descriptor to the specified size.
template<class T> inline void shasta::MemoryMapped::Vector<T>::truncate(int fileDescriptor, size_t fileSize)
{
    const int ftruncateReturnCode = ::ftruncate(fileDescriptor, fileSize);
    if(ftruncateReturnCode == -1) {
        ::close(fileDescriptor);
        throw runtime_error(
            "The following error occurred during ftruncate to size " + to_string(fileSize) +
            ":\nError " + to_string(errno) + ": " +  string(strerror(errno)) +
            "\nThe most likely cause for this error is insufficient memory. "
            "Run on a larger machine.");
    }
}

// Map to memory the given file descriptor for the specified size.
template<class T> inline void* shasta::MemoryMapped::Vector<T>::map(int fileDescriptor, size_t fileSize, bool writeAccess)
{
    void* pointer = ::mmap(0, fileSize, PROT_READ | (writeAccess ? PROT_WRITE : 0), MAP_SHARED, fileDescriptor, 0);
    if(pointer == reinterpret_cast<void*>(-1LL)) {
        ::close(fileDescriptor);
        if(errno == ENOMEM) {
            throw runtime_error("Memory allocation failure "
                "during mmap call for MemoryMapped::Vector.\n"
                "This assembly requires more memory than available.\n"
                "Rerun on a larger machine.");
        } else {
            throw runtime_error("Error " + boost::lexical_cast<string>(errno)
                + " during mremap call for MemoryMapped::Vector: " + string(strerror(errno)));
        }
    }
    return pointer;
}

// Find the size of the file corresponding to an open file descriptor.
template<class T> inline size_t shasta::MemoryMapped::Vector<T>::getFileSize(int fileDescriptor)
{
    struct stat fileInformation;
    const int fstatReturnCode = ::fstat(fileDescriptor, &fileInformation);
    if(fstatReturnCode == -1) {
        ::close(fileDescriptor);
        throw runtime_error("Error during fstat.");
    }
    return fileInformation.st_size;
}



// Create a new mapped vector with n objects.
// The last argument specifies the required capacity.
// Actual capacity will be a bit larger due to rounding up to the next page boundary.
// The vector is stored in a memory mapped file with the specified name.
template<class T> inline void shasta::MemoryMapped::Vector<T>::createNew(
    const string& name,
    size_t pageSize,
    size_t n,
    size_t requiredCapacity)
{
    SHASTA_ASSERT(pageSize==4096 || pageSize==2*1024*1024);

    if(name.empty()) {
        createNewAnonymous(pageSize, n, requiredCapacity);
        return;
    }

    try {
        // If already open, should have called close first.
        SHASTA_ASSERT(!isOpen);

        // Create the header.
        requiredCapacity = std::max(requiredCapacity, n);
        const Header headerOnStack(n, requiredCapacity, pageSize);
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
        for(size_t i=0; i<n; i++) {
            new(data+i) T();
        }

        // Indicate that the mapped vector is open with write access.
        isOpen = true;
        isOpenWithWriteAccess = true;
        fileName = name;

    } catch(std::exception& e) {
        cout << e.what() << endl;
        throw runtime_error("Error creating " + name);
    }

}



template<class T> inline void shasta::MemoryMapped::Vector<T>::createNewAnonymous(
    size_t pageSize,
    size_t n,
    size_t requiredCapacity)
{
    try {
        // If already open, should have called close first.
        SHASTA_ASSERT(!isOpen);

        // Create the header.
        requiredCapacity = std::max(requiredCapacity, n);
        const Header headerOnStack(n, requiredCapacity, pageSize);
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
            if(errno == ENOMEM) {
                throw runtime_error("Memory allocation failure "
                    "during mmap call for MemoryMapped::Vector.\n"
                    "This assembly requires more memory than available.\n"
                    "Rerun on a larger machine.");
            } else {
                throw runtime_error("Error " + boost::lexical_cast<string>(errno)
                    + " during mremap call for MemoryMapped::Vector: " + string(strerror(errno)));
            }
        }

        // Figure out where the data and the header go.
        header = static_cast<Header*>(pointer);
        data = reinterpret_cast<T*>(header+1);

        // Store the header.
        *header = headerOnStack;

        // Call the default constructor on the data.
        for(size_t i=0; i<n; i++) {
            new(data+i) T();
        }

        // Indicate that the mapped vector is open with write access.
        isOpen = true;
        isOpenWithWriteAccess = true;
        fileName = "";

    } catch(std::exception& e) {
        throw; // Nothing to contribute here.
    }
}


// Open a previously created vector with read-only or read-write access.
template<class T> inline void shasta::MemoryMapped::Vector<T>::accessExisting(const string& name, bool readWriteAccess)
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
template<class T> inline void shasta::MemoryMapped::Vector<T>::accessExistingReadOnly(const string& name)
{
    accessExisting(name, false);
}
template<class T> inline
    void shasta::MemoryMapped::Vector<T>::accessExistingReadWrite(
        const string& name)
{
    accessExisting(name, true);
}

// Attempt to open a previously created vector with read write access.
// If this fails, create the vector, empty, and access it with read write access.
template<class T> inline void shasta::MemoryMapped::Vector<T>::accessExistingReadWriteOrCreateNew(
    const string& name,
    size_t pageSize)
{
    try {
        accessExistingReadWrite(name);
        SHASTA_ASSERT(pageSize == header->pageSize);
    } catch(...) {
        createNew(name, pageSize);
    }

}


// Sync the mapped memory to disk.
template<class T> inline void shasta::MemoryMapped::Vector<T>::syncToDisk()
{
    SHASTA_ASSERT(isOpen);
    const int msyncReturnCode = ::msync(header, header->fileSize, MS_SYNC);
    if(msyncReturnCode == -1) {
        throw runtime_error("Error " + to_string(errno) + " during msync for " + fileName
            + ": " + ::strerror(errno) + ". Filesize is " + to_string(header->fileSize) + ".");
    }
}

// Unmap the memory.
template<class T> inline void shasta::MemoryMapped::Vector<T>::unmap()
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



template<class T> inline void shasta::MemoryMapped::Vector<T>::unmapAnonymous()
{
    SHASTA_ASSERT(isOpen);

    const int munmapReturnCode = ::munmap(header, header->fileSize);
    if(munmapReturnCode == -1) {
        throw runtime_error("Error " + boost::lexical_cast<string>(errno)
            + " unmapping MemoryMapped::Vector: " + string(strerror(errno)));
    }

    // Mark it as not open.
    isOpen = false;
    isOpenWithWriteAccess = false;
    header = 0;
    data = 0;
    fileName = "";

}


// Sync the mapped memory to disk, then unmap it.
template<class T> inline void shasta::MemoryMapped::Vector<T>::close()
{
    SHASTA_ASSERT(isOpen);
    syncToDisk();
    unmap();
}

// Close it and remove the supporting file.
template<class T> inline void shasta::MemoryMapped::Vector<T>::remove()
{
    if(fileName.empty()) {
        unmapAnonymous();
    } else {
        const string savedFileName = fileName;
        close();    // This forgets the fileName.
        filesystem::remove(savedFileName);
    }
}



// Resize works as for std::vector.
template<class T> inline void shasta::MemoryMapped::Vector<T>::resize(size_t newSize)
{
    SHASTA_ASSERT(isOpenWithWriteAccess);

    if(fileName.empty()) {
        resizeAnonymous(newSize);
        return;
    }

    const size_t oldSize = size();
    if(newSize == oldSize) {

        // No change in length - nothing to do.

    }
    if(newSize < oldSize) {

        // The vector is shrinking.
        // Just call the destructor on the elements that go away.
        for(size_t i=newSize; i<oldSize; i++) {
            (data+i)->~T();
        }
        header->objectCount = newSize;

    } else {

        // The vector is getting longer.
        if(newSize <= capacity()) {

            // No reallocation needed.
            header->objectCount = newSize;

            // Call the constructor on the elements we added.
            for(size_t i=oldSize; i<newSize; i++) {
                new(data+i) T();
            }


        } else {

            // The vector is growing beyond the current capacity.
            // We need to resize the mapped file.
            // Note that we don't have to copy the existing vector elements.

            // Save the page size.
            const size_t pageSize = header->pageSize;

            // Save the file name and close it.
            const string name = fileName;
            close();

            // Create a header corresponding to increased capacity.
            const Header headerOnStack(newSize, size_t(1.5*double(newSize)), pageSize);


            // Resize the file as necessary.
            const int fileDescriptor = openExisting(name, true);
            truncate(fileDescriptor, headerOnStack.fileSize);

            // Remap it.
            void* pointer = 0;
            try {
                pointer = map(fileDescriptor, headerOnStack.fileSize, true);
            } catch(const runtime_error& e) {
                throw runtime_error("An error occurred while resizing MemoryMapped::Vector "
                    + name + ":\n" +
                    e.what());
            }

            ::close(fileDescriptor);

            // Figure out where the data and the header are.
            header = static_cast<Header*>(pointer);
            data = reinterpret_cast<T*>(header+1);

            // Store the header.
            *header = headerOnStack;

            // Indicate that the mapped vector is open with write access.
            isOpen = true;
            isOpenWithWriteAccess = true;
            fileName = name;

            // Call the constructor on the elements we added.
            for(size_t i=oldSize; i<newSize; i++) {
                new(data+i) T();
            }
        }
    }

}



// Resize works as for std::vector.
template<class T> inline void
    shasta::MemoryMapped::Vector<T>::resizeAnonymous(size_t newSize)
{
    const size_t oldSize = size();
    if(newSize == oldSize) {

        // No change in length - nothing to do.

    }
    if(newSize < oldSize) {

        // The vector is shrinking.
        // Just call the destructor on the elements that go away.
        for(size_t i=newSize; i<oldSize; i++) {
            (data+i)->~T();
        }
        header->objectCount = newSize;

    } else {

        // The vector is getting longer.
        if(newSize <= capacity()) {

            // No reallocation needed.
            header->objectCount = newSize;

            // Call the constructor on the elements we added.
            for(size_t i=oldSize; i<newSize; i++) {
                new(data+i) T();
            }


        } else {

            // The vector is growing beyond the current capacity.
            // We need to resize the mapped file.
            // Note that we don't have to copy the existing vector elements.

            // Save the page size.
            const size_t pageSize = header->pageSize;

            // Create a header corresponding to increased capacity.
            const Header headerOnStack(newSize, size_t(1.5*double(newSize)), pageSize);



            // Remap it.
            // We can only use remap for Linux, and for 4K pages.
            bool useMremap = false;
            void* pointer = 0;
#ifdef __linux__
            useMremap = (pageSize == 4096);
#endif
            if(useMremap) {
#ifdef __linux__
                pointer = ::mremap(header, header->fileSize, headerOnStack.fileSize, MREMAP_MAYMOVE);
                if(pointer == reinterpret_cast<void*>(-1LL)) {
                    if(errno == ENOMEM) {
                        throw runtime_error("Memory allocation failure "
                            " during mremap call for MemoryMapped::Vector.\n"
                            "This assembly requires more memory than available.\n"
                            "Rerun on a larger machine.");
                    } else {
                        throw runtime_error("Error " + boost::lexical_cast<string>(errno)
                            + " during mremap call for MemoryMapped::Vector: " + string(strerror(errno)));
                    }
                }
#endif
            } else {

                // We cannot use mremap. We have to create a new mapping
                // and copy the data.
                int flags = MAP_PRIVATE | MAP_ANONYMOUS;
#ifdef __linux__
                if(pageSize == 2*1024*1024) {
                    flags |= MAP_HUGETLB | MAP_HUGE_2MB;
                }
#endif
                void* newPointer = ::mmap(0, headerOnStack.fileSize,
                    PROT_READ | PROT_WRITE, flags,
                    -1, 0);
                if(newPointer == reinterpret_cast<void*>(-1LL)) {
                    if(errno == ENOMEM) {
                        throw runtime_error("Memory allocation failure "
                            "during mmap call for MemoryMapped::Vector.\n"
                            "This assembly requires more memory than available.\n"
                            "Rerun on a larger machine.");
                    } else {
                        throw runtime_error("Error " + boost::lexical_cast<string>(errno)
                            + " during mremap call for MemoryMapped::Vector: " + string(strerror(errno)));
                    }
                }
                std::copy(
                    reinterpret_cast<char*>(header),
                    reinterpret_cast<char*>(header) + header->fileSize,
                    static_cast<char*>(newPointer));
                ::munmap(header, header->fileSize);
                pointer = newPointer;
            }



            // Figure out where the data and the header are.
            header = static_cast<Header*>(pointer);
            data = reinterpret_cast<T*>(header+1);

            // Store the header.
            *header = headerOnStack;

            // Indicate that the mapped vector is open with write access.
            isOpen = true;
            isOpenWithWriteAccess = true;
            fileName = "";

            // Call the constructor on the elements we added.
            for(size_t i=oldSize; i<newSize; i++) {
                new(data+i) T();
            }
        }
    }
}



template<class T> inline void shasta::MemoryMapped::Vector<T>::reserve()
{
    SHASTA_ASSERT(isOpenWithWriteAccess);
    reserve(size());
}


template<class T> inline void shasta::MemoryMapped::Vector<T>::reserve(size_t capacity)
{
    SHASTA_ASSERT(isOpenWithWriteAccess);
    SHASTA_ASSERT(capacity >= size());
    if(capacity == header->capacity) {
        return;
    }

    if(fileName.empty()) {
        reserveAnonymous(capacity);
        return;
    }

    // Save what we need and close it.
    const size_t currentSize = size();
    const string name = fileName;
    const size_t pageSize = header->pageSize;
    close();

    // Create a header corresponding to increased capacity.
    const Header headerOnStack(currentSize, capacity, pageSize);

    // Resize the file as necessary.
    const int fileDescriptor = openExisting(name, true);
    truncate(fileDescriptor, headerOnStack.fileSize);

    // Remap it.
    void* pointer = map(fileDescriptor, headerOnStack.fileSize, true);
    ::close(fileDescriptor);

    // Figure out where the data and the header are.
    header = static_cast<Header*>(pointer);
    data = reinterpret_cast<T*>(header+1);

    // Store the header.
    *header = headerOnStack;

    // Indicate that the mapped vector is open with write access.
    isOpen = true;
    isOpenWithWriteAccess = true;
    fileName = name;
}



template<class T> inline
    void shasta::MemoryMapped::Vector<T>::reserveAnonymous(
    size_t capacity)
{

    // Save what we need and close it.
    const size_t currentSize = size();
    const string name = fileName;
    const size_t pageSize = header->pageSize;

    // Create a header corresponding to increased capacity.
    const Header headerOnStack(currentSize, capacity, pageSize);


    // Remap it.
    // We can only use remap for Linux, and for 4K pages.
    bool useMremap = false;
#ifdef __linux__
    useMremap = (pageSize == 4096);
#endif
    void* pointer = 0;
    if(useMremap) {
#ifdef __linux__
        pointer = ::mremap(header, header->fileSize, headerOnStack.fileSize, MREMAP_MAYMOVE);
        if(pointer == reinterpret_cast<void*>(-1LL)) {
            if(errno == ENOMEM) {
                throw runtime_error("Memory allocation failure "
                    " during mremap call for MemoryMapped::Vector.\n"
                    "This assembly requires more memory than available.\n"
                    "Rerun on a larger machine.");
            } else {
                throw runtime_error("Error " + boost::lexical_cast<string>(errno)
                    + " during mremap call for MemoryMapped::Vector: " + string(strerror(errno)));
            }
        }
#endif
    } else {

        // We cannot use mremap. We have to create a new mapping
        // and copy the data.
        int flags = MAP_PRIVATE | MAP_ANONYMOUS;
#ifdef __linux__
        if(pageSize == 2*1024*1024) {
            flags |= MAP_HUGETLB | MAP_HUGE_2MB;
        }
#endif
        void* newPointer = ::mmap(0, headerOnStack.fileSize,
            PROT_READ | PROT_WRITE, flags,
            -1, 0);
        if(newPointer == reinterpret_cast<void*>(-1LL)) {
            if(errno == ENOMEM) {
                throw runtime_error("Memory allocation failure "
                    "during mmap call for MemoryMapped::Vector.\n"
                    "This assembly requires more memory than available.\n"
                    "Rerun on a larger machine.");
            } else {
                throw runtime_error("Error " + boost::lexical_cast<string>(errno)
                    + " during mremap call for MemoryMapped::Vector: " + string(strerror(errno)));
            }
        }
        std::copy(
            reinterpret_cast<char*>(header),
            reinterpret_cast<char*>(header) + header->fileSize,
            static_cast<char*>(newPointer));
        ::munmap(header, header->fileSize);
        pointer = newPointer;
    }



    // Figure out where the data and the header are.
    header = static_cast<Header*>(pointer);
    data = reinterpret_cast<T*>(header+1);

    // Store the header.
    *header = headerOnStack;

    // Indicate that the mapped vector is open with write access.
    isOpen = true;
    isOpenWithWriteAccess = true;
    fileName = "";
}



template<class T> inline void shasta::MemoryMapped::Vector<T>::unreserve()
{
    reserve(size());
}



// Make a copy of the Vector.
template<class T> inline void shasta::MemoryMapped::Vector<T>::makeCopy(
    Vector<T>& copy, const string& newName) const
    {
    copy.createNew(newName, size());
    std::copy(begin(), end(), copy.begin());
}

// Return a hash function of the stored data.
// Can be used to check for integrity.
template<class T> inline uint64_t shasta::MemoryMapped::Vector<T>::hash() const
{
    // The second argument to MurmurHash64A is a 4-byte integer.
    // Sooner or later we will have to deal with this.
    // For now we just check that there is no overflow.
    const uint64_t byteCount = size()*sizeof(T);
    SHASTA_ASSERT(byteCount <= uint64_t(std::numeric_limits<int>::max()));
    return MurmurHash64A(begin(), int(byteCount), 231);
}

template<class T> inline void shasta::MemoryMapped::Vector<T>::rename(const string& newFileName)
{
    SHASTA_ASSERT(isOpen);

    if(fileName.empty()) {
        SHASTA_ASSERT(newFileName.empty());
    } else {
        const string oldFileName = fileName;
        const bool writeAccess = isOpenWithWriteAccess;
        close();
        filesystem::move(oldFileName, newFileName);
        accessExisting(newFileName, writeAccess);
    }
}


#endif
