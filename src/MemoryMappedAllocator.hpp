#ifndef SHASTA_MEMORY_MAPPED_ALLOCATOR_HPP
#define SHASTA_MEMORY_MAPPED_ALLOCATOR_HPP

// A simple allocator compatible with standard containers which
// allocates memory from a MemoryMapped::Vector.
// Memory is always allocated at the end of the MemoryMapped::Vector
// and never freed, except when the allocator is destroyed.
// This could help performance in some situations,
// if the MemoryMapped::Vector is on 2 MB pages.

#include "MemoryMappedVector.hpp"

namespace shasta {
    namespace MemoryMapped {

        // The allocator class compatible with standard containers.
        template<class T> class Allocator;

        // Low level class used by above class Allocator.
        class ByteAllocator;

        class BadAllocation {};

        // Test function.
        void testMemoryMappedAllocator();
    }
}




class shasta::MemoryMapped::ByteAllocator {
public:

    // Create a ByteAllocator with this number of bytes,
    // which will never be increased.
    // The allocator will throw BadAllocation if
    // this space is insufficient.
    ByteAllocator(const string& name, uint64_t pageSize, uint64_t n) :
        allocatedByteCount(0),
        allocatedBlockCount(0)
    {
        data.createNew(name, pageSize);
        data.resize(n); // Never resize after this!
    }

    ~ByteAllocator()
    {
        SHASTA_ASSERT(allocatedBlockCount == 0);
        data.remove();
    }



    char* allocate(uint64_t n, size_t objectSize)
    {
        // Figure out the number of bytes.
        uint64_t byteCount = n * objectSize;

        // Make sure to always allocate a multiple of 8 bytes.
        const uint64_t remainder = byteCount & 7;
        if(remainder) {
            byteCount += (8-remainder);
        }
        SHASTA_ASSERT((byteCount & 7) == 0);

        // Figure out the new number of allocated bytes.
        const uint64_t newAllocatedByteCount = allocatedByteCount + byteCount;

        // If not enough space, throw.
        if(newAllocatedByteCount > data.size()) {
            throw BadAllocation();
        }

        // All good.
        char* p = data.begin() + allocatedByteCount;
        allocatedByteCount = newAllocatedByteCount;
        ++allocatedBlockCount;
        /* cout << "Requested " << n * objectSize << ", allocated " << byteCount <<
            ", total allocated " << allocatedByteCount << endl; */
        return p;
    }


    // Memory only gets deallocated when the allocator is destroyed.
    void deallocate()
    {
        // cout << "Deallocate" << endl;
        --allocatedBlockCount;
    }

private:
    Vector<char> data;
    uint64_t allocatedByteCount;
    uint64_t allocatedBlockCount;
};



template<class T> class shasta::MemoryMapped::Allocator {
public:

    // Construct an Allocator<T> that will use a given
    // ByteAllocator.
    Allocator<T>(ByteAllocator& byteAllocator) :
        byteAllocator(byteAllocator) {}

    // The "extended" copy constructor is required.
    // It creates an Allocator<T> that uses the same ByteAllocator.
    template<typename U> friend class Allocator;
    template<class U> Allocator<T>(const Allocator<U>& that) :
        byteAllocator(that.byteAllocator) {}

    // Required types.
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;

    // Allocate a block of n objects of type T.
    pointer allocate (uint64_t n)
    {
        return reinterpret_cast<pointer>(byteAllocator.allocate(n, sizeof(T)));
    }

    // Deallocate a block.
    // The memory does not actually get freed until the ByteAllocator
    // is destroyed.
    void deallocate(pointer p, uint64_t n)
    {
        byteAllocator.deallocate();
    }

    // Other required functions.
    pointer address(reference t)
    {
        return &t;
    }
    const_pointer address(const_reference t)
    {
        return &t;
    }
    template <class U, class... Args> void construct(U* p, Args&&... args)
    {
        ::new((void*)p) U (std::forward<Args>(args)...);
    }
    template <class U> void destroy(U* p)
    {
        p->~U();
    }
    template <class U> struct rebind {
        typedef Allocator<U> other;
    };



private:
    ByteAllocator& byteAllocator;
};

#endif
