// Class to describe a vector of vectors stored contiguously in mapped memory.
// A table of contents (toc) contains indexes pointing to the first element of each vector.

#ifndef SHASTA_MEMORY_MAPPED_VECTOR_OF_VECTORS_HPP
#define SHASTA_MEMORY_MAPPED_VECTOR_OF_VECTORS_HPP

// Shasta.
#include "MemoryMappedVector.hpp"
#include "span.hpp"

// Standard libraries.
#include "algorithm.hpp"
#include "utility.hpp"
#include "vector.hpp"

// Forward declarations.
namespace shasta {
    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }
}



template<class T, class Int> class shasta::MemoryMapped::VectorOfVectors {
public:

    void createNew(const string& nameArgument, size_t pageSizeArgument)
    {
        name = nameArgument;
        pageSize = pageSizeArgument;

        if(nameArgument.empty()) {
            toc.createNew("", pageSize);
            data.createNew("", pageSize);
        } else {
            toc.createNew(name + ".toc", pageSize);
            data.createNew(name + ".data", pageSize);
        }

        toc.push_back(0);
    }



    void accessExisting(const string& nameArgument, bool readWriteAccess)
    {
        name = nameArgument;
        pageSize = 0;   // We cannot use count.
        toc.accessExisting(name + ".toc", readWriteAccess);
        data.accessExisting(name + ".data", readWriteAccess);
    }
    void accessExistingReadOnly(const string& name)
    {
        accessExisting(name, false);
    }

    void accessExistingReadWrite(const string& name)
    {
        accessExisting(name, true);
    }

    // Attempt to open a previously created name with read write access.
    // If this fails, create the vector, empty, and access it with read write access.
    void accessExistingReadWriteOrCreateNew(
        const string& name,
        size_t pageSize)
    {
        try {
            accessExistingReadWrite(name);
        } catch(...) {
            createNew(name, pageSize);
        }
    }

    void clear()
    {
        toc.clear();
        toc.push_back(0);
        data.clear();
        if(count.isOpen) {
            count.clear();
        }
    }

    void remove()
    {
        toc.remove();
        data.remove();
        if(count.isOpen) {
            count.remove();
        }
    }

    size_t size() const
    {
        return toc.size() - 1;
    }
    size_t totalSize() const
    {
        return data.size();
    }
    void close()
    {
        toc.close();
        data.close();
    }
    bool empty() const
    {
        return toc.size() == 1;
    }

    T* begin()
    {
        return data.begin();
    }
    const T* begin() const
    {
        return data.begin();
    }

    T* end()
    {
        return data.end();
    }
    const T* end() const
    {
        return data.end();
    }


    // Return size/begin/end of the i-th vector.
    size_t size(size_t i) const
    {
        return toc[i+1] - toc[i];
    }
    T* begin(Int i)
    {
        return data.begin() + toc[i];
    }
    const T* begin(Int i) const
    {
        return data.begin() + toc[i];
    }
    T* end(Int i)
    {
        return data.begin() + toc[i+1];
    }
    const T* end(Int i) const
    {
        return data.begin() + toc[i+1];
    }

   // Add an empty vector at the end.
    void appendVector()
    {
        const Int tocBack = toc.back();
        toc.push_back(tocBack);
    }

    // Add a T at the end of the last vector.
    void append(const T& t)
    {
        // SHASTA_ASSERT(!empty());
        ++toc.back();
        data.push_back(t);
    }



    // Add a non-empty vector at the end.
    template<class Iterator> void appendList(Iterator begin, Iterator end)
    {
        // First, append an empty vector.
        appendVector();

        // Then, append all the elements to the last vector.
        for(Iterator it=begin; it!=end; ++it) {
            append(*it);
        }
    }
    template<class Iterator> void appendVector(Iterator beginIt, Iterator endIt)
    {
        appendVector(endIt - beginIt);
        copy(beginIt, endIt, begin(size()-1));
    }
    void appendVector(const vector<T>& v)
    {
        // appendVector(v.begin(), v.end());
        appendVector(v.size());
        copy(v.begin(), v.end(), begin(size()-1));
    }



    // Add a non-empty vector at the end, without initializing it.
    void appendVector(Int n)
    {
        toc.push_back(toc.back() + n);
        data.resize(toc.back());
    }


    // Append vectors at the end.
    void appendVectors(const vector< vector<T> >& v) {
        for(const auto x: v) {
            appendVector(x);
        }
    }



    // Operator[] returns a span object containing all elements
    // of the vector with the requested index.
    span<T> operator[](Int i)
    {
        return span<T>(begin(i), end(i));
    }
    span<const T> operator[](Int i) const
    {
        return span<const T>(begin(i), end(i));
    }

    // Front and back, in const and non-const versions.
    span<T> front() {
        SHASTA_ASSERT(size() > 0);
        return (*this)[0];
    }
    span<const T> front() const {
        SHASTA_ASSERT(size() > 0);
        return (*this)[0];
    }
    span<T> back() {
        SHASTA_ASSERT(size() > 0);
        return (*this)[size() - 1];
    }
    span<const T> back() const {
        SHASTA_ASSERT(size() > 0);
        return (*this)[size() - 1];
    }



    // Function to construct the VectorOfVectors in two passes.
    // In pass 1 we count the number of entries in each of the vectors.
    // In pass 2 we store the entries.
    // This can be easily turned into multithreaded code
    // if atomic memory access primitives are used.
    void beginPass1(Int n);
    void incrementCount(Int index, Int m=1);  // Called during pass 1.
    void incrementCountMultithreaded(Int index, Int m=1);  // Called during pass 1.
    void beginPass2();
    void store(Int index, const T&);            // Called during pass 2.
    void storeMultithreaded(Int index, const T&);            // Called during pass 2.
    void endPass2(bool check = true, bool free=true);

    // Touch the memory in order to cause the
    // supporting pages of virtual memory to be loaded in real memory.
    size_t touchMemory() const
    {
        return toc.touchMemory() + data.touchMemory();
    }

    bool isOpen() const
    {
        return toc.isOpen && data.isOpen;
    }
    bool isOpenWithWriteAccess() const
    {
        return toc.isOpenWithWriteAccess && data.isOpenWithWriteAccess;
    }

    // Given a global index k in a VectorOfVectors v,
    // find i and j such that v[i][j] (aka v.begin(i)[j]) is the same
    // (stored at the same position) as v.begin()[k].
    // This requires a binary search in the toc.
    pair<Int, Int> find(Int k) const;

    void rename(const string& newName)
    {
        toc.rename(newName + ".toc");
        data.rename(newName + ".data");
    }

    string getName() const
    {
        return name;
    }

private:
    Vector<Int> toc;
    Vector<Int> count;
    Vector<T> data;
    string name;
    size_t pageSize;
};



template<class T, class Int>
    void shasta::MemoryMapped::VectorOfVectors<T, Int>::beginPass1(Int n)
{

    if(!count.isOpen) {
        if(name.empty()) {
            count.createNew("", pageSize);
        } else {
            count.createNew(name + ".count", pageSize);
        }
    }
    count.reserveAndResize(n);
    fill(count.begin(), count.end(), Int(0));
}



template<class T, class Int>
    void shasta::MemoryMapped::VectorOfVectors<T, Int>::beginPass2()
{
    const Int n = Int(count.size());
    toc.reserveAndResize(n+1);
    toc[0] = 0;
    for(Int i=0; i<n; i++) {
        toc[i+1] = toc[i] + count[i];
    }
    const size_t  dataSize = toc.back();
    data.reserveAndResize(dataSize);
}



template<class T, class Int>
    void shasta::MemoryMapped::VectorOfVectors<T, Int>::endPass2(
        bool check, bool free)
{
    // Verify that all counts are now zero.
    if(check) {
        const Int n = Int(count.size());
        for(Int i=0; i<n; i++) {
            SHASTA_ASSERT(count[i] == 0);;
        }
    }

    // Free the memory of the count vector.
    if(free) {
        count.remove();
    } else {
        // This leaves the memory allocated.
        // Faster if we want to reuse the VectorOfVectors.
        count.clear();
    }
}



template<class T, class Int>
    void shasta::MemoryMapped::VectorOfVectors<T, Int>::incrementCount(Int index, Int m)
{
    count[index] += m;
}
template<class T, class Int>
    void shasta::MemoryMapped::VectorOfVectors<T, Int>::
    incrementCountMultithreaded(Int index, Int m)
{
    __sync_fetch_and_add(&count[index], m);
}


template<class T, class Int>
    void shasta::MemoryMapped::VectorOfVectors<T, Int>::store(Int index, const T& t)
{
    (*this)[index][--count[index]] = t;
}
template<class T, class Int>
    void shasta::MemoryMapped::VectorOfVectors<T, Int>::storeMultithreaded(Int index, const T& t)
{
    const Int i = __sync_sub_and_fetch(&count[index], 1);
    (*this)[index][i] = t;
}



// Given a global index k in a VectorOfVectors v,
// find i and j such that v[i][j] (aka v.begin(i)[j]) is the same
// (stored at the same position) as v.begin()[k].
// This requires a binary search in the toc.
template<class T, class Int>
    std::pair<Int, Int> shasta::MemoryMapped::VectorOfVectors<T, Int>::find(Int k) const
{
    const auto it = std::upper_bound(toc.begin(), toc.end(), k) - 1;
    const Int i = it - toc.begin();
    const Int j = k - *it;
    SHASTA_ASSERT(i < size());
    SHASTA_ASSERT(j < size(i));
    return make_pair(i, j);
}

#endif
