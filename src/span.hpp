#ifndef SHASTA_SPAN_HPP
#define SHASTA_SPAN_HPP

// A span class similar to std::span in C++20.
// We currently compile using the C++14 standard, so we cannot use std::span.

#include "algorithm.hpp"
#include "cstddef.hpp"
#include "iostream.hpp"
#include "iterator.hpp"
#include "stdexcept.hpp"
#include "string.hpp"

namespace shasta {
    template<class T> class MemoryAsContainer;

    // Output a MemoryAsContainer<char> as a string.
    inline ostream& operator<<(ostream&, const MemoryAsContainer<char>&);

    // Convert a MemoryAsContainer<char> to an integer.
    // This cannot be done using std::atol because the
    // MemoryAsContainer<char> is not null terminated.
    uint64_t atoul(const MemoryAsContainer<char>&);

    // Convert a MemoryAsContainer<char> to an std::string.
    string convertToString(const MemoryAsContainer<char>&);

    // MemoryAsContainer is similar to std::span in C++20.
    // Add this template typedef to ease the transition.
    template<class T> using span = MemoryAsContainer<T>;
}



template<class T> class shasta::MemoryAsContainer {
public:

    MemoryAsContainer(T* begin, T* end) :
        dataBegin(begin),
        dataEnd(end)
    {
    }

    MemoryAsContainer() : dataBegin(0), dataEnd(0) {}

    size_t size() const
    {
        return dataEnd - dataBegin;
    }
    bool empty() const
    {
        return dataBegin == dataEnd;
    }
    T* begin() const
    {
        return dataBegin;
    }
    T* end() const
    {
        return dataEnd;
    }
    T& operator[](size_t i) const
    {
        return dataBegin[i];
    }

    T& front()
    {
        SHASTA_ASSERT(dataBegin);
        SHASTA_ASSERT(dataEnd);
        SHASTA_ASSERT(!empty());
        return *dataBegin;
    }

    T& back()
    {
        SHASTA_ASSERT(dataBegin);
        SHASTA_ASSERT(dataEnd);
        SHASTA_ASSERT(!empty());
        return *(dataEnd - 1);
    }

    bool operator==(const MemoryAsContainer<T>& that) const
    {
        if(size() == that.size()) {
            return std::equal(begin(), end(), that.begin());
        } else {
            return false;
        }
    }
    bool operator!=(const MemoryAsContainer<T>& that) const
    {
        return not((*this) == that);
    }
private:
    T* dataBegin;
    T* dataEnd;
};



// Output a MemoryAsContainer<char> as a string.
inline std::ostream& shasta::operator<<(
    std::ostream& s,
    const shasta::MemoryAsContainer<char>&  m)
{
    copy(m.begin(), m.end(), ostream_iterator<char>(s));
    return s;
}


// Convert a MemoryAsContainer<char> to an integer.
// This cannot be done using std::atol because the
// MemoryAsContainer<char> is not null terminated.
// The string can only contain numeric characters.
inline uint64_t shasta::atoul(const MemoryAsContainer<char>& s)
{
    uint64_t n = 0;
    for(uint64_t i=0; ; i++) {
        const char c = s[i];
        if(not std::isdigit(c)) {
            throw runtime_error("Non-digit found in " + convertToString(s));
        }

        n += (c- '0');

        if(i ==  s.size()-1) {
            return n;
        }

        n *= 10;
    }
}



// Convert a MemoryAsContainer<char> to an std::string.
inline std::string shasta::convertToString(const MemoryAsContainer<char>& m)
{
    return string(m.begin(), m.size());
}



#endif
