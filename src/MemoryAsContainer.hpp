// Simple class to make a range objects in memory appear as a container.

#ifndef SHASTA_MEMORY_AS_CONTAINER_HPP
#define SHASTA_MEMORY_AS_CONTAINER_HPP

#include "cstddef.hpp"

namespace ChanZuckerberg {
    namespace shasta {
    template<class T> class MemoryAsContainer;
    }
}



template<class T> class ChanZuckerberg::shasta::MemoryAsContainer {
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
private:
    T* dataBegin;
    T* dataEnd;
};



#endif
