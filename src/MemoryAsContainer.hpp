// Simple class to make a range objects in memory appear as a container.

#ifndef CZI_NANOPORE2_MEMORY_AS_CONTAINER_HPP
#define CZI_NANOPORE2_MEMORY_AS_CONTAINER_HPP

#include "cstddef.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {
    template<class T> class MemoryAsContainer;
    }
}



template<class T> class ChanZuckerberg::Nanopore2::MemoryAsContainer {
public:

    MemoryAsContainer(T* begin, T* end) :
        dataBegin(begin),
        dataEnd(end)
    {
    }

    size_t size() const
    {
        return dataEnd - dataBegin;
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

private:
    T* dataBegin;
    T* dataEnd;
};



#endif
