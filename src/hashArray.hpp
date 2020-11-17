#ifndef SHASTA_HASH_ARRAY_HPP
#define SHASTA_HASH_ARRAY_HPP

#include "MurmurHash2.hpp"
#include "array.hpp"

namespace shasta {
    template<class T> class HashTuple;
}

// T can be an array, pair, or tuple of types not containing pointers or padding.
template<class T> class shasta::HashTuple {
public:
    uint64_t operator()(const T& v) const
    {
        return MurmurHash64A(&v, sizeof(v), 15741);
    }
};


#endif
