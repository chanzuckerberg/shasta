#ifndef SHASTA_PREFIX_LENGTH_HPP
#define SHASTA_PREFIX_LENGTH_HPP

#include "algorithm.hpp"

namespace shasta {

    template<class Container> uint64_t commonPrefixLength(
        const Container& x,
        const Container& y)
    {
        if(x.size() <= y.size()) {
            return std::mismatch(x.begin(), x.end(), y.begin()).first - x.begin();
        } else {
            return std::mismatch(y.begin(), y.end(), x.begin()).first - y.begin();
        }
    }



    template<class Container> uint64_t commonSuffixLength(
        const Container& x,
        const Container& y)
    {
        if(x.size() <= y.size()) {
            return std::mismatch(x.rbegin(), x.rend(), y.rbegin()).first - x.rbegin();
        } else {
            return std::mismatch(y.rbegin(), y.rend(), x.rbegin()).first - y.rbegin();
        }
    }
}



#endif
