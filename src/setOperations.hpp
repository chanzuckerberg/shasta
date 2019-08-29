#ifndef SHASTA_SET_OPERATIONS_HPP
#define SHASTA_SET_OPERATIONS_HPP


#include "cstdint.hpp"


namespace shasta {

    // Given two sorted ranges without duplicate elements
    // representing two sets, compute the size of
    // their intersection - that is, the number of common elements.
    template<class Iterator> uint64_t intersectionSize(
        Iterator begin0, Iterator end0,
        Iterator begin1, Iterator end1)
    {
        uint64_t n = 0;
        Iterator it0 = begin0;
        Iterator it1 = begin1;
        while(it0!=end0 && it1!=end1) {
            if(*it0 < *it1) {
                ++it0;
            } else if(*it1 < *it0) {
                ++it1;
            } else {
                ++n;
                ++it0;
                ++it1;
            }
        }
        return n;
    }

}



#endif
