#ifndef SHASTA_PREFIX_LENGTH_HPP
#define SHASTA_PREFIX_LENGTH_HPP

#include "algorithm.hpp"
#include "vector.hpp"

namespace shasta {

    template<class T> uint64_t commonPrefixLength(
        const vector<T>& x,
        const vector<T>& y)
    {
        typename vector<T>::const_iterator it;

        if(x.size() <= y.size()) {
            tie(it, ignore) = mismatch(x.begin(), x.end(), y.begin());
            return it - x.begin();
        } else {
            tie(it, ignore) = mismatch(y.begin(), y.end(), x.begin());
            return it - y.begin();
        }
    }



    template<class T> uint64_t commonSuffixLength(
        const vector<T>& x,
        const vector<T>& y)
    {
        typename vector<T>::const_reverse_iterator it;

        if(x.size() <= y.size()) {
            tie(it, ignore) = mismatch(x.rbegin(), x.rend(), y.rbegin());
            return it - x.rbegin();
        } else {
            tie(it, ignore) = mismatch(y.rbegin(), y.rend(), x.rbegin());
            return it - y.rbegin();
        }
    }
}



#endif
