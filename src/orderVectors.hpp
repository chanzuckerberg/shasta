#ifndef SHASTA_ORDER_VECTORS_HPP
#define SHASTA_ORDER_VECTORS_HPP

#include "vector.hpp"

// Classes to sort vectors by size.

namespace shasta {

    template<class T> class OrderVectorsByIncreasingSize;
    template<class T> class OrderVectorsByDecreasingSize;

}



template<class T> class shasta::OrderVectorsByIncreasingSize {
public:
     bool operator()(const vector<T>& x, const vector<T>& y) const
    {
         return x.size() < y.size();
    }
};



template<class T> class shasta::OrderVectorsByDecreasingSize {
public:
     bool operator()(const vector<T>& x, const vector<T>& y) const
    {
         return x.size() > y.size();
    }
};

#endif

