#ifndef SHASTA_SPLIT_RANGE_HPP
#define SHASTA_SPLIT_RANGE_HPP

// Given an half-open range [begin, end), of length n=end-begin,
// divide it in m slices that are as uniform in length as
// possible and return [sliceBegin, sliceEnd) for the
// i-th slice, where i is in [0, m).

// Here, [x, y) indicate the half-open range that
// begins at x (included) and ends at y (excluded).

// We write n = a*m + b, with a=n/m and b=n%m.
// The first b intervals are of size a+1
// and the remaining (m-b) intervals are of size a.

// Check that the total size of all ranges is correct:
// b*(a+1) + (m-b)*a = b*a + b + m*a -b*a = m*a + b = n.

#include "cstddef.hpp"
#include "iostream.hpp"
#include "utility.hpp"

namespace shasta {

    inline pair<size_t, size_t> splitRange(
        size_t begin,
        size_t end,
        size_t m,
        size_t i
        );

    inline void testSplitRange();

}



inline std::pair<size_t, size_t> shasta::splitRange(
    size_t begin,
    size_t end,
    size_t m,
    size_t i
    )
{
    SHASTA_ASSERT(m > 0);
    const size_t n = end - begin;
    const size_t a = n / m;
    const size_t b = n % m;

    pair<size_t, size_t> p;
    if(i < b) {
        p.first = begin + i * (a + 1);
        p.second = p.first + (a + 1);
    } else {
        p.first = begin + b * (a+1) + (i-b)*a;
        p.second = p.first + a;
    }
    return p;
}



inline void shasta::testSplitRange()
{
    while(true) {
        cout << "Enter begin, end, m:" << endl;
        size_t begin;
        size_t end;
        size_t m;
        cin >> begin >> end >> m;
        for(size_t i=0; i<m; i++) {
            const auto p = splitRange(begin, end, m, i);
            cout << i << ": " << p.first << " " << p.second << endl;
        }
    }

}



#endif
