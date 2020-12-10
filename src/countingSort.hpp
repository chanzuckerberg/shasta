#ifndef SHASTA_COUNTING_SORT_HPP
#define SHASTA_COUNTING_SORT_HPP

// Counting sort of objects with integer key
// (sorting by key only - value is ignored).
// This can be faster than other sorting methods when
// the maximum key value is less than the number of
// objects to be sorted, but requires extra work areas.

#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    template<class K, class V, class Allocator> void countingSort(
        vector<pair<K, V>, Allocator>& v,
        vector<K>& count,
        vector<pair<K, V>, Allocator >& w);
}


template<class K, class V, class Allocator> void shasta::countingSort(
    vector<pair<K, V>, Allocator>& v,
    vector<K>& count,
    vector<pair<K, V>, Allocator>& w)
{
    // Count occurrences of keys.
    count.clear();
    for(const auto& p: v) {
        const K& k = p.first;
        if(count.size() <= k) {
            count.resize(k+1, 0);
        }
        ++count[k];
    }
    /*
    cout << "***Count before accumulating:" << endl;
    for(K k=0; k<count.size(); k++) {
        cout << k << " " << count[k] << endl;
    }
    */

    // Accumulate into count.
    K sum = 0;
    for(K& c: count) {
        sum+= c;
        c = sum;
    }
    /*
    cout << "***Count after accumulating:" << endl;
    for(K k=0; k<count.size(); k++) {
        cout << k << " " << count[k] << endl;
    }
    */

    // Gather into w.
    w.resize(v.size());
    for(const auto& p: v) {
        const K& k = p.first;
        // cout << "***F " << k << " " << count[k] << endl;
        w[--count[k]] = p;
    }
    // cout << "***D" << endl;

    // Copy back into v.
    copy(w.begin(), w.end(), v.begin());
    // cout << "***E" << endl;

}

#endif
