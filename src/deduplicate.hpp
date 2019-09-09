#ifndef SHASTA_DEDUPLICATE_HPP
#define SHASTA_DEDUPLICATE_HPP

#include "algorithm.hpp"
#include "vector.hpp"

namespace shasta {

    // Remove duplicate elements in a vector.
    template<class T> void deduplicate(vector<T>& v)
    {
        sort(v.begin(), v.end());
        v.resize(unique(v.begin(), v.end()) - v.begin());
    }


    // Remove duplicate elements in a vector and count occurrences of each
    template<class T, class Int> void deduplicateAndCount(
        vector<T>& v,
        vector<Int>& count)
    {
        // Clear the count vector.
        count.clear();

        // If the given vector is empty, return now.
        if(v.empty()) {
            return;
        }

        // Sort the vector.
        sort(v.begin(), v.end());

        // Add elements, keeping track of the number
        // of occurrences of each.
        typename vector<T>::iterator output = v.begin();
        typename vector<T>::iterator input = v.begin();
        while(input != v.end()) {

            // Store this element.
            *output = *input;
            ++output;

            // Count how many there are.
            typename vector<T>::iterator it = input;
            while(it!=v.end() && *it==*input) {
                ++it;
            }

            // Store the count.
            count.push_back(Int(it - input));

            // Update our output iterator.
            input = it;

        }
        v.resize(count.size());
    }


    inline void testDeduplicateAndCount()
    {
        vector<int> v = {7, 4, 5, 7, 4, 18, 2, 4};
        vector<int> count;
        deduplicateAndCount(v, count);
        SHASTA_ASSERT(v.size() == count.size());
        for(uint64_t i=0; i<v.size(); i++) {
            cout << v[i] << " " << count[i] << endl;
        }
    }
}

#endif
