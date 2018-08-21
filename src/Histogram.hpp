#ifndef CZI_SHASTA_HISTOGRAM_HPP
#define CZI_SHASTA_HISTOGRAM_HPP

#include "algorithm.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        class Histogram;
    }
}



// A simple histogram class that stores the histogram in a vector.
class ChanZuckerberg::shasta::Histogram : public vector<size_t> {
public:

    void increment(size_t i)
    {
        const size_t requiredSize = i + 1ULL;
        if(size() < requiredSize) {
            resize(requiredSize, 0ULL);
        }
        ++(*this)[i];
    }

    size_t sum() const
    {
        return std::accumulate(begin(), end(), 0ULL);
    }

    // Return the histogram position with the maximum frequency.
    size_t bestValue() const
    {
        return std::max_element(begin(), end()) - begin();
    }

    // Return the best frequency stored in the histogram.
    size_t bestFrequency() const
    {
        return *std::max_element(begin(), end());
    }
};


#endif
