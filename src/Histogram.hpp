#ifndef SHASTA_HISTOGRAM_HPP
#define SHASTA_HISTOGRAM_HPP

#include "algorithm.hpp"
#include "vector.hpp"
#include <numeric>
#include <cstdint>
#include <fstream>

namespace shasta {
    class Histogram;
    class IterativeHistogram;
    void testIterativeHistogram();
}

using std::ostream;


// A simple histogram class that stores the histogram in a vector.
class shasta::Histogram : public vector<size_t> {
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


// A histogram class that calculates bin size and finds the appropriate bin to increment given an x value
class shasta::IterativeHistogram{
public:
    /// Attributes ///
    const double start;
    const double stop;
    const size_t nBins;
    const double binSize;
    vector<uint64_t> histogram;
    bool unboundedLeft;
    bool unboundedRight;

    /// Methods ///
    IterativeHistogram(
            double start,
            double stop,
            size_t nBins,
            bool unboundedLeft=false,
            bool unboundedRight=false);

    void update(double x);
    void getNormalizedHistogram(vector<double>& normalizedHistogram);
    void writeToHtml(ostream& html, uint64_t sizePx);

private:
    /// Methods ///
    int64_t findIndex(double x);
};


#endif
