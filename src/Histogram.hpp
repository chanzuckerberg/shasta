#ifndef SHASTA_HISTOGRAM_HPP
#define SHASTA_HISTOGRAM_HPP

#include "algorithm.hpp"
#include <deque>
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include <numeric>
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {
    class Histogram;
    class Histogram2;
    void testIterativeHistogram();
    void writeHistogramsToHtml(
        ostream& html,
        Histogram2& histogramA,
        Histogram2& histogramB,
        uint64_t sizePx,
        int32_t precision);
}


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


// A histogram class that calculates bin size and finds the appropriate bin to increment given an x value.
// When "dynamicBounds" is set to true, no input value is considered out of bounds, because the histogram will
// dynamically resize to reach new values, in both directions on the number line. However, this means that
// "unboundedLeft" and "unboundedRight" are ignored. Otherwise, "unbounded..." params will modify the edge behavior
// of the first or last bin so that any out of range item will be funneled into the terminal bin.
// There is currently no initialization parameter to directly specify binSize, so when using dynamicBounds,
// start, stop, and binCount must be specified as a means of indirectly specifying binSize.
class shasta::Histogram2{
public:
    /// Attributes ///
    double start;
    double stop;
    uint64_t binCount;
    const double binSize;
    std::deque<uint64_t> histogram;
    bool unboundedLeft;
    bool unboundedRight;
    bool dynamicBounds;

    /// Methods ///
    Histogram2(
            double start,
            double stop,
            uint64_t binCount,
            bool unboundedLeft=false,
            bool unboundedRight=false,
            bool dynamicBounds=false);

    void update(double x);
    void getNormalizedHistogram(vector<double>& normalizedHistogram);
    void writeToHtml(ostream& html, uint64_t sizePx, int32_t precision);
    void writeToCsv(ostream& csv, int32_t precision);
    double thresholdByCumulativeProportion(double fraction);
    pair<string,string> getBoundStrings(uint64_t binIndex, int32_t precision);
    uint64_t getSum();
    int64_t findIndex(double x);
};


#endif
