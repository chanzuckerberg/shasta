#ifndef SHASTA_PEAK_FINDER_HPP
#define SHASTA_PEAK_FINDER_HPP


#include "stdexcept.hpp"
#include "iostream.hpp"
#include "utility.hpp"
#include "vector.hpp"

using std::runtime_error;
using std::ostream;
using std::pair;
using std::vector;
using std::string;
using std::cout;

namespace shasta {
    class PeakFinder;
}


class shasta::PeakFinder {
private:
    class Peak {
    public:
        /// Attributes ///

        // Start is the x value corresponding to the highest point in the peak
        uint64_t start;

        // Stop is the x value corresponding to lowest point in the peak before it was merged.
        uint64_t stop;

        // During iteration, left and right are the x bounds of each growing peak
        uint64_t left;
        uint64_t right;

        // Has this peak been taken over by a taller adjacent peak?
        bool isMerged;

        // Persistence is the height from the top of the peak to its base, where it was merged with another peak
        uint64_t persistence;

        /// Methods ///
        Peak(uint64_t start);
    };


public:
    /// Attributes ///

    // This vector will be filled with peaks during iteration
    vector<Peak> peaks;


    /// Methods ///

    // Iterate the y values and identify all the peaks/valleys
    void findPeaks(const vector<uint64_t>& y);

    void sortByPersistence();
    uint64_t findXCutoff(const vector<uint64_t>& y);
    uint64_t calculateArea(const vector<uint64_t>& y, uint64_t xMin, uint64_t xMax);

};



#endif //SHASTA_PEAK_FINDER_HPP
