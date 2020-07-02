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
    class Peak;
    class PeakFinder;
}

class shasta::Peak {
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


ostream& operator<<(ostream& o, shasta::Peak& peak){
    o << peak.start << ',' << peak.stop << ',' << peak.isMerged << ',' << peak.left << ',' << peak.right << ',' << peak.persistence;

    return o;
}


shasta::Peak::Peak(uint64_t start){
    this->left = start;
    this->right = start;
    this->start = start;
    this->stop = 0;
    this->isMerged = false;
}


class shasta::PeakFinder {
public:
    /// Attributes ///

    // This vector will be filled with peaks during iteration
    vector<Peak> peaks;


    /// Methods ///

    // Iterate the y values and identify all the peaks/valleys
    void findPeaks(vector<uint64_t>& y);

    void sortByPersistence();
    uint64_t findXCutoff(vector<uint64_t>& y);
    uint64_t calculateArea(vector<uint64_t>& y, uint64_t x_min, uint64_t x_max);
};



#endif //SHASTA_PEAK_FINDER_HPP
