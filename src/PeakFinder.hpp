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


ostream& operator<<(ostream& o, shasta::Peak& peak){
    o << peak.start << ',' << peak.stop << ',' << peak.isMerged << ',' << peak.left << ',' << peak.right << ',' << peak.persistence;

    return o;
}


class shasta::PeakFinder {
public:
    /// Attributes ///

    // This vector will be filled with peaks during iteration
    vector<Peak> peaks;


    /// Methods ///

    // Iterate the y values and identify all the peaks/valleys
    void findPeaks(const vector<uint64_t>& y);

    void sortByPersistence();
    uint64_t findXCutoff(const vector<uint64_t>& y);
    uint64_t calculateArea(const vector<uint64_t>& y, uint64_t x_min, uint64_t x_max);
};



#endif //SHASTA_PEAK_FINDER_HPP
