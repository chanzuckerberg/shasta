#ifndef SHASTA_PEAK_FINDER_HPP
#define SHASTA_PEAK_FINDER_HPP


/***********************************************************************************************************************

 PeakFinder uses topographic prominence define peaks as sets of adjacent points surrounding a local maximum. During
 iteration, PeakFinder generates a hierarchy of peaks, for which each peak is a subpeak of the nearest taller peak.
 The global maximum becomes the parent peak of all local maximums.


 Peak domains:
 1  2  3
 ----------
 O        - ########################################  *
 X        - #####################
 X        - ###########
 X        - #######
 X        - #####
 X  X     - ####
 X  X  X  - ######
 X  X  O  - ########  *
 X  X  X  - ######
 X  X     - ########
 X  X     - ##########
 X  X     - ############
 X  O     - #############  *
 X  X     - #############
 X  X     - ###########
 X  X     - #########
 X  X     - #####
 X        - ##
 X        - #
 X        - #


 In shasta, the marker coverage distribution has a guaranteed global maximum at 1, and any number of smaller maxima
 afterwards. The 1 max is the collection of markers that contain errors, for which no or few adjacent reads had an
 identical marker. The second most prominent peak is expected to be the true distribution of coverage inside the
 marker graph, usually falling around 20-30 for a 60X dataset. The valley between the error peak and the true
 coverage peak is the target cutoff for this algorithm.

 Because there may be noise in the distribution, local minima in this valley must be ignored, which is the motivation
 for using topographic prominence. Note how the domain of the true coverage peak (2) in the diagram contains the
 domain of peak 3, which is noise. The left boundary of peak 2's domain is the target cutoff.

 The implementation is based on this page: https://www.sthu.org/blog/13-perstopology-peakdetection/index.html
 The complexity for annotating all points in the distribution is O(N) + the cost of sorting.

***********************************************************************************************************************/


#include "stdexcept.hpp"
#include "iostream.hpp"
#include "utility.hpp"
#include "vector.hpp"


namespace shasta {
    class PeakFinder;
    class PeakFinderException;
}


class shasta::PeakFinderException{};


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
