#include "PeakFinder.hpp"
#include <algorithm>

using namespace shasta;
using std::sort;
using std::cerr;


PeakFinder::Peak::Peak(uint64_t start):
        start(start),
        stop(0),
        left(start),
        right(start),
        isMerged(false),
        persistence(0)
{}


void PeakFinder::findPeaks(const vector<uint64_t>& y){
    size_t xSize = y.size();

    // Create a vector of values, in which each value can be used to fetch the parent peak for that position in the
    // y. If that position does not yet belong to a peak, it has index -1
    vector<int64_t> peakIndex(xSize, -1);

    // Create a vector of x values which point to y values
    vector<size_t> indexes(xSize, 0);
    for (size_t i=0; i<indexes.size() ; i++) {
        indexes[i] = i;
    }

    // Sort the x values based on their corresponding y
    sort(indexes.begin(), indexes.end(),
         [&](const size_t& a, const size_t& b) {
             // Equal y values are sorted by lowest x value
             if (y[a] == y[b]){
                 return (a < b);
             }
             else {
                 return (y[a] > y[b]);
             }
         }
    );

    // Merge peaks
    for (auto& i: indexes){
        bool hasLeftPeak = false;
        if (i > 0){
            hasLeftPeak = (peakIndex[i - 1] >= 0);
        }

        bool hasRightPeak = false;
        if (i < xSize - 1){
            hasRightPeak = (peakIndex[i + 1] >= 0);
        }

        // No adjacent peaks
        if (not hasLeftPeak and not hasRightPeak){
            Peak p(i);
            peaks.push_back(p);
            peakIndex[i] = int64_t(peaks.size() - 1);
        }
        // Left side is peak
        else if (hasLeftPeak and not hasRightPeak){
            Peak& leftPeak = peaks[size_t(peakIndex[i - 1])];
            leftPeak.right = i;
            peakIndex[i] = peakIndex[i - 1];
        }
        // Right side is peak
        else if (not hasLeftPeak and hasRightPeak){
            Peak& rightPeak = peaks[size_t(peakIndex[i + 1])];
            rightPeak.left = i;
            peakIndex[i] = peakIndex[i + 1];
        }
        // Both sides are peaks
        else{
            Peak& leftPeak = peaks[size_t(peakIndex[i - 1])];
            uint64_t xLeft = leftPeak.start;
            uint64_t yLeft = y[xLeft];

            Peak& rightPeak = peaks[size_t(peakIndex[i + 1])];
            uint64_t xRight = rightPeak.start;
            uint64_t yRight = y[xRight];

            // Right peak is bigger
            if (yRight > yLeft){
                rightPeak.left = leftPeak.left;
                peakIndex[i] = peakIndex[i + 1];

                // Update the bounds of the dead peak so it stops (inclusive) at the merge point
                leftPeak.right = i;

                // Update boundary peak index
                peakIndex[leftPeak.left] = peakIndex[i + 1];
                peakIndex[leftPeak.right] = peakIndex[i + 1];

                // Update the weaker peak to reflect that is has terminated
                leftPeak.stop = i;
                leftPeak.isMerged = true;
                leftPeak.persistence = y[rightPeak.start] - y[i];
            }
            // Left peak is bigger, or same size
            else{
                leftPeak.right = rightPeak.right;
                peakIndex[i] = peakIndex[i - 1];

                // Update the bounds of the dead peak so it stops (inclusive) at the merge point
                rightPeak.left = i;

                // Update boundary peak index
                peakIndex[rightPeak.right] = peakIndex[i - 1];
                peakIndex[rightPeak.left] = peakIndex[i - 1];

                // Update the weaker peak to reflect that is has terminated
                rightPeak.stop = i;
                rightPeak.isMerged = true;
                rightPeak.persistence = y[rightPeak.start] - y[i];
            }
        }
    }

    // Set the persistence for the tallest peak as its total height from 0
    peaks[0].persistence = y[peaks[0].start];
}


void PeakFinder::sortByPersistence(){
    // Sort the peaks based on their persistence
    sort(peaks.begin(), peaks.end(),
         [&](const Peak& a, const Peak& b) {
             // Equal values are sorted by lowest x value
             if (a.persistence == b.persistence){
                 return (a.start < b.start);
             }
             else {
                 return (a.persistence > b.persistence);
             }
         }
    );
}


uint64_t PeakFinder::calculateArea(const vector<uint64_t>& y, uint64_t xMin, uint64_t xMax){
    uint64_t total = 0;

    for (size_t i=xMin; i <= xMax; i++){
        total += y[i];
    }

    return total;
}


uint64_t PeakFinder::findXCutoff(const vector<uint64_t>& y){
    sortByPersistence();

    // find the total AUC
    uint64_t totalArea = calculateArea(y, 1, y.size() - 1);

    // Second most prominent peak
    Peak& p = peaks[1];

    // Find the AUC inside the bounds of the peak (from y=0 and up)
    uint64_t peakArea = calculateArea(y, p.left, p.right);

    double percentArea = 100 * (double(peakArea) / double(totalArea));

    uint64_t xCutoff;

    // Check if second most prominent peak is reasonable size (not a false peak)
    if (percentArea > 1){   //TODO: readjust this proportion for histograms with 0=0 frequency
        xCutoff = p.left;
    }
    else{
        throw runtime_error("ERROR: No significant cutoff found in marker frequency histogram");
    }

    return xCutoff;
}

