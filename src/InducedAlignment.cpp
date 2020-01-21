// Shasta.
#include "InducedAlignment.hpp"
#include "PngImage.hpp"
using namespace shasta;

// Standard library.
#include "algorithm.hpp"



void InducedAlignment::writePngImage(
    uint32_t markerCount0,
    uint32_t markerCount1,
    const string& fileName) const
{
    // Create the image, which gets initialized to black.
    const int n0 = int(markerCount0);
    const int n1 = int(markerCount1);
    PngImage image(n0, n1);

    // Write a grid.
    vector<int> gridSpacing;
    vector< array<int, 3> > gridRgb;
    gridSpacing.push_back(   10); gridRgb.push_back({ 15,  15,  15});  // Grey
    gridSpacing.push_back(   50); gridRgb.push_back({ 30,  30,  30});  // Grey
    gridSpacing.push_back(  100); gridRgb.push_back({ 90,  90,  90});  // Grey
    gridSpacing.push_back(  500); gridRgb.push_back({160, 160, 160});  // Grey
    gridSpacing.push_back( 1000); gridRgb.push_back({255, 255, 255});  // White
    gridSpacing.push_back( 5000); gridRgb.push_back({255, 120, 255});  // Purple
    gridSpacing.push_back(10000); gridRgb.push_back({255, 255,  60});  // Yellow
    gridSpacing.push_back(50000); gridRgb.push_back({255, 255, 120});  // Yellow
    for(size_t i=0; i<gridSpacing.size(); i++) {
        const int spacing = gridSpacing[i];
        const array<int, 3>& rgb = gridRgb[i];
        for(int i0=0; i0<n0; i0+=spacing) {
            for(int i1=0; i1<n1; i1++) {
                image.setPixel(i0, i1, rgb[0], rgb[1], rgb[2]);
            }
        }
        for(int i1=0; i1<n1; i1+=spacing) {
            for(int i0=0; i0<n0; i0++) {
                image.setPixel(i0, i1, rgb[0], rgb[1], rgb[2]);
            }
        }
    }

    // Write the induced alignment.
    for(const InducedAlignmentData& d: data) {
        image.setPixel(d.ordinal0, d.ordinal1, 0, 255, 0);
    }

    // Write it out.
    image.write(fileName);
}



// Evaluate the quality of an induced alignment.
// Returns true if the induced alignment satisfies the specified criteria.
bool InducedAlignment::evaluate(
    uint32_t markerCount0,
    uint32_t markerCount1,
    const InducedAlignmentCriteria& inducedAlignmentCriteria) const
{
    SHASTA_ASSERT(not data.empty());

    // Compute the average and standard deviation of the offset
    // between the first and second ordinal.
    int64_t sum1 = 0;
    int64_t sum2 = 0;
    for(const auto& d: data) {
        const int64_t offset = int64_t(d.ordinal0) - int64_t(d.ordinal1);
        sum1 += offset;
        sum2 += offset * offset;
    }
    const double n = double(data.size());
    const double offset = double(sum1) / n;
    const double sigma = (data.size()==1) ? 0. : sqrt((double(sum2) - n*offset*offset) / (n-1.));
    // cout << "Offset: average " << offset << ", sigma " << sigma << endl;
    if(uint32_t(sigma) > inducedAlignmentCriteria.maxOffsetSigma) {
        /*
        cout << "Offset sigma is too large. Induced alignment follows" << endl;
        for(const auto& d: data) {
            cout << d.vertexId << " " << d.ordinal0 << " " << d.ordinal1 << " " << int(d.ordinal0)-int(d.ordinal1) << endl;
        }
        cout << "Offset sigma is too large." << endl;
        */
        return false;
    }

    // Compute ordinal sums and sort them.
    vector<uint32_t> ordinalSum;
    ordinalSum.reserve(data.size());
    for(const auto& d: data) {
        ordinalSum.push_back(d.ordinal0 + d.ordinal1);
    }
    std::sort(ordinalSum.begin(), ordinalSum.end());

    // Compute minimum and maximum ordinal sum for a perfect alignment
    // with this offset.
    const double minOrdinalSum = abs(offset);
    const double maxOrdinalSum = min(
        double(2*markerCount1) + offset,
        double(2*markerCount0) - offset);
    /*
    cout << "Minimum ordinal sum: ideal " << minOrdinalSum <<
        ", actual " << ordinalSum.front() <<
        ", deviation " << double(ordinalSum.front()) - minOrdinalSum << endl;
    cout << "Maximum ordinal sum: ideal " << maxOrdinalSum <<
        ", actual " << ordinalSum.back() <<
        ", deviation " << double(ordinalSum.back()) - maxOrdinalSum << endl;
    */

    if(abs(double(ordinalSum.front()) - minOrdinalSum) >
        double(2*inducedAlignmentCriteria.maxTrim)) {
        // cout << "Too much trim on left." << endl;
        return false;
    }
    if(abs(double(ordinalSum.back()) - maxOrdinalSum) >
        double(2*inducedAlignmentCriteria.maxTrim)) {
        // cout << "Too much trim on right." << endl;
        return false;
    }

    // Check for gaps.
    for(uint64_t i=1; i<ordinalSum.size(); i++) {
        if(ordinalSum[i] - ordinalSum[i-1] > 2*inducedAlignmentCriteria.maxSkip) {
            // cout << "Large gap." << endl;
            return false;
        }
    }

    // If getting here, this is a good induced alignment according to
    // the given criteria.
    // cout << "Good induced alignment." << endl;
    return true;
}

