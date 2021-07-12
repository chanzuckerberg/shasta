// Shasta.
#include "InducedAlignment.hpp"
#include "PngImage.hpp"
using namespace shasta;

// Standard library.
#include "algorithm.hpp"



void InducedAlignment::writePngImage(
    uint32_t markerCount0,
    uint32_t markerCount1,
    bool useCompressedOrdinals,
    const string& fileName) const
{
    // Create the image, which gets initialized to black.
    const int n0 = useCompressedOrdinals ? compressedMarkerCount[0]: int(markerCount0);
    const int n1 = useCompressedOrdinals ? compressedMarkerCount[1]: int(markerCount1);
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
        if(useCompressedOrdinals) {
            image.setPixel(d.compressedOrdinal0, d.compressedOrdinal1, 0, 255, 0);
        } else {
            image.setPixel(d.ordinal0, d.ordinal1, 0, 255, 0);
        }
    }

    // Write it out.
    image.write(fileName);
}




