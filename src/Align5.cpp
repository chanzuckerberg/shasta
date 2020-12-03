#include "Align5.hpp"
#include "orderPairs.hpp"
#include "PngImage.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace Align5;

#include "tuple.hpp"



void shasta::align5(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug)
{
    switch(options.m) {
    case 1:
        align5<1>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    case 2:
        align5<2>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    case 3:
        align5<3>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    case 4:
        align5<4>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    default:
        SHASTA_ASSERT(0);
    }
}



// Version templated on m, the number of markers that define
// a "feature" used in the alignment.
template<uint64_t m> void shasta::align5(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug)
{
    Align5::Aligner<m> graph(markers0, markers1,
        options, alignment, alignmentInfo,
        debug);
}



template<uint64_t m> shasta::Align5::Aligner<m>::Aligner(
    const MarkerSequence& markerSequence0,
    const MarkerSequence& markerSequence1,
    const Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug) :
    nx(uint32_t(markerSequence0.size())),
    ny(uint32_t(markerSequence1.size())),
    deltaX(int32_t(options.deltaX)),
    deltaY(int32_t(options.deltaY))
{
    if(debug) {
        cout << timestamp << "Align5 begins." << endl;
        cout << timestamp << "Input sequences have " <<
            nx << " and " << ny << " markers." << endl;
    }

    // Check that we are in the templated version consistent with
    // the options.
    SHASTA_ASSERT(options.m == m);

    // Create features sorted by feature.
    if(debug) {
        cout << timestamp << "Creating sorted features." << endl;
    }
    sortFeatures(markerSequence0, sortedFeatures0);
    sortFeatures(markerSequence1, sortedFeatures1);

    // Create cells.
    if(debug) {
        cout << timestamp << "Creating cells." << endl;
    }
    createCells();
    if(debug) {
        cout << timestamp << "Counting occupied cells." << endl;
        uint64_t occupiedCellCount = 0;
        for(const Cell& cell: cells) {
            if(cell.featureCount > 0) {
                ++occupiedCellCount;
            }
        }
        cout << timestamp << "Found " << occupiedCellCount <<
            " occupied cells out of " << cellCountX * cellCountY << " total." << endl;
        cout << timestamp << "Writing cells." << endl;
        writeCells("Align5-Cells.png");
    }

    if(debug) {
        cout << timestamp << "Writing alignment matrix in feature space." << endl;
        writeAlignmentMatrixInFeatureSpace("Align5-AlignmentMatrixInFeatureSpace.png");
    }

    // Create markers sorted by KmerId.
    if(debug) {
        cout << timestamp << "Creating sorted markers." << endl;
    }
    sortMarkers(markerSequence0, sortedMarkers0);
    sortMarkers(markerSequence1, sortedMarkers1);

    if(debug) {
        cout << timestamp << "Writing alignment matrix in marker space." << endl;
        writeAlignmentMatrixInMarkerSpace("Align5-AlignmentMatrixInMarkerSpace.png");
    }

    if(debug) {
        cout << timestamp << "Align5 ends." << endl;
    }
}



template<uint64_t m> void shasta::Align5::Aligner<m>::sortFeatures(
    const MarkerSequence& markerSequence,
    vector< pair<Feature, uint32_t> >& sortedFeatures)
{
    sortedFeatures.clear();

    // If the sequence is too short, it has no features.
    if(markerSequence.size() < m) {
        return;
    }
    sortedFeatures.reserve(markerSequence.size() - (m-1));

    // Start with the first m kmerIds.
    Feature feature;
    uint64_t nextOrdinal = 0;
    for(uint64_t i=0; i<m; i++) {
        feature[i] = markerSequence[nextOrdinal++].kmerId;
    }


    // Add the features.
    for(uint32_t ordinal=0; ; ordinal++) {
        sortedFeatures.push_back(make_pair(feature, ordinal));

        // Check if done.
        if(nextOrdinal == markerSequence.size()) {
            break;
        }

        // Shift by one.
        for(uint64_t i=1; i<m; i++) {
            feature[i-1] = feature[i];
        }
        feature.back() = markerSequence[nextOrdinal++].kmerId;
    }
    SHASTA_ASSERT(sortedFeatures.size() == markerSequence.size() - (m-1));

    // Sort by feature.
    sort(sortedFeatures.begin(), sortedFeatures.end(),
        OrderPairsByFirstOnly<Feature, uint32_t>());
}



template<uint64_t m> void shasta::Align5::Aligner<m>::sortMarkers(
    const MarkerSequence& markerSequence,
    vector< pair<KmerId, uint32_t> >& sortedMarkers)
{
    sortedMarkers.resize(markerSequence.size());
    for(uint32_t i=0; i<markerSequence.size(); i++) {
        sortedMarkers[i] = make_pair(markerSequence[i].kmerId, i);
    }
    sort(sortedMarkers.begin(), sortedMarkers.end(),
        OrderPairsByFirstOnly<KmerId, uint32_t>());
}



template<uint64_t m> void shasta::Align5::Aligner<m>::
    writeAlignmentMatrixInMarkerSpace(const string& fileName) const
{
    PngImage image(nx, ny);
    writeCheckerboard(image);
    image.writeGrid(   10,  15,  15,  15);      // Grey
    image.writeGrid(   50,  30,  30,  30);      // Grey
    image.writeGrid(  100,  90,  90,  90);      // Grey
    image.writeGrid(  500, 160, 160, 160);      // Grey
    image.writeGrid( 1000, 255, 255, 255);      // White
    image.writeGrid( 5000, 255, 120, 255);      // Purple
    image.writeGrid(10000, 255, 255,  60);      // Yellow
    image.writeGrid(50000, 255, 255, 120);      // Yellow



    // Joint loop over the markers, looking for common k-mer ids.
    auto begin0 = sortedMarkers0.begin();
    auto begin1 = sortedMarkers1.begin();
    auto end0 = sortedMarkers0.end();
    auto end1 = sortedMarkers1.end();

    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->first < it1->first) {
            ++it0;
        } else if(it1->first < it0->first) {
            ++it1;
        } else {

            // We found a common k-mer id.
            const KmerId kmerId = it0->first;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer in kmers0 and kmers1.
            auto it0Begin = it0;
            auto it1Begin = it1;
            auto it0End = it0Begin;
            auto it1End = it1Begin;
            while(it0End!=end0 && it0End->first == kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->first == kmerId) {
                ++it1End;
            }

            // Loop over pairs in the streaks.
            for(auto jt0=it0Begin; jt0!=it0End; ++jt0) {
                for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const uint32_t x = jt0->second;
                    const uint32_t y = jt1->second;
                    image.setPixel(x, y, 255, 0, 0);
                }
            }

            // Continue the joint loop over k-mers.
            it0 = it0End;
            it1 = it1End;
        }

    }


    image.write(fileName);
}



template<uint64_t m> void shasta::Align5::Aligner<m>::
    writeAlignmentMatrixInFeatureSpace(const string& fileName) const
{
    PngImage image(nx, ny);
    writeCheckerboard(image);
    image.writeGrid(   10,  15,  15,  15);      // Grey
    image.writeGrid(   50,  30,  30,  30);      // Grey
    image.writeGrid(  100,  90,  90,  90);      // Grey
    image.writeGrid(  500, 160, 160, 160);      // Grey
    image.writeGrid( 1000, 255, 255, 255);      // White
    image.writeGrid( 5000, 255, 120, 255);      // Purple
    image.writeGrid(10000, 255, 255,  60);      // Yellow
    image.writeGrid(50000, 255, 255, 120);      // Yellow



    // Joint loop over the features, looking for common k-mer ids.
    auto begin0 = sortedFeatures0.begin();
    auto begin1 = sortedFeatures1.begin();
    auto end0 = sortedFeatures0.end();
    auto end1 = sortedFeatures1.end();

    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->first < it1->first) {
            ++it0;
        } else if(it1->first < it0->first) {
            ++it1;
        } else {

            // We found a common k-mer id.
            const Feature feature = it0->first;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer in kmers0 and kmers1.
            auto it0Begin = it0;
            auto it1Begin = it1;
            auto it0End = it0Begin;
            auto it1End = it1Begin;
            while(it0End!=end0 && it0End->first == feature) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->first == feature) {
                ++it1End;
            }

            // Loop over pairs in the streaks.
            for(auto jt0=it0Begin; jt0!=it0End; ++jt0) {
                for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const uint32_t x = jt0->second;
                    const uint32_t y = jt1->second;
                    image.setPixel(x, y, 255, 0, 0);
                }
            }

            // Continue the joint loop over features.
            it0 = it0End;
            it1 = it1End;
        }

    }


    image.write(fileName);
}



template<uint64_t m> void
    shasta::Align5::Aligner<m>::writeCheckerboard(PngImage& image) const
{
    for(int32_t y=0; y<ny; y++) {
        for(int32_t x=0; x<nx; x++) {
            uint32_t iX, iY;
            tie(iX, iY) = getCellIndexes(x, y);
            if(((iX + iY) %2) == 0) {
                image.setPixel(x, y, 0, 48, 0);
            }
        }
    }
}



// Given ordinals x, y, return coordinates iX, iY of the containing cell.
template<uint64_t m> pair<uint32_t, uint32_t>
    shasta::Align5::Aligner<m>::getCellIndexes(uint32_t x, uint32_t y) const
{
    const uint32_t X = x + y;
    const uint32_t Y = nx + y - x - 1;
    const uint32_t iX = X / deltaX;
    const uint32_t iY = Y / deltaY;
    if(iX >= cellCountX) {
        cout << x << " " << y << " " << X << " " << Y << " " << iX << " " << iY << endl;
    }
    SHASTA_ASSERT(iX < cellCountX);
    SHASTA_ASSERT(iY < cellCountY);
    return make_pair(iX, iY);
}


template<uint64_t m> void shasta::Align5::Aligner<m>::createCells()
{
    // Number of X and Y values.
    const uint32_t n = nx + ny - 1;

    // Number of cells in X and Y direction.
    cellCountX = (n-1) / deltaX + 1;
    cellCountY = (n-1) / deltaY + 1;

    // Resize the cell vector.
    cells.clear();
    cells.resize(cellCountX * cellCountY);

    // Count the number of matrix elements of in each cell
    // (for the alignment matrix in feature space).
    // Joint loop over the features, looking for common k-mer ids.
    auto begin0 = sortedFeatures0.begin();
    auto begin1 = sortedFeatures1.begin();
    auto end0 = sortedFeatures0.end();
    auto end1 = sortedFeatures1.end();

    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->first < it1->first) {
            ++it0;
        } else if(it1->first < it0->first) {
            ++it1;
        } else {

            // We found a common feature.
            const Feature feature = it0->first;


            // This feature could appear more than once in each of the sequences,
            // so we need to find the streak of this feature.
            auto it0Begin = it0;
            auto it1Begin = it1;
            auto it0End = it0Begin;
            auto it1End = it1Begin;
            while(it0End!=end0 && it0End->first == feature) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->first == feature) {
                ++it1End;
            }

            // Loop over pairs in the streaks.
            for(auto jt0=it0Begin; jt0!=it0End; ++jt0) {
                for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const uint32_t x = jt0->second;
                    const uint32_t y = jt1->second;
                    const uint32_t X = x + y;
                    const uint32_t Y = nx + y - x - 1;
                    const uint32_t iX = X / deltaX;
                    const uint32_t iY = Y / deltaY;
                    Cell& cell = getCell(iX, iY);
                    cell.featureCount++;
                    cell.minX = min(cell.minX, X);
                    cell.maxX = max(cell.maxX, X);
                    cell.minY = min(cell.minY, Y);
                    cell.maxY = max(cell.maxY, Y);
                }
            }

            // Continue the joint loop over features.
            it0 = it0End;
            it1 = it1End;
        }

    }

}



template<uint64_t m> void shasta::Align5::Aligner<m>::writeCells(const string& fileName) const
{
    const uint32_t markersPerPixel = 4;
    const uint32_t n = nx + ny - 1;
    PngImage image(n/markersPerPixel, n/markersPerPixel);
    //  cout << "Cells image size " << n/markersPerPixel << endl;

    for(uint32_t iX=0; iX<cellCountX; iX++) {
        for(uint32_t iY=0; iY<cellCountY; iY++) {
            const Cell& cell = getCell(iX, iY);

            // If the cell is empty, skip it.
            if(cell.featureCount == 0) {
                continue;
            }

            // Display the bounding box.
            const uint32_t minX = cell.minX / markersPerPixel;
            const uint32_t maxX = cell.maxX / markersPerPixel;
            const uint32_t minY = cell.minY / markersPerPixel;
            const uint32_t maxY = cell.maxY / markersPerPixel;
            /*
            cout <<
                iX << " " <<
                iY << " " <<
                minX << " " <<
                maxX << " " <<
                minY << " " <<
                maxY << " " <<
                endl;
            */
            for(uint32_t i=minX; i<=maxX; i++) {
                for(uint32_t j=minY; j<=maxY; j++) {
                    image.setPixel(i, j, 255, 255, 255);
                }
            }
        }

    }

    image.write(fileName);
}
