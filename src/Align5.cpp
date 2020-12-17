// Shasta
#include "Align5.hpp"
#include "Alignment.hpp"
#include "countingSort.hpp"
#include "hashArray.hpp"
#include "orderPairs.hpp"
#include "PngImage.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace Align5;

// Seqan.
#include <seqan/align.h>

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "chrono.hpp"
#include "fstream.hpp"
#include <map>
#include <stack>
#include "tuple.hpp"
#include <unordered_map>



void shasta::align5(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const Options& options,
    MemoryMapped::ByteAllocator& byteAllocator,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug)
{
    Align5::Aligner graph(markers0, markers1,
        options, byteAllocator, alignment, alignmentInfo,
        debug);
}



Aligner::Aligner(
    const MarkerSequence& markerSequence0,
    const MarkerSequence& markerSequence1,
    const Options& options,
    MemoryMapped::ByteAllocator& byteAllocator,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug) :
    nx(uint32_t(markerSequence0.size())),
    ny(uint32_t(markerSequence1.size())),
    deltaX(int32_t(options.deltaX)),
    deltaY(int32_t(options.deltaY)),
    byteAllocator(byteAllocator)
{
    // Parameters to expose when code stabilizes.
    const uint32_t minEntryCountPerCell = 10;
    const uint32_t maxDistanceFromBoundary = deltaX / 2;



    if(debug) {
        cout << timestamp << "Align5 begins." << endl;
        cout << timestamp << "Input sequences have " <<
            nx << " and " << ny << " markers." << endl;
    }

    // Create markers sorted by KmerId.
    if(debug) {
        cout << timestamp << "Creating sorted markers." << endl;
    }
    storeMarkers(markerSequence0, markers[0], sortedMarkers[0]);
    storeMarkers(markerSequence1, markers[1], sortedMarkers[1]);

    // Create the sparse representation of the alignment matrix.
    if(debug) {
        cout << timestamp << "Creating the alignment matrix." << endl;
    }
    // const auto t0 = steady_clock::now();
    createAlignmentMatrix();
    /*
    const auto t1 = steady_clock::now();
    if(debug) {
        cout << "Computation of alignment matrix took " << seconds(t1-t0) << " s." << endl;
    }
    */

    // Gather well populated cells.
    if(debug) {
        cout << timestamp << "Creating cells." << endl;
    }
    createCells(minEntryCountPerCell, maxDistanceFromBoundary);

    // Forward search.
    // Find cells that are forward accessible from the left/top.
    if(debug) {
        cout << timestamp << "Forward search begins." << endl;
    }
    forwardSearch();

    // Backward search.
    if(debug) {
        cout << timestamp << "Backward search begins." << endl;
    }
    backwardSearch();

    // Group active cells by connected component.
    if(debug) {
        cout << timestamp << "Grouping active cells by connected component." << endl;
    }
    findActiveCellsConnectedComponents();

    // Compute a banded alignment for each connected component of
    // active cells.
    if(debug) {
        cout << timestamp << "Computing namded alignments." << endl;
    }
    computeBandedAlignments(debug);


    if(debug) {
        cout << timestamp << "Writing cells." << endl;
        writeCellsPng("Align5-Cells.png");
        writeCellsCsv("Align5-Cells.csv");
    }

    if(debug) {
        cout << timestamp << "Writing the alignment matrix." << endl;
        writeAlignmentMatrixPng("Align5-AlignmentMatrix.png", maxDistanceFromBoundary);
        writeAlignmentMatrixCsv("Align5-AlignmentMatrix.csv");
    }

    if(debug) {
        cout << timestamp << "Align5 ends." << endl;
    }
}


void Aligner::storeMarkers(
    const MarkerSequence& markerSequence,
    vector<KmerId>& markers,
    vector< pair<KmerId, uint32_t> >& sortedMarkers)
{
    const uint64_t n = markerSequence.size();
    markers.resize(n);
    sortedMarkers.resize(n);
    for(uint32_t i=0; i<n; i++) {
        const KmerId kmerId = markerSequence[i].kmerId;
        markers[i] = kmerId;
        sortedMarkers[i] = make_pair(markerSequence[i].kmerId, i);
    }
    sort(sortedMarkers.begin(), sortedMarkers.end(),
        OrderPairsByFirstOnly<KmerId, uint32_t>());
}



// Return (X,Y) given (x,y).
Coordinates Aligner::getXY(Coordinates xy) const
{
    return Coordinates(
        xy.first + xy.second,
        nx + xy.second - xy.first - 1
        );
}



// Convert an arbitrary (X,Y) to (x,y).
// If the point is outside the alignment matrix,
// we can end up with negative values.
SignedCoordinates Aligner::getxy(Coordinates XY) const
{
    const int32_t X = int32_t(XY.first);
    const int32_t Y = int32_t(XY.second);
    const int32_t x = (X - Y + int32_t(nx) - 1) / 2;
    const int32_t y = (X + Y - int32_t(nx) + 1) / 2;
    return SignedCoordinates(x, y);
}



void Aligner::createAlignmentMatrix()
{
    alignmentMatrix.clear();

    // Joint loop over the sorted markers, looking for common markers.
    auto begin0 = sortedMarkers[0].begin();
    auto begin1 = sortedMarkers[1].begin();
    auto end0 = sortedMarkers[0].end();
    auto end1 = sortedMarkers[1].end();

    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->first < it1->first) {
            ++it0;
        } else if(it1->first < it0->first) {
            ++it1;
        } else {

            // We found a common KmerId.
            const KmerId kmerId = it0->first;


            // This KmerId could appear more than once in each of the sequences,
            // so we need to find the streak of this KmerId.
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
                const uint32_t x = jt0->second;
                for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const uint32_t y = jt1->second;
                    const Coordinates xy(x, y);
                    const Coordinates iXY = getCellIndexesFromxy(xy);
                    const uint32_t iY = iXY.second;
                    if(alignmentMatrix.size() <= iY) {
                        alignmentMatrix.resize(iY+1,
                            AlignmentMatrixEntryVector(0, AlignmentMatrixAllocator(byteAllocator)));
                    }
                    const uint32_t iX = iXY.first;
                    alignmentMatrix[iY].push_back(make_pair(iX, xy));
                }
            }

            // Continue the joint loop over KmerId's.
            it0 = it0End;
            it1 = it1End;
        }

    }

    // For each iY value, sort by iX.
    // Except for short vectors, counting sort is faster.
    vector<uint32_t> count;
    AlignmentMatrixEntryVector w(byteAllocator);
    for(auto& v: alignmentMatrix) {
        if(false /*v.size() < 5*/) {
            sort(v.begin(), v.end(), OrderPairsByFirstOnly<uint32_t, Coordinates>());
        } else {
            countingSort(v, count, w);
        }
    }

}



void Aligner::writeAlignmentMatrixCsv(const string& fileName) const
{
    uint64_t entryCount = 0;
    ofstream csv(fileName);
    csv << "iX,iY,X,Y,x,y\n";
    for(uint32_t iY=0; iY<alignmentMatrix.size(); iY++) {
        for(const auto& v: alignmentMatrix[iY]) {
            const uint32_t iX = v.first;
            const Coordinates& xy = v.second;
            const Coordinates XY = getXY(xy);
            csv << iX << ",";
            csv << iY << ",";
            csv << XY.first << ",";
            csv << XY.second << ",";
            csv << xy.first << ",";
            csv << xy.second << "\n";
            SHASTA_ASSERT(Coordinates(iX, iY) == getCellIndexesFromxy(xy));
            SHASTA_ASSERT(Coordinates(iX, iY) == getCellIndexesFromXY(XY));
            ++entryCount;
        }
    }
    cout << "The alignment matrix contains " << entryCount << " entries." << endl;
}



void Aligner::writeAlignmentMatrixPng(
    const string& fileName,
    uint32_t maxDistanceFromBoundary) const
{
    PngImage image(nx, ny);

    writeCheckerboard(image, maxDistanceFromBoundary);

    image.writeGrid(   10,  15,  15,  15);      // Grey
    image.writeGrid(   50,  30,  30,  30);      // Grey
    image.writeGrid(  100,  90,  90,  90);      // Grey
    image.writeGrid(  500, 160, 160, 160);      // Grey
    image.writeGrid( 1000, 255, 255, 255);      // White
    image.writeGrid( 5000, 255, 120, 255);      // Purple
    image.writeGrid(10000, 255, 255,  60);      // Yellow
    image.writeGrid(50000, 255, 255, 120);      // Yellow

    for(uint32_t iY=0; iY<alignmentMatrix.size(); iY++) {
        for(const auto& v: alignmentMatrix[iY]) {
            const Coordinates& xy = v.second;
            const uint32_t x = xy.first;
            const uint32_t y = xy.second;
            image.setPixel(x, y, 255, 0, 0);
        }
    }

    image.write(fileName);
}



void Aligner::writeCheckerboard(
    PngImage& image,
    uint32_t maxDistanceFromBoundary) const
{
    Coordinates xy;
    uint32_t& x = xy.first;
    uint32_t& y = xy.second;

    Coordinates iXY;
    uint32_t& iX = iXY.first;
    uint32_t& iY = iXY.second;

    for(y=0; y<ny; y++) {
        for(x=0; x<nx; x++) {
            iXY = getCellIndexesFromxy(xy);
            /*
            const bool isNearLeftOrTop =
                (cellDistanceFromLeft(iXY) < maxDistanceFromBoundary) or
                (cellDistanceFromTop(iXY)  < maxDistanceFromBoundary);
            const bool isNearRightOrBottom =
                (cellDistanceFromRight(iXY)  < maxDistanceFromBoundary) or
                (cellDistanceFromBottom(iXY) < maxDistanceFromBoundary);
            */
            const Cell* cell = findCell(iXY);
            const bool isEvenCell = (((iX + iY) % 2) == 0);

            int r = 0;
            int g = 0;
            int b = 0;
            if(cell) {
                if(cell->isActive()) {
                    g = 255;
                } else if(cell->isForwardAccessible) {
                    b = 255;
                } else {
                    r = 128;
                    g = 128;
                    b = 128;
                }
            } else {
                if(isEvenCell) {
                    r = 64;
                }
            }

            image.setPixel(x, y, r, g, b);
        }
    }
}



void Aligner::createCells(
    uint32_t minEntryCountPerCell,
    uint32_t maxDistanceFromBoundary)
{
    // Start with nothing.
    cells.clear();
    cells.resize(alignmentMatrix.size());

    // Loop over iY values.
    for(uint32_t iY=0; iY<alignmentMatrix.size(); iY++) {

        // Access the vectors for this value of iY.
        const auto& iYAlignmentMatrix = alignmentMatrix[iY];
        vector< pair<uint32_t, Cell> >& iYCells = cells[iY];

        // Each vector in the alignment matrix is sorted by iX, so we can scan it,
        // creating a new cell each time we encounter a new value of iX
        // (but only if the cell is sufficiently populated).
        for(auto it=iYAlignmentMatrix.begin(); it!=iYAlignmentMatrix.end(); /* Increment later */) {

            // When getting here, we are starting a new cell.
            const auto& entry = *it;
            const uint32_t iX = entry.first;

            // Now count the the alignment matrix entries in this cell.
            const auto it0 = it;
            while(true) {
                ++it;
                if(it == iYAlignmentMatrix.end()) {
                    break;
                }
                if(it->first != iX) {
                    break;
                }
            }

            // If too few entries in this cell, don't store it.
            if(it - it0 < minEntryCountPerCell) {
                continue;
            }

            // Construct the cell.
            const Coordinates iXY(iX, iY);
            Cell cell;
            cell.isNearLeftOrTop =
                (cellDistanceFromLeft(iXY) < maxDistanceFromBoundary) or
                (cellDistanceFromTop(iXY)  < maxDistanceFromBoundary);
            cell.isNearRightOrBottom =
                (cellDistanceFromRight(iXY)  < maxDistanceFromBoundary) or
                (cellDistanceFromBottom(iXY) < maxDistanceFromBoundary);

            // Store this cell.
            iYCells.push_back(make_pair(iX, cell));
        }

    }
}



void Aligner::writeCellsCsv(
    const string& fileName) const
{
    uint64_t cellCount = 0;
    ofstream csv(fileName);
    csv << "iX,iY,minX,maxX,minY,maxY,sizeX,sizeY\n";
    for(uint32_t iY=0; iY<cells.size(); iY++) {
        for(const auto& p: cells[iY]) {
            const uint32_t iX = p.first;
            csv << iX << ",";
            csv << iY << "\n";
            ++cellCount;
        }
    }
    cout << "There are " << cellCount << " cells." << endl;
}



void Aligner::writeCellsPng(
    const string& fileName) const
{
    // The size of the square in XY space.
    const uint32_t sizeXY = nx + ny -1;

    // Number of markers per pixel - to control the size of the picture.
    const uint32_t markersPerPixel = 5;

    // Create the image.
    const uint32_t imageSize = sizeXY / markersPerPixel;
    PngImage image(imageSize, imageSize);

    // Write the checkerboard in XY space.
    for(uint32_t j=0; j<imageSize; j++) {
        const uint32_t Y = j * markersPerPixel;
        for(uint32_t i=0; i<imageSize; i++) {
            const uint32_t X = i * markersPerPixel;
            const Coordinates iXY = getCellIndexesFromXY(Coordinates(X, Y));
            const uint32_t iX = iXY.first;
            const uint32_t iY = iXY.second;
            if(((iX+iY) %2) == 0) {
                image.setPixel(i, j, 0, 50, 0);
            }
        }
    }

    // Write the cells.
    for(uint32_t iY=0; iY<cells.size(); iY++) {
        for(const auto& p: cells[iY]) {
            const uint32_t iX = p.first;
            const Cell& cell = p.second;
            SHASTA_ASSERT(iX < sizeXY);
            SHASTA_ASSERT(iY < sizeXY);
            const uint32_t iMin = ( iX    * deltaX)  / markersPerPixel;
            const uint32_t iMax = ((iX+1) * deltaX)  / markersPerPixel;
            const uint32_t jMin = ( iY    * deltaY)  / markersPerPixel;
            const uint32_t jMax = ((iY+1) * deltaY)  / markersPerPixel;

            int r = 0;
            int g = 0;
            int b = 0;
            if(cell.isActive()) {
                g = 255;
            } else if(cell.isForwardAccessible) {
                b = 255;
           } else {
                r = 128;
                g = 128;
                b = 128;
            }

            for(uint32_t j=jMin; j<jMax; j++) {
                for(uint32_t i=iMin; i<iMax; i++) {
                    if((i < imageSize) and (j < imageSize)) {
                        image.setPixel(i, j, r, g, b);
                    }
                }
            }
         }
    }

    // Write the image.
    image.write(fileName);
}



// Return the distance of a cell from the left boundary
// of the alignment matrix, or 0 if the cell
// is partially or entirely to the left of that boundary.
uint32_t Aligner::cellDistanceFromLeft(const Coordinates& iXY) const
{
    const uint32_t iX = iXY.first;
    const uint32_t iY = iXY.second;

    // The bottom left corner in (XY) space
    // is closest to the left boundary of the alignment matrix.
    const Coordinates XY = {iX * deltaX, (iY + 1) * deltaY};

    // Convert to (xy).
    const SignedCoordinates xy = getxy(XY);
    const int32_t x = xy.first;

    if(x < 0) {
        return 0;
    } else {
        return uint32_t(x);
    }

}



// Return the distance of a cell from the right boundary
// of the alignment matrix, or 0 if the cell
// is partially or entirely to the right of that boundary.
uint32_t  Aligner::cellDistanceFromRight(const Coordinates& iXY) const
{
    const uint32_t iX = iXY.first;
    const uint32_t iY = iXY.second;

    // The top right corner in (XY) space
    // is closest to the right boundary of the alignment matrix.
    const Coordinates XY = {(iX + 1) * deltaX, iY * deltaY};

    // Convert to (xy).
    const SignedCoordinates xy = getxy(XY);
    const int32_t x = xy.first;

    if(x >= int32_t(nx) - 1) {
        return 0;
    } else {
        return uint32_t(nx - 1 - uint32_t(x));
    }
}



// Return the distance of a cell from the top boundary
// of the alignment matrix, or 0 if the cell
// is partially or entirely above that boundary.
uint32_t Aligner::cellDistanceFromTop(const Coordinates& iXY) const
{
    const uint32_t iX = iXY.first;
    const uint32_t iY = iXY.second;

    // The top left corner in (XY) space
    // is closest to the left boundary of the alignment matrix.
    const Coordinates XY = {iX * deltaX, iY * deltaY};

    // Convert to (xy).
    const SignedCoordinates xy = getxy(XY);
    const int32_t y = xy.second;

    if(y < 0) {
        return 0;
    } else {
        return uint32_t(y);
    }

}



// Return the distance of a cell from the bottom boundary
// of the alignment matrix, or 0 if the cell
// is partially or entirely below that boundary.
uint32_t Aligner::cellDistanceFromBottom(const Coordinates& iXY) const
{
    const uint32_t iX = iXY.first;
    const uint32_t iY = iXY.second;

    // The bottom right corner in (XY) space
    // is closest to the left boundary of the alignment matrix.
    const Coordinates XY = {(iX + 1) * deltaX, (iY + 1) * deltaY};

    // Convert to (xy).
    const SignedCoordinates xy = getxy(XY);
    const int32_t y = xy.second;

    if(y >= int32_t(ny) - 1) {
        return 0;
    } else {
        return uint32_t(ny - 1 - uint32_t(y));
    }

}



// Find a cell with given (iX,iY).
Aligner::Cell* Aligner::findCell(const Coordinates& iXY)
{
    const uint32_t iX = iXY.first;
    const uint32_t iY = iXY.second;
    if(iY >= cells.size()) {
        return 0;
    }
    vector< pair<uint32_t, Cell> >& v = cells[iY];

    // Look for a cell with this iX.
    auto it = std::lower_bound(v.begin(), v.end(),
        make_pair(iX, Cell()),
        OrderPairsByFirstOnly<uint32_t, Cell>());
    if(it == v.end()) {
        return 0;
    }
    if(it->first != iX) {
        return 0;
    }
    return &(it->second);
}



// const version of the above
const Aligner::Cell* Aligner::findCell(const Coordinates& iXY) const
{
    const uint32_t iX = iXY.first;
    const uint32_t iY = iXY.second;
    if(iY >= cells.size()) {
        return 0;
    }
    const vector< pair<uint32_t, Cell> >& v = cells[iY];

    // Look for a cell with this iX.
    auto it = std::lower_bound(v.begin(), v.end(),
        make_pair(iX, Cell()),
        OrderPairsByFirstOnly<uint32_t, Cell>());
    if(it == v.end()) {
        return 0;
    }
    if(it->first != iX) {
        return 0;
    }
    return &(it->second);
}



// DFS in cell space, starting from cells near the left/top.
// We use DFS instead of BFS for better locality of memory access.
void Aligner::forwardSearch()
{
    uint64_t n = 0;

    // Initialize the stack of undiscovered iXY.
    std::stack<Coordinates> s;
    for(uint32_t iY=0; iY<cells.size(); iY++) {
        vector< pair<uint32_t, Cell> >& v = cells[iY];
        for(auto& p: v) {
            Cell& cell = p.second;
            if(cell.isNearLeftOrTop) {
                cell.isForwardAccessible = 1;
                ++n;
                const uint32_t iX = p.first;
                s.push(Coordinates(iX, iY));
            }
        }
    }

    // DFS.
    vector<Coordinates> children;
    while(not s.empty()) {
        const Coordinates& iXY0 = s.top();
        const uint32_t iX0 = iXY0.first;
        const uint32_t iY0 = iXY0.second;
        s.pop();

        // Loop over possible children.
        for(int32_t dY=-1; dY<=1; dY++) {
            const int32_t iY1Signed = int32_t(iY0) + dY;
            if(iY1Signed < 0) {
                continue;
            }
            const uint32_t iY1 = uint32_t(iY1Signed);
            for(uint32_t dX=0; dX<=1; dX++) {
                const uint32_t iX1 = iX0 + dX;
                const Coordinates iXY1(iX1, iY1);
                Cell* cell1 = findCell(iXY1);
                if(cell1 and not cell1->isForwardAccessible) {
                    cell1->isForwardAccessible = 1;
                    ++n;
                    s.push(iXY1);
                }
            }
        }
    }
    cout << "Forward search found " << n << " cells." << endl;
}



// Backward DFS in cell space, starting from cells near the
// right/bottom that are also forward accessible from the left/top.
// We use DFS instead of BFS for better locality of memory access.
void Aligner::backwardSearch()
{
    uint64_t n = 0;

    // Initialize the stack of undiscovered iXY.
    std::stack<Coordinates> s;
    for(uint32_t iY=0; iY<cells.size(); iY++) {
        vector< pair<uint32_t, Cell> >& v = cells[iY];
        for(auto& p: v) {
            Cell& cell = p.second;
            if(cell.isNearRightOrBottom and cell.isForwardAccessible) {
                cell.isBackwardAccessible = 1;
                ++n;
                const uint32_t iX = p.first;
                s.push(Coordinates(iX, iY));
            }
        }
    }

    // DFS.
    vector<Coordinates> children;
    while(not s.empty()) {
        const Coordinates& iXY0 = s.top();
        const uint32_t iX0 = iXY0.first;
        const uint32_t iY0 = iXY0.second;
        s.pop();

        // Loop over possible parents.
        for(int32_t dY=-1; dY<=1; dY++) {
            const int32_t iY1Signed = int32_t(iY0) + dY;
            if(iY1Signed < 0) {
                continue;
            }
            const uint32_t iY1 = uint32_t(iY1Signed);
            for(int32_t dX=-1; dX<=0; dX++) {
                const int32_t iX1Signed = int32_t(iX0) + dX;
                if(iX1Signed < 0) {
                    continue;
                }
                const uint32_t iX1 = uint32_t(iX1Signed);
                const Coordinates iXY1(iX1, iY1);
                Cell* cell1 = findCell(iXY1);
                if(cell1 and not cell1->isBackwardAccessible) {
                    cell1->isBackwardAccessible = 1;
                    ++n;
                    s.push(iXY1);
                }
            }
        }
    }
    cout << "Backward search found " << n << " cells." << endl;
}



// Group active cells in connected component.
void Aligner::findActiveCellsConnectedComponents()
{
    // Gather all the active cells.
    // Store a contiguously numbered id for each of them.
    // The id will be used for the connected component computation below.
    uint32_t nextCellId = 0;
    std::unordered_map<Coordinates, uint32_t, HashTuple<Coordinates> > activeCells;
    for(uint32_t iY=0; iY<cells.size(); iY++) {
        const vector< pair<uint32_t, Cell> >& iYCells = cells[iY];
        for(const pair<uint32_t, Cell>& p: iYCells) {
            const Cell& cell = p.second;
            if(cell.isActive()) {
                const uint32_t iX = p.first;
                activeCells.insert(make_pair(Coordinates(iX, iY), nextCellId++));
            }
        }
    }
    const uint32_t activeCellCount = nextCellId;



    // Compute the connected components.
    vector<uint32_t> rank(activeCellCount);
    vector<uint32_t> parent(activeCellCount);
    boost::disjoint_sets<uint32_t*, uint32_t*> disjointSets(&rank[0], &parent[0]);
    for(uint32_t i=0; i<activeCellCount; i++) {
        disjointSets.make_set(i);
    }
    for(const auto& p: activeCells) {
        const Coordinates& iXY0 = p.first;
        const uint32_t iX0 = iXY0.first;
        const uint32_t iY0 = iXY0.second;
        const uint32_t cellId0 = p.second;

        // Loop over possible neighbors.
        for(int32_t dY=-1; dY<=1; dY++) {
            const int32_t iY1Signed = int32_t(iY0) + dY;
            if(iY1Signed < 0) {
                continue;
            }
            const uint32_t iY1 = uint32_t(iY1Signed);
            for(int32_t dX=-1; dX<=1; dX++) {
                if((dX == 0) and (dY == 0)) {
                    continue;
                }
                const int32_t iX1Signed = int32_t(iX0) + dX;
                if(iX1Signed < 0) {
                    continue;
                }
                const uint32_t iX1 = uint32_t(iX1Signed);
                const auto it = activeCells.find(Coordinates(iX1, iY1));
                if(it == activeCells.end()) {
                    continue;
                }
                const uint32_t cellId1 = it->second;
                disjointSets.union_set(cellId0, cellId1);
            }
        }
    }


    // Gather the cells in each connected component.
    std::map<uint32_t, vector<Coordinates> > connectedComponents;
    for(const auto& p: activeCells) {
        const Coordinates& iXY = p.first;
        const uint32_t cellId = p.second;
        const uint32_t componentId = disjointSets.find_set(cellId);
        connectedComponents[componentId].push_back(iXY);
    }
    cout << "Found " << connectedComponents.size() << " connected components of sizes:";
    activeCellsConnectedComponents.clear();
    for(const auto& p: connectedComponents) {
        cout << " " << p.second.size();
        activeCellsConnectedComponents.push_back(p.second);
    }
    cout << endl;
}



// Compute a banded alignment for each connected component of
// active cells.
void Aligner::computeBandedAlignments(bool debug) const
{

    // Loop over connected components of active cells.
    vector< pair<bool, bool> > alignment;
    for(const vector<Coordinates>& component: activeCellsConnectedComponents) {
        if(debug) {
            cout << "Connected component with " << component.size() << " active cells." << endl;
        }

        // Compute the iY range.
        uint32_t iYMin = std::numeric_limits<uint32_t>::max();
        uint32_t iYMax = 0;
        for(const Coordinates& iXY: component) {
            const uint32_t iY = iXY.second;
            iYMin = min(iYMin, iY);
            iYMax = max(iYMax, iY);
        }
        if(debug) {
            cout << "iY range " << iYMin << " " << iYMax << endl;
        }

        // Compute the corresponding y range.
        const uint32_t YMin = iYMin * deltaY;
        const uint32_t YMax = (iYMax+1) * deltaY - 1;
        if(debug) {
            cout << "Y range " << YMin << " " << YMax << endl;
        }

        // Compute the corresponding band.
        const int32_t bandMin = int32_t(nx) -1 - int32_t(YMax);
        const int32_t bandMax = int32_t(nx) -1 - int32_t(YMin);
        const int32_t bandWidth = bandMax - bandMin + 1;
        const int32_t bandCenter = (bandMin + bandMax) / 2;
        const int32_t nominalAlignmentLength =
            min(int32_t(nx), int32_t(ny) + bandCenter) -
            max(0, bandCenter);
        if(debug) {
            cout << "Band " << bandMin << " " << bandMax << endl;
            cout << "Band width " << bandWidth << endl;
            cout << "Band center " << bandCenter << endl;
            cout << "Alignment length at band center " << nominalAlignmentLength << endl;
        }

        // Compute an alignment with this band.
        computeBandedAlignment(bandMin, bandMax, debug);
    }

}



// Compute a banded alignment with a given band.
bool Aligner::computeBandedAlignment(
    int32_t bandMin,
    int32_t bandMax,
    bool debug) const
{
    // Some seqan types and constants we need.
    using namespace seqan;
    using TSequence = String<KmerId>;
    using TStringSet = StringSet<TSequence>;
    using TDepStringSet = StringSet< TSequence, Dependent<> >;
    using TAlignGraph = Graph< seqan::Alignment<TDepStringSet> >;
    const uint32_t seqanGapValue = 45;

    if(debug) {
        cout << timestamp << "Banded alignment computation begins." << endl;
    }

    // Fill in the seqan sequences.
    // Add 100 to kMerIds to prevent collision from the seqan gap value.
    array<TSequence, 2> sequences;
    for(uint64_t i=0; i<2; i++) {
        for(const KmerId kmerId: markers[i]) {
            appendValue(sequences[i], kmerId + 100);
        }
    }

    TStringSet sequencesSet;
    appendValue(sequencesSet, sequences[0]);
    appendValue(sequencesSet, sequences[1]);

    // Compute the banded alignment.
    TAlignGraph graph(sequencesSet);
    const int score = globalAlignment(
        graph,
        Score<int, Simple>(matchScore, mismatchScore, gapScore),
        AlignConfig<true, true, true, true>(),
        bandMin, bandMax,
        LinearGaps());
    if(score == seqan::MinValue<int>::VALUE) {
        cout << "SeqAn banded alignment computation failed." << endl;
        return false;
    } else if(debug) {
        cout << "Alignment score is " << score << endl;
    }

    TSequence align;
    convertAlignment(graph, align);
    const int totalAlignmentLength = int(seqan::length(align));
    SHASTA_ASSERT((totalAlignmentLength % 2) == 0);    // Because we are aligning two sequences.
    const int alignmentLength = totalAlignmentLength / 2;
    if(debug) {
        cout << "Alignment length " << alignmentLength << endl;
    }


    // Fill in the marker alignment.
    Alignment alignment;
    uint32_t ordinal0 = 0;
    uint32_t ordinal1 = 0;
    for(int i=0;
        i<alignmentLength and ordinal0<markers[0].size() and ordinal1<markers[1].size(); i++) {
        if( align[i] != seqanGapValue and
            align[i + alignmentLength] != seqanGapValue and
            markers[0][ordinal0] == markers[1][ordinal1]) {
            alignment.ordinals.push_back(array<uint32_t, 2>{ordinal0, ordinal1});
        }
        if(align[i] != seqanGapValue) {
            ++ordinal0;
        }
        if(align[i + alignmentLength] != seqanGapValue) {
            ++ordinal1;
        }
    }

    // Create the AlignmentInfo.
    const AlignmentInfo alignmentInfo(
        alignment, uint32_t(markers[0].size()), uint32_t(markers[1].size()));
    pair<uint32_t, uint32_t> trim = alignmentInfo.computeTrim();
    cout << "Aligned marker count " << alignmentInfo.markerCount << endl;
    cout << "Aligned marker fraction " <<
        alignmentInfo.alignedFraction(0) << " " <<
        alignmentInfo.alignedFraction(1) << endl;
    cout << "maxSkip " << alignmentInfo.maxSkip << endl;
    cout << "maxDrift " << alignmentInfo.maxDrift << endl;
    cout << "Trim " << trim.first << " " << trim.second << endl;

    if(debug) {
        cout << timestamp << "Banded alignment computation ends." << endl;
    }
    return true;
}

