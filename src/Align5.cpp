#include "Align5.hpp"
#include "orderPairs.hpp"
#include "PngImage.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace Align5;

#include "fstream.hpp"
#include <stack>
#include "tuple.hpp"



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
    sortMarkers(markerSequence0, sortedMarkers0);
    sortMarkers(markerSequence1, sortedMarkers1);

    // Create the sparse representation of the alignment matrix.
    if(debug) {
        cout << timestamp << "Creating the alignment matrix." << endl;
    }
    createAlignmentMatrix();

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
    if(debug) {
        cout << timestamp << "Backward search begins." << endl;
    }
    backwardSearch();

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



void Aligner::sortMarkers(
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



// Return (X,Y) given (x,y).
Coordinates Aligner::getXY(Coordinates xy) const
{
    return Coordinates(
        xy.first + xy.second,
        nx + xy.second - xy.first - 1
        );
}



// Return (iX,iY) given (X,Y).
Coordinates Aligner::getCellIndexesFromXY(Coordinates XY) const
{
    return Coordinates(
        XY.first  / deltaX,
        XY.second / deltaY
        );
}


// Return (iX,iY) given (x,y).
Coordinates Aligner::getCellIndexesFromxy(Coordinates xy) const
{
    const Coordinates XY = getXY(xy);
    return getCellIndexesFromXY(XY);
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
                for(auto jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const uint32_t x = jt0->second;
                    const uint32_t y = jt1->second;
                    const Coordinates iXY = getCellIndexesFromxy(Coordinates(x, y));
                    const uint32_t iX = iXY.first;
                    const uint32_t iY = iXY.second;
                    if(alignmentMatrix.size() <= iY) {
                        alignmentMatrix.resize(iY+1,
                            AlignmentMatrixEntryVector(0, AlignmentMatrixAllocator(byteAllocator)));
                    }
                    alignmentMatrix[iY].push_back(make_pair(iX, Coordinates(x, y)));
                }
            }

            // Continue the joint loop over KmerId's.
            it0 = it0End;
            it1 = it1End;
        }

    }

    for(auto& v: alignmentMatrix) {
        sort(v.begin(), v.end(), OrderPairsByFirstOnly<uint32_t, Coordinates>());
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
                if(cell->isForwardAccessible and cell->isBackwardAccessible) {
                    g = 255;
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
            if(cell.isForwardAccessible and cell.isBackwardAccessible) {
                g = 255;
            } else {
                r = 128;
                g = 128;
                b = 128;
            }

            for(uint32_t j=jMin; j<jMax; j++) {
                for(uint32_t i=iMin; i<iMax; i++) {
                    SHASTA_ASSERT(i < imageSize);
                    SHASTA_ASSERT(j < imageSize);
                    image.setPixel(i, j, r, g, b);
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




