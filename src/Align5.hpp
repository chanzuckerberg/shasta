#ifndef SHASTA_ALIGN5_HPP
#define SHASTA_ALIGN5_HPP

/*******************************************************************************

Marker alignment of two sequences markerSequence0 and markerSequence1,
each defined as a sequence of marker KmerId's.

We call x or y the index (or position or ordinal)
of a marker in markerSequence0 or markerSequence1 respectively, so:
    - markerSequence0[x] is the marker at position x of markerSequence0
    - markerSequence1[y] is the marker at position y of markerSequence1

The number of markers in markerSequence0 is nx
and the number of markers in markerSequence1 is ny.
For any two positions x and y the following hold:
0 <= x <= nx-1
0 <= y <= ny-1

The coordinates in the alignment matrix (x, y)
are represented with x along the horizontal axis and increasing toward the right,
and y along the vertical axis and increasing toward
the bottom. The alignment matrix element at position (x,y)
exists if markerSequence0[x]=markerSequence1[y] when
working with sequences of markers and if
featureSequence0[x]=featureSequence1[y] when working with sequences
of features.

We also use coordinates X and Y defined as:
X = x + y
Y = y + (nx - 1 - x)

It can be verified that:
0 <= X <= nx + ny - 2
0 <= Y <= nx + ny - 2
So the total number of distinct values of X and Y is nx + ny - 1.

X is a coordinate along the diagonal of the alignment matrix,
and Y is orthogonal to it and identifies the diagonal.
In (X,Y) coordinates the alignment matrix is a subset of
the square of size nx + ny -1. The alignment matrix is
rotated by 45 degrees relative to this square.

The inverse transformation is
x = (X - Y + nx - 1) / 2
y = (X + Y - nx + 1) / 2

We use a sparse representation of the alignment matrix
in which non-zero alignment matrix entries are stored
organized by cell in a rectangular arrangement if cells
of size (deltaX, deltaY) in (X,Y) space.

*******************************************************************************/

#include "Marker.hpp"
#include "MemoryMappedAllocator.hpp"
#include "span.hpp"

#include "array.hpp"
#include <limits>
#include <unordered_map>
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    class Alignment;
    class AlignmentInfo;
    class PngImage;

    namespace Align5 {
        class Aligner;
        class MatrixEntry;
        class Options;

        // This is used to store (x,y), (X,Y), or (iX, iY).
        using Coordinates = pair<uint32_t, uint32_t>;

        // When converting an arbitrary (X,Y) to (x,y)
        // we can end up with negative values.
        using SignedCoordinates = pair<uint32_t, uint32_t>;
    }

    namespace MemoryMapped {
        class ByteAllocator;
    }

    void align5(
        const span<const CompressedMarker>&,
        const span<const CompressedMarker>&,
        const Align5::Options&,
        MemoryMapped::ByteAllocator&,
        Alignment&,
        AlignmentInfo&,
        bool debug);
}



class shasta::Align5::Options {
public:
    uint64_t deltaX;
    uint64_t deltaY;
    int64_t matchScore;
    int64_t mismatchScore;
    int64_t gapScore;
};



class shasta::Align5::MatrixEntry {
public:
    Coordinates xy;
    MatrixEntry() {}
    MatrixEntry(const Coordinates& xy) : xy(xy) {}
};



class shasta::Align5::Aligner {
public:

    using MarkerSequence = span<const CompressedMarker>;

    // The constructor does all the work.
    Aligner(
        const MarkerSequence&,
        const MarkerSequence&,
        const Options&,
        MemoryMapped::ByteAllocator&,
        Alignment&,
        AlignmentInfo&,
        bool debug);

private:

    // Number of markers (not features) in the two sequences being aligned.
    uint32_t nx;
    uint32_t ny;

    // Cell sizes in the X and Y direction.
    uint32_t deltaX;
    uint32_t deltaY;

    // For each sequence, vectors of pairs (markerId, ordinal)
    // sorted by marker id.
    vector< pair<KmerId, uint32_t> > sortedMarkers0;
    vector< pair<KmerId, uint32_t> > sortedMarkers1;
    static void sortMarkers(
        const MarkerSequence&,
        vector< pair<KmerId, uint32_t> >& sortedMarkers);

    // The alignment matrix, in a sparse representation organized by
    // cells in (X,Y) space.
    // For each iY, we store pairs(iX, xy) sorted by iX.
    // Even though this requires sorting, it is more efficient
    // than using a hash table, due to the better memory access pattern.
    using AlignmentMatrixEntry = pair<uint32_t, Coordinates>; // (iX, xy)
    using AlignmentMatrixAllocator = MemoryMapped::Allocator<AlignmentMatrixEntry>;
    using AlignmentMatrixEntryVector = vector<AlignmentMatrixEntry, AlignmentMatrixAllocator>; // For one iY
    using AlignmentMatrix = vector<AlignmentMatrixEntryVector>; // Indexed by iY.
    AlignmentMatrix alignmentMatrix;
    void createAlignmentMatrix();
    void writeAlignmentMatrixCsv(const string& fileName) const;
    void writeAlignmentMatrixPng(
        const string& fileName,
        uint32_t maxDistanceFromBoundary) const;
    void writeCheckerboard(
        PngImage&,
        uint32_t maxDistanceFromBoundary) const;



    // Cells in (X,Y) space.
    // Stored similarly to alignmentMatrix above: for each iY,
    // we store pairs (iX, Cell) sorted by iX.
    class Cell {
    public:
        uint8_t isNearLeftOrTop         : 1;
        uint8_t isNearRightOrBottom     : 1;
        uint8_t isForwardAccessible     : 1;
        uint8_t isBackwardAccessible    : 1;
        Cell()
        {
            *reinterpret_cast<uint8_t*>(this) = 0;
        }
        bool isActive() const
        {
            return (isForwardAccessible == 1) and (isBackwardAccessible == 1);
        }
    };
    static_assert(sizeof(Cell)==1, "Unexpected size of Align5::Aligner::Cell.");
    vector< vector< pair<uint32_t, Cell> > > cells;
    void createCells(
        uint32_t minEntryCountPerCell,
        uint32_t maxDistanceFromBoundary);
    void writeCellsCsv(const string& fileName) const;
    void writeCellsPng(const string& fileName) const;

    // Return the distance of a cell from the left boundary
    // of the alignment matrix, or 0 if the cell
    // is partially or entirely to the left of that boundary.
    uint32_t cellDistanceFromLeft(const Coordinates& iXY) const;

    // Return the distance of a cell from the right boundary
    // of the alignment matrix, or 0 if the cell
    // is partially or entirely to the right of that boundary.
    uint32_t cellDistanceFromRight(const Coordinates& iXY) const;

    // Return the distance of a cell from the top boundary
    // of the alignment matrix, or 0 if the cell
    // is partially or entirely above that boundary.
    uint32_t cellDistanceFromTop(const Coordinates& iXY) const;

    // Return the distance of a cell from the bottom boundary
    // of the alignment matrix, or 0 if the cell
    // is partially or entirely below that boundary.
    uint32_t cellDistanceFromBottom(const Coordinates& iXY) const;

    // Find a cell with given (iX,iY).
    Cell* findCell(const Coordinates& iXY);
    const Cell* findCell(const Coordinates& iXY) const;



    // Coordinate transformations.

    // Return (X,Y) given (x,y).
    Coordinates getXY(Coordinates xy) const;

    // Return (iX,iY) given (X,Y).
    Coordinates getCellIndexesFromXY(Coordinates XY) const
    {
        return Coordinates(
            XY.first  / deltaX,
            XY.second / deltaY
            );
    }

    // Return (iX,iY) given (x,y).
    Coordinates getCellIndexesFromxy(Coordinates xy) const
    {
        const Coordinates XY = getXY(xy);
        return getCellIndexesFromXY(XY);
    }


    // Convert an arbitrary (X,Y) to (x,y).
    // If the point is outside the alignment matrix,
    // we can end up with negative values.
    SignedCoordinates getxy(Coordinates XY) const;

    // Searches in cell space.
    void forwardSearch();
    void backwardSearch();



    // Group active cells in connected component.
    // Each connected component also generates a diagonal range.
    vector< vector<Coordinates> > activeCellsConnectedComponents;
    vector< pair<int32_t, int32_t > > diagonalRanges;
    void findActiveCellsConnectedComponents();



    MemoryMapped::ByteAllocator& byteAllocator;


};



#endif
