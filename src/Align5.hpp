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

We also consider Feature's which are sequences of m markers
in each of the two sequences, where m is the template parameter of
class Aligner.

So for example consider an input sequence
consisting of the following marker KmerId's:
45 58 106 17
If m=2, the sequence of Feature's representing this sequence is
(45,58) (58,106), (106,17).

The two sequences of Feature's corresponding to
markerSequence0 and markerSequence1
are featureSequence0 and featureSequence1.

Note that the sequences of Feature's are shorter (by m-1)
than the original marker sequences.

The alignment matrix in feature space is sparse because of the
large alphabet. For example, with default options there are
about 8000 marker KmerId's and therefore about 64000000
distinct features for m=2.

The coordinates in the alignment matrix in marker or feature space
are x and y, with x
represented along the horizontal axis and increasing toward the right,
and y represented along the vertical axis and increasing toward
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
So the total number of disntict values of X and Y is nx + ny - 1.

X is a coordinate along the diagonal of the alignment matrix,
and Y is orthogonal to it and identifies the diagonal.
In (X,Y) coordinates the alignment matrix is a subset of
the square of size nx + ny -1. The alignment matrix is
rotated by 45 degrees relative to this square.

Consider two successive points (x0,y0) and (x1, y1) in
a hypothetical alignment in Feature space,
and corresponding (X0,Y0) and (X1,Y1).
For a good alignment, the "skip" along the diagonal, |X1-X0|,
must be sufficiently small, and the "diagonal drift",
|Y1-Y0| must also be sufficiently small.
We require
|X1-X0| <= deltaX
|Y1-Y0| <= deltaY
Note that deltaX is similar in spirit to Shasta maxSkip
when working in marker space (as opposed to Feature space).
Similarly, deltaY is similar to Shasta maxDrift.

We construct a sparse representation of this sparse
alignment matrix with a special structure that, given an
alignment matrix entry at (x0,y0), allows efficient
look up of entries (x1,y1) which satisfy the
deltaX and deltaY criteria.

To achieve this, we use rectangular cells in (X,Y) space of size
(deltaX, deltaY). If (iX,iY) are the indices of a cell
containing (x0,y0), all possible (x1,y1) must
be in in the 6 cells
(iX-1, iY  ) (iX,iY  ) (iX+1, iY  )
(iX-1, iY+1) (iX,iY+1) (iX+1, iY+1)


*******************************************************************************/

#include "hashArray.hpp"
#include "Marker.hpp"
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
        template<uint64_t m> class Aligner;
        class Options;

    }

    void align5(
        const span<const CompressedMarker>&,
        const span<const CompressedMarker>&,
        const Align5::Options&,
        Alignment&,
        AlignmentInfo&,
        bool debug);

    template<uint64_t m> void align5(
        const span<const CompressedMarker>&,
        const span<const CompressedMarker>&,
        const Align5::Options&,
        Alignment&,
        AlignmentInfo&,
        bool debug);

}



class shasta::Align5::Options {
public:
    uint64_t m;
    uint64_t deltaX;
    uint64_t deltaY;
    int64_t matchScore;
    int64_t mismatchScore;
    int64_t gapScore;
};



template<uint64_t m> class shasta::Align5::Aligner {
public:

    using MarkerSequence = span<const CompressedMarker>;

    // The constructor does all the work.
    Aligner(
        const MarkerSequence&,
        const MarkerSequence&,
        const Options&,
        Alignment&,
        AlignmentInfo&,
        bool debug);

private:

    // Number of markers (not features in the two sequences being aligned.
    int32_t nx;
    int32_t ny;

    // Cell sizes in the X and Y direction.
    int32_t deltaX;
    int32_t deltaY;

    // For each sequence, vectors of pairs (Feature, ordinal)
    // sorted by feature.
    using Feature = array<KmerId, m>;
    vector< pair<Feature, uint32_t> > sortedFeatures0;
    vector< pair<Feature, uint32_t> > sortedFeatures1;
    static void sortFeatures(
        const MarkerSequence&,
        vector< pair<Feature, uint32_t> >&);

    // For each sequence, vectors of pairs (markerId, ordinal)
    // sorted by marker id.
    vector< pair<KmerId, uint32_t> > sortedMarkers0;
    vector< pair<KmerId, uint32_t> > sortedMarkers1;
    static void sortMarkers(
        const MarkerSequence&,
        vector< pair<KmerId, uint32_t> >& sortedMarkers);

    // Write the alignment matrix in marker space or feature space
    // to a png image.
    void writeAlignmentMatrixInMarkerSpace(const string& fileName) const;
    void writeAlignmentMatrixInFeatureSpace(const string& fileName) const;
    void writeCheckerboard(PngImage&) const;

    // This is used to store (x,y), (X,Y), or (iX, iY).
    using Coordinates = pair<uint32_t, uint32_t>;

    // An element of the alignment matrix in feature space.
    class AlignmentMatrixEntry {
    public:
        Coordinates xy;
        AlignmentMatrixEntry(const Coordinates& xy) : xy(xy) {}
    };

    // Cells in (X,Y) space.
    class Cell {
    public:
        vector<AlignmentMatrixEntry> alignmentMatrixEntries;
    };
    std::unordered_map<Coordinates, Cell, HashTuple<Coordinates> > cells;
    void createCells();
    void writeCells(const string& fileName) const;

    // Given ordinals x, y, return coordinates iX, iY of the containing cell.
    pair<uint32_t, uint32_t> getCellIndexes(uint32_t x, uint32_t y) const;
};



#endif
