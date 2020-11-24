#ifndef SHASTA_ALIGN4_HPP
#define SHASTA_ALIGN4_HPP

/*******************************************************************************

Marker alignment of two sequences, each defined as a sequence of
marker KmerId's.

We work with Feature's defined as a sequence of m markers,
where m is the template parameter of class Align4.

Each of the two input marker sequences sequence0 and sequence1
is treated as a sequence of Features's.
So for example consider an input sequence
consisting of the following marker KmerId's:
45 58 106 17
If m=2, the sequence of Feature's representing this sequence is
(45,58) (58,106), (106,17).

We call x the index or position in the first sequence sequence0
and y the index or position in the second sequence sequence1.
sequence[x] is the Feature at position x of sequence0
and sequence1[y] is the feature at position y of sequence1.
In other Shasta code, x and y are often called ordinals or marker
ordinals - but here they refer to Feature's, not markers.
The number of features in sequence0 is nx
and the number of features in sequence1 is ny.
For any two positions x and y the following hold:
0 <= x <= nx-1
0 <= y <= ny-1

Note that the sequences of Feature's are shorter (by m-1)
than the original marker sequences.

The alignment matrix in feature space is sparse because of the
large alphabet. For example, with default options there are
about 8000 marker KmerId's and therefore about 64000000
distinct features.

The coordinates in the alignment matrix are x and y, with x
represented along the horizontal axis and increasing toward the right,
and y represented along the vertical axis and increasing toward
the bottom. The alignment matrix element at position (x,y)
exists if sequence0[x]=sequence1[y].

We also use coordinates X and Y defined as:
X = x + y
Y = y + (nx - 1 - x)

It can be verified that:
0 <= X <= nx + ny - 2
0 <= Y <= nx + ny - 2

X is a coordinate along the diagonal of the alignment matrix,
and Y is orthogonal to it and identifies the diagonal.
In (X,Y) coordinates the alignment matrix is a subset of
the square of size nx + ny -1. The alignment matrix is
rotated by 45 degrees relative to this square.

Consider two successive points (x0,y0) and (x1, y1) in
a hypothetical alignment, and corresponding (X0,Y0) and (X1,Y1).
For a good alignment, the "skip" along the diagonal, |X1-X0|,
must be sufficiently small, and the "diagonal drift",
|Y1-Y0| must also be sufficiently small.
We require
|X1-X0| <= deltaX
|Y1-Y0| <= deltaY
Note that deltaX is similar in spirit to Shasta maxSkip
when working in marker space (9as opposed to Feature space).
Similarly, deltaY is similaro to Shasta maxDrift.

We construct a sparse representation of this sparse
alignment matrix with a special structure that, given an
alignment matrix entry at (x0,y0), allows efficient
look up of entries (x1,y1) which satisfy the
deltaX and deltaY criteria.

To achieve this, we use rectangular cells in (X,Y) space of size
(deltaX, deltaY). If (iX,iY) are the indices of a cell
containing (x0,y0), to locate all possible (x1,y1) we must
only look in the same cell (iX, iY) plus all of its 8
immediate neighbors (including diagonal neighbors).

The AlignmentMatrix is stored as an unordered_multimap with key (iX, iY)
and values consisting of Feature's.

To construct candidate alignments, we consider each element
in the alignment matrix as a vertex in a directed graph.
There is an edge between the vertices at positions (x0,y0) and (x1,y1)
if x1>=x0 and y1>=y0 and the above conditions on deltaX, deltaY
are satisfied.

Each candidate alignment must have at least one vertex near the top or left
and one vertex near the bottom or right of the alignment matrix.

To avoid doing a BFS on the entire alignment matrix, we do a BFS
starting only with vertices near the top or left.
We only keep connected components that also contain at
least one vertex neat the bottom or right.

For performance, we don't explicitly construct edges of the graph.
We find neighbors as needed.

*******************************************************************************/



#include "hashArray.hpp"
#include "Marker.hpp"
#include "span.hpp"

#include <boost/graph/adjacency_list.hpp>

#include "array.hpp"
#include <unordered_map>
#include <set>
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {
    template<uint64_t m> class Align4;
    class Align4Options;
    class Alignment;
    class AlignmentInfo;

    void align4(
        const span<const CompressedMarker>&,
        const span<const CompressedMarker>&,
        const Align4Options&,
        Alignment&,
        AlignmentInfo&,
        bool debug);

    template<uint64_t m> void align4(
        const span<const CompressedMarker>&,
        const span<const CompressedMarker>&,
        const Align4Options&,
        Alignment&,
        AlignmentInfo&,
        bool debug);
}



class shasta::Align4Options {
public:
    uint64_t m;
    uint64_t deltaX;
    uint64_t deltaY;
    int64_t matchScore;
    int64_t mismatchScore;
    int64_t gapScore;
};



template<uint64_t m> class shasta::Align4 {
public:

    using Sequence = span<const CompressedMarker>;

    // The constructor does all the work.
    Align4(
        const Sequence&,
        const Sequence&,
        const Align4Options&,
        Alignment&,
        AlignmentInfo&,
        bool debug);

private:

    // A Feature is a sequence of m markers.
    using Feature = array<KmerId, m>;

    // A FeatureMap gives the x (or y) where each Feature occurs in one of the
    // sequences being aligned.
    using FeatureMap = std::unordered_multimap<Feature, uint32_t, HashTuple<Feature> >;
    static void fillFeatureMap(const Sequence&, FeatureMap&);

    // Cell sizes in the X and Y direction.
    int32_t deltaX;
    int32_t deltaY;



    // An entry of the alignment matrix.
    using Coordinates = pair<int32_t, int32_t>;
    class AlignmentMatrixEntry {
    public:

        // The vertex_descriptor in the Graph defined below.
        // This is only valid after createGraph was called.
        // This happens after unreachable vertices are removed.
        uint64_t v = std::numeric_limits<uint64_t>::max();
        bool hasValidVertex() const {
            return v != std::numeric_limits<uint64_t>::max();
        }

        Coordinates xy;
        Coordinates XY;

        // Flags for elements near the boundary of the alignment matrix.

        // Set if x < deltaX.
        uint8_t isNearLeft : 1;

        // Set if y < deltaX.
        uint8_t isNearTop : 1;

        // Set if nx-1-x < deltaX.
        uint8_t isNearRight : 1;

        // Set if ny-1-y < deltaX.
        uint8_t isNearBottom : 1;

        // Set if there is a forward path to this feature starting
        // near the top or left of the alignment matrix.
        uint8_t isForwardReachableFromTopOrLeft : 1;

        // Set if there is a backward path to this feature starting
        // near the bottom or right of the alignment matrix.
        uint8_t isBackwardReachableFromBottomOrRight : 1;

        // Flag used during the BFS.
        bool wasDiscovered = false;

        bool isNeighbor(const AlignmentMatrixEntry&,
            int32_t deltaX,
            int32_t deltaY) const;

        // Return true if (*this) is a child of that.
        bool isChild(
            const AlignmentMatrixEntry& that,
            int32_t deltaX,
            int32_t deltaY) const;

        // Return true if (*this) is a parent of that.
        bool isParent(
            const AlignmentMatrixEntry& that,
            int32_t deltaX,
            int32_t deltaY) const;
    };


    // In the alignment matrix, the key has the (iX,iY) cell coordinates.
    using AlignmentMatrix = std::unordered_multimap<Coordinates, AlignmentMatrixEntry, HashTuple<Coordinates> >;
    AlignmentMatrix alignmentMatrix;
    void fillAlignmentMatrix(
        const FeatureMap& featureMap0,
        const Sequence& sequence1,
        int32_t nx,
        int32_t ny);
    void writeMatrixCsv(const string& fileName);
    void writeMatrixPng(
        uint32_t nx, uint32_t ny,
        const string& fileName);
    void clearDiscoveredFlags();

    // Compute the reachability flags in the alignment matrix.
    void computeReachability();

    // Remove alignment matrix entries that are not reachable in both directions.
    void removeUnreachable();

    // Find neighbors and flag them as discovered.
    void findAndFlagUndiscoveredNeighbors(
        typename AlignmentMatrix::iterator,
        vector<typename AlignmentMatrix::iterator>&
    );
    void findChildren(
        typename AlignmentMatrix::iterator,
        vector<typename AlignmentMatrix::iterator>&
    );
    void findAndFlagUndiscoveredChildren(
        typename AlignmentMatrix::iterator,
        vector<typename AlignmentMatrix::iterator>&
    );
    void findAndFlagUndiscoveredParents(
        typename AlignmentMatrix::iterator,
        vector<typename AlignmentMatrix::iterator>&
    );

    // Order AlignmentMatrix iterators by increasing X, then Y.
    class OrderByIncreasingXThenY {
    public:
        bool operator() (
            const typename AlignmentMatrix::iterator& it0,
            const typename AlignmentMatrix::iterator& it1) const
        {
            const AlignmentMatrixEntry& entry0 = it0->second;
            const AlignmentMatrixEntry& entry1 = it1->second;
            return entry0.XY < entry1.XY;
        }
    };



    // After unreachable vertices are removed from the alignment matrix,
    // we create the equivalent Boost graph.
    // Each vertex contains an AlignmentMatrixEntry.
    class AlignmentMatrixEntry;
    using Graph = boost::adjacency_list<
        boost::vecS,
        boost::vecS,
        boost::bidirectionalS,
        typename AlignmentMatrix::iterator,
        boost::property<boost::edge_weight_t, uint64_t> >;
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;
    Graph graph;
    void createGraph();
    void findShortestPaths(bool debug);
};



#endif
