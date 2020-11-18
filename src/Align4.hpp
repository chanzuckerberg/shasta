#ifndef SHASTA_ALIGN4_HPP
#define SHASTA_ALIGN4_HPP


#include "hashArray.hpp"
#include "Marker.hpp"
#include "span.hpp"

#include "array.hpp"
#include <unordered_map>
#include "utility.hpp"



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
        bool debug,
        ostream& html);

    template<uint64_t m> void align4(
        const span<const CompressedMarker>&,
        const span<const CompressedMarker>&,
        const Align4Options&,
        Alignment&,
        AlignmentInfo&,
        bool debug,
        ostream& html);
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
        bool debug,
        ostream& html);

private:

    // A Feature is a sequence of m markers.
    using Feature = array<KmerId, m>;

    // A FeatureMap gives the ordinals where each Feature occurs in one of the
    // sequences being aligned.
    using FeatureMap = std::unordered_multimap<Feature, uint32_t, HashTuple<Feature> >;
    static void fillFeatureMap(const Sequence&, FeatureMap&);



    // The alignment matrix stores pairs (ordinal0, ordinal1) giving
    // starting ordinals for each common feature between the two sequences.
    // For efficient look ups, these pairs of ordinals are stored
    // in rectangular cells in diagonal coordinates (X, Y)
    // x = ordinal in sequence 0 (nx values starting at 0).
    // y = ordinal in sequence 1 (ny values starting at 0).
    // X = x + y
    // Y = y + (nx - 1 - x)
    // Note that with these definitions X and Y are never negative.
    // The cell size in the X direction is 2*maxSkip and
    // the cell size in the Y direction is maxDrift.
    // That way, when looking for successors of an alignment matrix element
    // in cell (iX,iY), we only have to look in 5 cells:
    // (iX, iY-1), (iX, iY+1), (iX+1, iY-1), (iX+1, iY, iX+1, iY+1).
    // The AlignmentMatrix is keyed by (iX, iY).
    using OrdinalPair = pair<uint32_t, uint32_t>;
    using Cell = pair<uint32_t, uint32_t>;
    using AlignmentMatrix = std::unordered_multimap<Cell, OrdinalPair, HashTuple<Cell> >;
    static void fillAlignmentMatrix(
        const FeatureMap& featureMap0,
        const Sequence& sequence1,
        uint64_t nx,
        uint64_t cellSizeX,
        uint64_t cellSizeY,
        AlignmentMatrix&
    );
    static void write(const AlignmentMatrix&, ostream& html);
};



#endif
