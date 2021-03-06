#ifndef SHASTA_INDUCED_ALIGNMENT_HPP
#define SHASTA_INDUCED_ALIGNMENT_HPP

// Shasta.
#include "MarkerGraph.hpp"

/*******************************************************************************

The marker graph induces an effective alignment between each pair
of oriented reads which can be obtained by following each of the oriented reads
in the marker graph. Aligned markers are those that are on the same vertex.

The induced alignment matrix of two oriented reads x and y
with nx and ny markers is an nx by ny matrix.
Element ij of the matrix is 1 if marker i< of x and marker <j of y
are on the same marker graph vertex and 0 otherwise.

*******************************************************************************/

// Shasta.
#include "MarkerGraph.hpp"

// Standard library.
#include "algorithm.hpp"
#include "array.hpp"
#include "string.hpp"
#include "tuple.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    class InducedAlignment;
    class InducedAlignmentData;
}



class shasta::InducedAlignmentData {
public:

    MarkerGraph::VertexId vertexId;

    // The marker ordinals in the two reads.
    uint32_t ordinal0;
    uint32_t ordinal1;

    // The compressed ordinals. These only count markers
    // that are associated with a marker graph vertex.
    uint32_t compressedOrdinal0 = std::numeric_limits<uint32_t>::max();
    uint32_t compressedOrdinal1 = std::numeric_limits<uint32_t>::max();

    InducedAlignmentData(
        MarkerGraph::VertexId vertexId,
        uint32_t ordinal0,
        uint32_t ordinal1
        ) :
        vertexId(vertexId),
        ordinal0(ordinal0),
        ordinal1(ordinal1)
        {}

    // For convenience order by the ordinals.
    // But this ordering does not have an meaning.
    bool operator<(const InducedAlignmentData& that) const
    {
        return tie(ordinal0, ordinal1) < tie(that.ordinal0, that.ordinal1);
    }
};



class shasta::InducedAlignment {
public:

    // A vector defining this induced alignment.
    vector<InducedAlignmentData> data;

    // The number marker of markers associated with a
    /// marker graph vertex, for each of the oriented reads
    // involved in this induced alignment.
    array<uint32_t, 2> compressedMarkerCount;

    void sort()
    {
        std::sort(data.begin(), data.end());
    }

    void writePngImage(
        uint32_t markerCount0,
        uint32_t markerCount1,
        bool useCompressedOrdinals,
        const string& fileName) const;

};



#endif
