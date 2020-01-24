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
    class InducedAlignmentCriteria;
}



// Criteria used to evaluate an InducedAlignment.
class shasta::InducedAlignmentCriteria {
public:
    uint32_t maxOffsetSigma;
    uint32_t maxTrim;
    uint32_t maxSkip;
};



class shasta::InducedAlignmentData {
public:

    MarkerGraph::VertexId vertexId;

    // The marker ordinals in the two reads.
    uint32_t ordinal0;
    uint32_t ordinal1;

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

    void sort()
    {
        std::sort(data.begin(), data.end());
    }

    void writePngImage(
        uint32_t markerCount0,
        uint32_t markerCount1,
        const string& fileName) const;

    // Evaluate the quality of an induced alignment.
    // Returns true if the induced alignment satisfies the specified criteria.
    bool evaluate(
        uint32_t markerCount0,
        uint32_t markerCount1,
        const InducedAlignmentCriteria&) const;
    bool evaluate(
        uint32_t markerCount0,
        uint32_t markerCount1,
        uint32_t leftTrim0,
        uint32_t rightTrim0,
        uint32_t leftTrim1,
        uint32_t rightTrim1,
        const InducedAlignmentCriteria&) const;
};



#endif
