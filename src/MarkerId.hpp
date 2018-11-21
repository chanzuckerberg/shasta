#ifndef CZI_SHASTA_MARKER_ID_HPP
#define CZI_SHASTA_MARKER_ID_HPP

#include "cstdint.hpp"
#include "Uint.hpp"

namespace ChanZuckerberg {
    namespace shasta {

        // Type used to globally identify a marker on an oriented read.
        // This is the global index of the marker in Assembler::markers.
        // For a human assembly with coverage 40X the total number
        // of markers is more than 20 billions (counting both strands),
        // so this needs to be uint64_t. There could, however, be situations
        // where uint32_t is sufficient.
        using MarkerId = uint64_t;

        // Types used to identify a vertex and edge of the global marker graph.
        using GlobalMarkerGraphVertexId = MarkerId;
        const GlobalMarkerGraphVertexId invalidGlobalMarkerGraphVertexId =
            std::numeric_limits<GlobalMarkerGraphVertexId>::max();
        using GlobalMarkerGraphEdgeId = GlobalMarkerGraphVertexId;
        const GlobalMarkerGraphEdgeId invalidGlobalMarkerGraphEdgeId =
            invalidGlobalMarkerGraphVertexId;


        // To save memory, store it using 5 bytes.
        // This allows for up to 2^40 = 1 Ti markers (both strands).
        // A human size run with 40x coverage and 10% markers
        // has around 25 G markers (both strands).
        using CompressedGlobalMarkerGraphVertexId = Uint40;
        const CompressedGlobalMarkerGraphVertexId
            invalidCompressedGlobalMarkerGraphVertexId = invalidGlobalMarkerGraphVertexId;

    }
}

#endif
