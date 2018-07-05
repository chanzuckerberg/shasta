#ifndef CZI_SHASTA_FIND_MARKER_ID_HPP
#define CZI_SHASTA_FIND_MARKER_ID_HPP

#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

#include "cstdint.hpp"
#include "tuple.hpp"
#include "utility.hpp"

namespace ChanZuckerberg {
    namespace shasta {

        // Given a global marker id in the global marker table,
        // return the corresponding OrientedReadId and ordinal.
        // This requires a binary search in the markers toc.
        inline pair<OrientedReadId, uint32_t> findMarkerId(
            MarkerId,
            const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers);

    }
}


inline std::pair<ChanZuckerberg::shasta::OrientedReadId, uint32_t>
    ChanZuckerberg::shasta::findMarkerId(
    MarkerId markerId,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers)
{
    OrientedReadId::Int orientedReadIdValue;
    uint32_t ordinal;
    tie(orientedReadIdValue, ordinal) = markers.find(markerId);
    return make_pair(OrientedReadId(orientedReadIdValue), ordinal);
}


#endif
