#ifndef CZI_NANOPORE2_FIND_MARKER_ID_HPP
#define CZI_NANOPORE2_FIND_MARKER_ID_HPP

#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

#include "cstdint.hpp"
#include "tuple.hpp"
#include "utility.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {

        // Given a global marker id in the global marker table,
        // return the corresponding OrientedReadId and ordinal.
        // This requires a binary search in the markers toc.
        inline pair<OrientedReadId, uint32_t> findMarkerId(
            MarkerId,
            const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers);

    }
}


inline std::pair<ChanZuckerberg::Nanopore2::OrientedReadId, uint32_t>
    ChanZuckerberg::Nanopore2::findMarkerId(
    MarkerId markerId,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers)
{
    OrientedReadId::Int orientedReadIdValue;
    uint32_t ordinal;
    tie(orientedReadIdValue, ordinal) = markers.find(markerId);
    return make_pair(OrientedReadId(orientedReadIdValue), ordinal);
}


#endif
