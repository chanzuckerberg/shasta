#ifndef SHASTA_ASSEMBLE_MARKER_GRAPH_PATH_HPP
#define SHASTA_ASSEMBLE_MARKER_GRAPH_PATH_HPP

#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "span.hpp"

namespace shasta {

    class AssembledSegment;

    void assembleMarkerGraphPath(
        uint64_t readRepresentation,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&,
        const span<const MarkerGraph::EdgeId>& markerGraphPath,
        bool storeCoverageData,
        AssembledSegment& assembledSegment);

}



#endif
