// Shasta.
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "InducedAlignment.hpp"
#include "orderPairs.hpp"
using namespace shasta;



// Compute an alignment between two oriented reads
// induced by the marker graph. See InducedAlignment.hpp for more
// information.
void Assembler::computeInducedAlignment(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    InducedAlignment& inducedAlignment
)
{
    // Get the Marker graph vertices for the two oriented reads.
    using VertexId = MarkerGraph::VertexId;
    vector< pair<uint32_t, VertexId> > vertices0;
    vector< pair<uint32_t, VertexId> > vertices1;
    getMarkerGraphVertices(orientedReadId0, vertices0);
    getMarkerGraphVertices(orientedReadId1, vertices1);

    // Sort them by vertex id.
    sort(vertices0.begin(), vertices0.end(),
        OrderPairsBySecondOnly<uint32_t, VertexId>());
    sort(vertices1.begin(), vertices1.end(),
        OrderPairsBySecondOnly<uint32_t, VertexId>());

    // Some iterators we need below.
    using Iterator =  vector< pair<uint32_t, VertexId> >::iterator;
    const Iterator begin0 = vertices0.begin();
    const Iterator begin1 = vertices1.begin();
    const Iterator end0 = vertices0.end();
    const Iterator end1 = vertices1.end();



    // Gather the vertices of the induced alignments.
    inducedAlignment.data.clear();
    Iterator it0 = begin0;
    Iterator it1 = begin1;
    while(it0 not_eq end0 and it1 not_eq end1) {
        const VertexId vertexId0 = it0->second;
        const VertexId vertexId1 = it1->second;
        if(vertexId0 < vertexId1) {
            ++it0;
            continue;
        }
        if(vertexId1 < vertexId0) {
            ++it1;
            continue;
        }

        // If getting here, we found a common vertex.
        SHASTA_ASSERT(vertexId0 == vertexId1);
        const VertexId vertexId = vertexId0;


        // See if there are more markers from these reads on these vertices.
        Iterator it0Begin = it0;
        Iterator it1Begin = it1;
        Iterator it0End = it0Begin;
        Iterator it1End = it1Begin;
        while(it0End!=end0 and it0End->second==vertexId) {
            ++it0End;
        }
        while(it1End!=end1 and it1End->second==vertexId) {
            ++it1End;
        }

        // When generating marker graph vertices, we don't allow two
        // markers from the same oriented read on a vertex.
        // So both streaks should have size 1.
        // However the code below does not make this assumption.
        SHASTA_ASSERT(it0End - it0Begin == 1);
        SHASTA_ASSERT(it1End - it1Begin == 1);

        // Loop over pairs in the streaks.
        for(Iterator jt0=it0Begin; jt0!=it0End; ++jt0) {
            for(Iterator jt1=it1Begin; jt1!=it1End; ++jt1) {
                inducedAlignment.data.push_back(
                    InducedAlignmentData(vertexId, jt0->first, jt1->first));
            }
        }

       // Position at the end of these streaks to continue processing.
        it0 = it0End;
        it1 = it1End;
    }

    inducedAlignment.sort();
}



// Compute induced alignments between an oriented read orientedReadId0
// and the oriented reads stored sorted in orientedReadIds1.
void Assembler::computeInducedAlignments(
    OrientedReadId orientedReadId0,
    const vector<OrientedReadId>& orientedReadIds1,
    vector<InducedAlignment>& inducedAlignments
)
{
    // Initialize one induced alignment for each
    // OrientedReadId in orientedReadIds1.
    inducedAlignments.clear();
    inducedAlignments.resize(orientedReadIds1.size());

    // Fast exit if orientedReadIds1 is empty.
    if(orientedReadIds1.empty()) {
        return;
    }

    // Check that orientedReadIds1 is sorted.
    SHASTA_ASSERT(std::is_sorted(orientedReadIds1.begin(), orientedReadIds1.end()));



    // Main loop over markers in orientedReadId0.
    const MarkerId firstMarkerId = markers.begin(orientedReadId0.getValue()) - markers.begin();
    const uint32_t markerCount = uint32_t(markers.size(orientedReadId0.getValue()));
    for(uint32_t ordinal0=0; ordinal0<markerCount; ordinal0++) {
        const MarkerId markerId0 = firstMarkerId + ordinal0;

        // Find the vertex that this marker is on.
        const MarkerGraph::CompressedVertexId compressedVertexId =
            markerGraph.vertexTable[markerId0];

        // If this marker is not on a marker graph vertex, skip.
        if(compressedVertexId == MarkerGraph::invalidCompressedVertexId) {
            continue;
        }

        // Loop over all markers on this vertex.
        const span<MarkerId> vertexMarkers =
            markerGraph.getVertexMarkerIds(compressedVertexId);
        for(const MarkerId markerId1: vertexMarkers) {

            // Skip the marker that we started from.
            if(markerId1 == markerId0) {
                continue;
            }

            // Find the oriented read on this marker.
            OrientedReadId orientedReadId1;
            uint32_t ordinal1;
            tie(orientedReadId1, ordinal1) = findMarkerId(markerId1);

            // Look for orientedReadId1 in the orientedReadIds1 vector.
            const vector<OrientedReadId>::const_iterator it =
                std::lower_bound(
                    orientedReadIds1.begin(),
                    orientedReadIds1.end(),
                    orientedReadId1);
            if((it == orientedReadIds1.end()) or (*it != orientedReadId1)) {
                continue;
            }

            // Add this marker to the corresponding induced alignment.
            const uint64_t index1 = it -orientedReadIds1.begin();
            const MarkerGraph::VertexId vertexId = compressedVertexId;
            inducedAlignments[index1].data.push_back(
                InducedAlignmentData(vertexId, ordinal0, ordinal1));
        }
    }

    // Sort the induced alignments.
    for(InducedAlignment& inducedAlignment: inducedAlignments) {
        inducedAlignment.sort();
    }

    // Fill in the compressed ordinals.
    SHASTA_ASSERT(orientedReadIds1.size() == inducedAlignments.size());
    for(uint64_t i=0; i<orientedReadIds1.size(); i++){
        const OrientedReadId orientedReadId1 = orientedReadIds1[i];
        InducedAlignment& inducedAlignment = inducedAlignments[i];
        fillCompressedOrdinals(orientedReadId0, orientedReadId1, inducedAlignment);
    }

}



// Fill in compressed ordinals of an InducedAlignment.
// Compressed ordinals are marker ordinals in which
// only markers associated with a marker graph vertex are counted.
void Assembler::fillCompressedOrdinals(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    InducedAlignment& inducedAlignment)
{
    array<OrientedReadId, 2> orientedReadIds = {orientedReadId0, orientedReadId1};



    // Create vectors to convert ordinals to compressed ordinals.
    // Indexed by the ordinal.
    // This also fills in the compressedMarkerCount in the induced alignment.
    array< vector<uint32_t>, 2> ordinalTable;
    for(uint64_t i=0; i<2; i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];

        const MarkerId firstMarkerId = markers.begin(orientedReadId.getValue()) - markers.begin();
        const uint32_t markerCount = uint32_t(markers.size(orientedReadId.getValue()));
        ordinalTable[i].resize(markerCount);

        uint32_t compressedOrdinal = 0;
        for(uint32_t ordinal=0; ordinal<markerCount; ordinal++) {
            ordinalTable[i][ordinal] = compressedOrdinal;
            const MarkerId markerId = firstMarkerId + ordinal;

            // Find the vertex that this marker is on.
            const MarkerGraph::CompressedVertexId compressedVertexId =
                markerGraph.vertexTable[markerId];

            // If this marker is on a marker graph vertex, increment
            // the compressed ordinal.
            if(compressedVertexId != MarkerGraph::invalidCompressedVertexId) {
                ++compressedOrdinal;
            }
        }

        // Store the number of compressed markers on this oriented read.
        inducedAlignment.compressedMarkerCount[i] = compressedOrdinal;
    }


    // Now we can fill in the compressed ordinal.
    for(InducedAlignmentData& data: inducedAlignment.data) {
        data.compressedOrdinal0 = ordinalTable[0][data.ordinal0];
        data.compressedOrdinal1 = ordinalTable[1][data.ordinal1];
    }

}
