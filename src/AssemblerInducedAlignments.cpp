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



// Find all pairs of incompatible reads that involve a given read.
// A pair of reads is incompatible if it has a "bad" induced alignment.
// See InducedAlignment.hpp for more information.
// Python-callable overload.
vector<OrientedReadPair> Assembler::findIncompatibleReadPairs(
    ReadId readId0,

    // If true, only consider ReadId's readid1<readId0.
    bool onlyConsiderLowerReadIds,

    // If true, skip pairs that are in the read graph.
    // Those are already known to have a good induced alignment
    // by construction.
    bool skipReadGraphEdges)
{
    vector<OrientedReadPair> incompatiblePairs;
    findIncompatibleReadPairs(
        readId0,
        onlyConsiderLowerReadIds,
        skipReadGraphEdges,
        incompatiblePairs);
    return  incompatiblePairs;
}



// Find all pairs of incompatible reads that involve a given read.
// A pair of reads is incompatible if it has a "bad" induced alignment.
// See InducedAlignment.hpp for more information.
void Assembler::findIncompatibleReadPairs(
    ReadId readId0,

    // If true, only consider ReadId's readId1<readId0.
    bool onlyConsiderLowerReadIds,

    // If true, skip pairs that are in the read graph.
    // Those are already known to have a good induced alignment
    // by construction.
    bool skipReadGraphEdges,

    // The incompatible pairs found.
    vector<OrientedReadPair>& incompatiblePairs)
{
    const bool debug = true;

    // Criteria used to evaluate induced alignments.
    // Get thsi values from comand line options when the code stabilizes.
    InducedAlignmentCriteria inducedAlignmentCriteria;
    inducedAlignmentCriteria.maxOffsetSigma = 50;
    inducedAlignmentCriteria.maxTrim = 100;
    inducedAlignmentCriteria.maxSkip = 100;

    // Check that we have what we need.
    // The code as written only supports the directed read graph.
    SHASTA_ASSERT(directedReadGraph.edges.isOpen);
    SHASTA_ASSERT(directedReadGraph.edgesBySource.isOpen());
    SHASTA_ASSERT(directedReadGraph.edgesByTarget.isOpen());


    // We need to find oriented reads that share at least one marker graph
    // vertex with this read (on the positive strand).
    // We will store them in this vector.
    vector<OrientedReadId> incompatibleCandidates;



    // To do this, we loop over all markers of this
    // read (on the positive strand)
    const OrientedReadId orientedReadId0(readId0, 0);
    const MarkerId firstMarkerId = markers.begin(orientedReadId0.getValue()) - markers.begin();
    const uint32_t markerCount = uint32_t(markers.size(orientedReadId0.getValue()));
    for(uint32_t ordinal=0; ordinal<markerCount; ordinal++) {
        const MarkerId markerId0 = firstMarkerId + ordinal;

        // Find the vertex that this marker is on.
        const MarkerGraph::CompressedVertexId compressedVertexId =
            markerGraph.vertexTable[markerId0];

        // If this marker is not on a marker graph vertex, skip.
        if(compressedVertexId == MarkerGraph::invalidCompressedVertexId) {
            continue;
        }

        // Loop over all markers on this vertex,
        // except the one on orientedReadId0 that we started from.
        const MemoryAsContainer<MarkerId> vertexMarkers =
            markerGraph.vertices[compressedVertexId];
        for(const MarkerId markerId1: vertexMarkers) {
            if(markerId1 == markerId0) {
                continue;
            }

            // Find the oriented read on this marker.
            OrientedReadId orientedReadId1;
            uint32_t ordinal1;
            tie(orientedReadId1, ordinal1) = findMarkerId(markerId1);
            if(orientedReadId1.getReadId() == readId0) {
                continue;
            }



            // See if this pair should be skipped because it corresponds
            // to a read graph edge.
            if(skipReadGraphEdges) {
                const DirectedReadGraph::VertexId v0 = orientedReadId0.getValue();
                const DirectedReadGraph::VertexId v1 = orientedReadId1.getValue();
                const bool forwardExists = directedReadGraph.findEdge(v0, v1)
                    != DirectedReadGraph::invalidEdgeId;
                const bool backwardExists = directedReadGraph.findEdge(v1, v0)
                    != DirectedReadGraph::invalidEdgeId;

                if(forwardExists or backwardExists) {
                    continue;
                }
            }




            // Add this oriented read to our incompatible candidates.
            if(onlyConsiderLowerReadIds) {

                // Only if readId1 < readId0.
                if(orientedReadId1.getReadId() < orientedReadId0.getReadId()) {
                    incompatibleCandidates.push_back(orientedReadId1);
                }

            } else {

                // Unconditionally.
                incompatibleCandidates.push_back(orientedReadId1);
            }
        }
    }
    deduplicate(incompatibleCandidates);
    if(debug) {
        cout << "Found " << incompatibleCandidates.size() <<
            " incompatible candidates." << endl;
    }




    // Compute the induced alignment with each of these incompatible candidates.
    incompatiblePairs.clear();
    InducedAlignment inducedAlignment;
    for(const OrientedReadId orientedReadId1: incompatibleCandidates) {
        cout << "Checking induced alignment of " <<
            orientedReadId0 << " and " << orientedReadId1 << endl;
        computeInducedAlignment(
            orientedReadId0,
            orientedReadId1,
            inducedAlignment);

        const uint32_t markerCount0 = uint32_t(markers.size(orientedReadId0.getValue()));
        const uint32_t markerCount1 = uint32_t(markers.size(orientedReadId1.getValue()));
        if(not inducedAlignment.evaluate(markerCount0, markerCount1, inducedAlignmentCriteria)) {
            incompatiblePairs.push_back(OrientedReadPair(
                readId0,
                orientedReadId1.getReadId(),
                orientedReadId1.getStrand() == 0));
        }

    }
}


