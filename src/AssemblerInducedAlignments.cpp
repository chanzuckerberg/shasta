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
            markerGraph.vertices[compressedVertexId];
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
    const bool debug = false;

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
        const span<MarkerId> vertexMarkers =
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
        /*
        cout << "Checking induced alignment of " <<
            orientedReadId0 << " and " << orientedReadId1 << endl;
        */
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



// Find all incompatible read pairs.
void Assembler::findAllIncompatibleReadPairs()
{
    vector<OrientedReadPair> allIncompatiblePairs;
    vector<OrientedReadPair> incompatiblePairs;

    for(ReadId readId=0; readId<readCount(); readId++) {
        findIncompatibleReadPairs(readId, true, true, incompatiblePairs);
        copy(incompatiblePairs.begin(), incompatiblePairs.end(),
            back_inserter(allIncompatiblePairs));
    }
    cout << "Found " << allIncompatiblePairs.size() << " incompatible pairs for " <<
        readCount() << " reads." << endl;

    ofstream csv("IncompatiblePairs.csv");
    for(const OrientedReadPair& p: allIncompatiblePairs) {
        csv << p.readIds[0] << ",";
        csv << p.readIds[1] << ",";
        csv << int(p.isSameStrand) << "\n";
    }

    ofstream dot("IncompatiblePairs.dot");
    dot << "graph G{\n";
    for(const OrientedReadPair& p: allIncompatiblePairs) {
        OrientedReadId orientedReadId0(p.readIds[0], 0);
        OrientedReadId orientedReadId1(p.readIds[1], p.isSameStrand ? 0 : 1);
        dot << "\"" << orientedReadId0 << "\"--\"" << orientedReadId1 << "\";\n";
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        dot << "\"" << orientedReadId0 << "\"--\"" << orientedReadId1 << "\";\n";
    }
    dot << "}\n";
}




// Evaluate an induced alignment.
// Contrary to InducedAlignment::evaluate, this takes into account
// markers that don't correspond to a marker graph vertex.
bool Assembler::evaluateInducedAlignment(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const InducedAlignment& inducedAlignment,
    const InducedAlignmentCriteria& criteria,
    vector<bool>& work0,
    vector<bool>& work1)
{
    // Sanity check.
    SHASTA_ASSERT(not inducedAlignment.data.empty());

    // Compute the average and standard deviation of the offset
    // between the first and second ordinal.
    int64_t sum1 = 0;
    int64_t sum2 = 0;
    for(const auto& d: inducedAlignment.data) {
        const int64_t offset = int64_t(d.ordinal0) - int64_t(d.ordinal1);
        sum1 += offset;
        sum2 += offset * offset;
    }
    const double n = double(inducedAlignment.data.size());
    const double offset = double(sum1) / n;
    const double sigma = (inducedAlignment.data.size()==1) ? 0. :
        sqrt((double(sum2) - n*offset*offset) / (n-1.));

    // If the standard deviation is too large, return false.
    if(uint32_t(sigma) > criteria.maxOffsetSigma) {
        return false;
    }

    // Access the markers for the two oriented reads.
    const span<CompressedMarker> markers0 = markers[orientedReadId0.getValue()];
    const span<CompressedMarker> markers1 = markers[orientedReadId1.getValue()];
    const int64_t markerCount0 = markers0.size();
    const int64_t markerCount1 = markers1.size();
    const int64_t firstMarkerId0 = markers.begin(orientedReadId0.getValue()) - markers.begin();
    const int64_t firstMarkerId1 = markers.begin(orientedReadId1.getValue()) - markers.begin();



    // Assuming this offset, the induced alignment covers a portion of each
    // of the two oriented reads.
    const int64_t iOffset = uint64_t(offset);
    int64_t ordinal0Begin, ordinal0End;
    int64_t ordinal1Begin, ordinal1End;
    if(iOffset >= 0) {
        ordinal0Begin = iOffset;
        ordinal0End = min(markerCount0, markerCount1 + iOffset);
        ordinal1Begin = 0;
        ordinal1End = ordinal0End - iOffset;
    } else {
        ordinal1Begin = -iOffset;
        ordinal1End = min(markerCount1, markerCount0 - iOffset);
        ordinal0Begin = 0;
        ordinal0End = ordinal1End + iOffset;
    }
    SHASTA_ASSERT(ordinal0Begin >= 0);
    SHASTA_ASSERT(ordinal0Begin <  markerCount0);
    SHASTA_ASSERT(ordinal1Begin >= 0);
    SHASTA_ASSERT(ordinal1Begin <  markerCount1);
    SHASTA_ASSERT(ordinal0End   >  0);
    SHASTA_ASSERT(ordinal0End   <=  markerCount0);
    SHASTA_ASSERT(ordinal1End   >  0);
    SHASTA_ASSERT(ordinal1End   <=  markerCount1);
    const int64_t m =  ordinal0End - ordinal0Begin;
    SHASTA_ASSERT(m == ordinal1End - ordinal1Begin);


    // Use work0 and work1 to keep tract of markers in these two ranges.
    work0.clear();
    work1.clear();
    work0.resize(m, false);
    work1.resize(m, false);

    // Set to true the ones that are in the induced alignment.
    for(const auto& d: inducedAlignment.data) {
        if(d.ordinal0 >= ordinal0Begin and d.ordinal0 < ordinal0End) {
            work0[d.ordinal0 - ordinal0Begin] = true;
        }
        if(d.ordinal1 >= ordinal1Begin and d.ordinal1 < ordinal1End) {
            work1[d.ordinal1 - ordinal1Begin] = true;
        }
    }

    // Set to true the ones that are not associated with a marker graph vertex.
    for(int64_t ordinal0=ordinal0Begin; ordinal0!=ordinal0End; ++ordinal0) {
        if(markerGraph.vertexTable[firstMarkerId0 + ordinal0] == MarkerGraph::invalidCompressedVertexId) {
            work0[ordinal0 - ordinal0Begin] = true;
        }
    }
    for(int64_t ordinal1=ordinal1Begin; ordinal1!=ordinal1End; ++ordinal1) {
        if(markerGraph.vertexTable[firstMarkerId1 + ordinal1] == MarkerGraph::invalidCompressedVertexId) {
            work1[ordinal1 - ordinal1Begin] = true;
        }
    }



    // Any entries still set to false in the work arrays correspond to markers
    // that were not in the induced alignment, yet have an associated
    // marker graph vertex.
    // If we find a long streak of these entries, this is a bad induced alignment and
    // we return false;
    uint64_t lastTrueIndex0 = 0;
    for(uint64_t index0=0; index0<work0.size(); index0++) {
        if(work0[index0]) {
            lastTrueIndex0 = index0;
        } else {
            if(index0 - lastTrueIndex0 > criteria.maxSkip) {
                return false;
            }
        }
    }
    if(work0.size() - lastTrueIndex0 > criteria.maxSkip) {
        return false;
    }
    uint64_t lastTrueIndex1 = 0;
    for(uint64_t index1=0; index1<work1.size(); index1++) {
        if(work1[index1]) {
            lastTrueIndex1 = index1;
        } else {
            if(index1 - lastTrueIndex1> criteria.maxSkip) {
                return false;
            }
        }
    }
    if(work1.size() - lastTrueIndex1 > criteria.maxSkip) {
        return false;
    }



    // If getting here, this is a good induced alignment.
    return true;
}
