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



// Evaluate an induced alignment.
// This is more sophisticated than InducedAlignment::evaluate,
// as it takes into account
// markers that don't correspond to a marker graph vertex.
bool Assembler::evaluateInducedAlignment(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const InducedAlignment& inducedAlignment,
    const InducedAlignmentCriteria& criteria,
    vector<uint64_t>& work)
{

    // Set the flag to control debug output.
#if 1
    const bool debug = false;
#else
    const bool debug =
        (orientedReadId0 == OrientedReadId(691, 0)  and
         orientedReadId1 == OrientedReadId(307, 0))
        or
        (orientedReadId0 == OrientedReadId(307, 0) and
         orientedReadId1 == OrientedReadId(691, 0));
#endif
   if(debug) {
       cout << "Evaluating induced alignment " << orientedReadId0 << " " << orientedReadId1 << endl;
   }

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
    const int64_t iOffset = int64_t(offset);
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

    if(debug) {
        cout << "iOffset "<< iOffset << endl;
        cout << "ordinal0Begin "<< ordinal0Begin << endl;
        cout << "ordinal0End "<< ordinal0End << endl;
        cout << "ordinal1Begin "<< ordinal1Begin << endl;
        cout << "ordinal1End "<< ordinal1End << endl;
        cout << "m "<< m << endl;
    }



    // For each of the two oriented reads, we now have
    // an expected alignment range based on the computed offset.
    // We locate large gaps of the alignment in this range.
    // If the large gaps contain markers associated to marker graph vertices,
    // this is a bad induced alignment.

    // Check large gaps on orientedReadId0.
    work.clear();
    for(const auto& d: inducedAlignment.data) {
        if(d.ordinal0 >= ordinal0Begin and d.ordinal0 < ordinal0End) {
            work.push_back(d.ordinal0);
        }
    }
    sort(work.begin(), work.end());
    for(uint64_t i=0; i<=work.size(); i++) {
        const uint64_t begin = (i==0) ? ordinal0Begin : work[i-1];
        const uint64_t end = (i==work.size()) ? ordinal0End : work[i]+1;
        if(end - begin > criteria.maxSkip) {
            // This is a long gap.
            for(uint64_t ordinal0=begin+1; ordinal0!=end; ordinal0++) {
                if(markerGraph.vertexTable[firstMarkerId0 + ordinal0] != MarkerGraph::invalidCompressedVertexId) {
                    if(debug) {
                        cout << "Bad induced alignment for " << orientedReadId0 <<
                            " at ordinal " << ordinal0 <<
                            ", gap " << begin << " " << end << endl;
                    }
                    return false;
                }
            }
        }
    }

    // Check large gaps on orientedReadId1.
    work.clear();
    for(const auto& d: inducedAlignment.data) {
        if(d.ordinal1 >= ordinal1Begin and d.ordinal1 < ordinal1End) {
            work.push_back(d.ordinal1);
        }
    }
    sort(work.begin(), work.end());
    for(uint64_t i=0; i<=work.size(); i++) {
        const uint64_t begin = (i==0) ? ordinal1Begin : work[i-1];
        const uint64_t end = (i==work.size()) ? ordinal1End : work[i]+1;
        if(end - begin > criteria.maxSkip) {
            // This is a long gap.
            for(uint64_t ordinal1=begin+1; ordinal1!=end; ordinal1++) {
                if(markerGraph.vertexTable[firstMarkerId1 + ordinal1] != MarkerGraph::invalidCompressedVertexId) {
                    if(debug) {
                        cout << "Bad induced alignment for " << orientedReadId1 <<
                            " at ordinal " << ordinal1 <<
                            ", gap " << begin << " " << end << endl;
                    }
                    return false;
                }
            }
        }
    }


#if 0
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
    if(debug) {
        for(int64_t i=0; i<m; i++) {
            cout << "Intermediate " <<
                i << " " <<
                i+ordinal0Begin << " " <<
                i+ordinal1Begin << " " <<
                int(work0[i]) << " " <<
                int(work1[i]) << endl;
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


    if(debug) {
        for(int64_t i=0; i<m; i++) {
            cout << "Final " <<
                i << " " <<
                i+ordinal0Begin << " " <<
                i+ordinal1Begin << " " <<
                int(work0[i]) << " " <<
                int(work1[i]) << endl;
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
#endif


    // If getting here, this is a good induced alignment.
    return true;
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
