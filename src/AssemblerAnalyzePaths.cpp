// Shasta.
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
// Shasta.
#include "seqan.hpp"
using namespace shasta;

// Standard library.
#include "array.hpp"
#include "fstream.hpp"
#include <map>
#include <set>



void Assembler::computePseudoPath(
    OrientedReadId orientedReadId,

    // The marker graph path computed using computeOrientedReadMarkerGraphPath.
    // This is computed by this function - it does not need to be filled in
    // in advance.
    vector<MarkerGraph::EdgeId>& path,
    vector< pair<uint32_t, uint32_t> >& pathOrdinals,

    // The pseudo-path computed by this function.
    PseudoPath& pseudoPath) const
{
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    using SegmentId = AssemblyGraph::EdgeId;

    // Compute the marker graph path.
    const uint64_t markerCount = markers.size(orientedReadId.getValue());
    if(markerCount < 2) {
        pathOrdinals.clear();
        path.clear();
    } else {
        computeOrientedReadMarkerGraphPath(
            orientedReadId,
            0, uint32_t(markerCount - 1),
            path, pathOrdinals);
        SHASTA_ASSERT(path.size() == pathOrdinals.size());
    }



    // Now compute the pseudo-path.
    pseudoPath.clear();
    pseudoPath.clear();
    PseudoPathEntry pseudoPathEntry;
    pseudoPathEntry.segmentId = std::numeric_limits<SegmentId>::max();
    for(uint64_t i=0; i<path.size(); i++) {
        const MarkerGraph::EdgeId markerGraphEdgeId = path[i];
        const pair<uint32_t, uint32_t>& ordinals = pathOrdinals[i];

        // Get the corresponding assembly graph segments.
        const span<const pair<SegmentId, uint32_t> > v =
            assemblyGraph.markerToAssemblyTable[markerGraphEdgeId];

        // If no segments, skip.
        if(v.size() == 0) {
            continue;
        }

        // If detangling was used, there can be more than one,
        // and we don't want this here.
        SHASTA_ASSERT(v.size() == 1);

        // There is only one segment.
        const SegmentId segmentId = v.front().first;
        const uint32_t positionInSegment = v.front().second;

        // Same as the previous.
        if(segmentId == pseudoPathEntry.segmentId) {
            pseudoPathEntry.lastOrdinal = ordinals.second;
            pseudoPathEntry.lastPosition = positionInSegment;
            ++pseudoPathEntry.markerGraphEdgeCount;
            continue;
        }

        // This is the next segment edge encountered
        // by this oriented read along its marker graph path.
        if(pseudoPathEntry.segmentId != std::numeric_limits<SegmentId>::max()) {
            pseudoPath.push_back(pseudoPathEntry);
        }
        pseudoPathEntry.segmentId = segmentId;
        pseudoPathEntry.firstOrdinal = ordinals.first;
        pseudoPathEntry.lastOrdinal = ordinals.second;
        pseudoPathEntry.firstPosition = positionInSegment;
        pseudoPathEntry.lastPosition = positionInSegment;
        pseudoPathEntry.markerGraphEdgeCount = 1;
    }

    // Add the last entry.
    if(pseudoPathEntry.segmentId != std::numeric_limits<SegmentId>::max()) {
        pseudoPath.push_back(pseudoPathEntry);
    }
}



// Write the pseudo-path of an oriented read to a csv file.
// The pseudo-path is the sequence of assembly graph edges
// (not necsssarily all adjacent, so not necessatily a path)
// encountered by the oriented read.
void Assembler::writePseudoPath(ReadId readId, Strand strand) const
{
    // Compute the pseudo path.
    const OrientedReadId orientedReadId(readId, strand);
    vector<MarkerGraph::EdgeId> path;
    vector< pair<uint32_t, uint32_t> > pathOrdinals;
    PseudoPath pseudoPath;
    computePseudoPath(orientedReadId, path, pathOrdinals, pseudoPath);

    // Write it out.
    ofstream csv("PseudoPath.csv");
    csv << "Segment id,First ordinal,Last ordinal,"
        "First position in segment,Last position in segment, Marker graph edge count\n";
    for(const PseudoPathEntry& pseudoPathEntry: pseudoPath) {
        csv << pseudoPathEntry.segmentId << ",";
        csv << pseudoPathEntry.firstOrdinal << ",";
        csv << pseudoPathEntry.lastOrdinal << ",";
        csv << pseudoPathEntry.firstPosition << ",";
        csv << pseudoPathEntry.lastPosition << ",";
        csv << pseudoPathEntry.markerGraphEdgeCount << "\n";
    }
}



// Get the vector of segments corresponding to a PseudoPath.
void Assembler::getPseudoPathSegments(
    const PseudoPath& pseudoPath,
    vector<AssemblyGraph::EdgeId>& segmentIds)
{
    segmentIds.clear();
    for(const PseudoPathEntry& pseudoPathEntry: pseudoPath) {
        segmentIds.push_back(pseudoPathEntry.segmentId);
    }
}



void Assembler::alignPseudoPaths(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1)
{
    using SegmentId = AssemblyGraph::EdgeId;
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;

    // Parameters that control the process below. EXPOSE WHEN CODE STABILIZES. *********
    const int matchScore = 1;
    const int mismatchScore = -1;
    const int gapScore = -1;

    // Gather the oriented read ids.
    const array<OrientedReadId, 2> orientedReadIds =
        {OrientedReadId(readId0, strand0), OrientedReadId(readId1, strand1)};
    cout << "Aligning pseudo-paths of " << orientedReadIds[0] <<
        " and " << orientedReadIds[1] << endl;


    // Compute the two pseudo-paths.
    vector<MarkerGraph::EdgeId> path;
    vector< pair<uint32_t, uint32_t> > pathOrdinals;
    PseudoPath pseudoPath;
    array<vector<SegmentId>, 2> pseudoPathSegments;
    for(uint64_t i=0; i<2; i++) {
            computePseudoPath(orientedReadIds[i], path, pathOrdinals,
                pseudoPath);
            getPseudoPathSegments(pseudoPath, pseudoPathSegments[i]);
        cout << "The pseudo-path of " << orientedReadIds[i] <<
            " has " << pseudoPathSegments[i].size() << " segments." << endl;
    }

    // Align them.
    vector< pair<bool, bool> > alignment;
    const uint64_t alignmentScore = shasta::seqanAlign(
        pseudoPathSegments[0].begin(), pseudoPathSegments[0].end(),
        pseudoPathSegments[1].begin(), pseudoPathSegments[1].end(),
        matchScore,
        mismatchScore,
        gapScore,
        true, true,
        alignment);
    cout << "Alignment score " << alignmentScore << endl;
    cout << "Alignment length " << alignment.size() << endl;



    // Write out the alignment.
    uint64_t position0 = 0;
    uint64_t position1 = 0;
    uint64_t weakMatchCount =0;
    uint64_t strongMatchCount =0;
    uint64_t mismatchCount =0;
    uint64_t gapCount =0;
    uint64_t leftUnalignedCount =0;
    uint64_t rightUnalignedCount =0;
    ofstream csv("PseudoPathsAlignment.csv");
    for(const auto& p: alignment) {
        if(p.first) {
            const SegmentId segment0 = pseudoPathSegments[0][position0];
            csv << segment0;
            }
        csv << ",";
        if(p.second) {
            const SegmentId segment1 = pseudoPathSegments[1][position1];
            csv << segment1;
            }
        csv << ",";

        // Write an annotation column.
        if(p.first and p.second) {
            if(pseudoPathSegments[0][position0] != pseudoPathSegments[1][position1]) {
                csv << "Mismatch";
                ++mismatchCount;
            } else {
                // Match.
                // Decide if it is a strong or weak match.
                const SegmentId segmentId = pseudoPathSegments[0][position0];
                const AssemblyGraph::Edge& edge = assemblyGraph.edges[segmentId];
                const AssemblyGraph::VertexId v0 = edge.source;
                const AssemblyGraph::VertexId v1 = edge.target;
                const auto out0 = assemblyGraph.outDegree(v0);
                const auto in1 = assemblyGraph.inDegree(v1);
                if(out0==1 and in1==1) {
                    csv << "Weak match";
                    ++weakMatchCount;
                } else {
                    csv << "Strong match";
                    ++strongMatchCount;
                }
            }
        } else if(position0 == 0 or position1==0) {
            csv << "Left unaligned portion";
            ++leftUnalignedCount;
        } else if(
            position0 == pseudoPathSegments[0].size() or
            position1 == pseudoPathSegments[1].size()) {
            csv << "Right unaligned portion";
            ++rightUnalignedCount;
        } else if(not (p.first and p.second)) {
            csv << "Gap";
            ++gapCount;
        }
        csv << "\n";

        if(p.first) {
            ++position0;
        }
        if(p.second) {
            ++position1;
        }
    }
    SHASTA_ASSERT(position0 == pseudoPathSegments[0].size());
    SHASTA_ASSERT(position1 == pseudoPathSegments[1].size());

    const uint64_t matchCount = weakMatchCount + strongMatchCount;
    cout << "Total match "<< matchCount << endl;
    cout << "Strong match "<< strongMatchCount << endl;
    cout << "Weak match "<< weakMatchCount << endl;
    cout << "Mismatch "<< mismatchCount << endl;
    cout << "Gap "<< gapCount << endl;
    cout << "Left unaligned "<< leftUnalignedCount << endl;
    cout << "Right unaligned "<< rightUnalignedCount << endl;
    cout << "Mismatch/match ratio " << double(mismatchCount)/double(matchCount) << endl;
}
