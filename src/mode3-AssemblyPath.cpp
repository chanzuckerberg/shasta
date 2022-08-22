// Shasta.
#include "mode3-AssemblyPath.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "ConsensusCaller.hpp"
#include "deduplicate.hpp"
#include "html.hpp"
#include "invalid.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "Reads.hpp"
#include "mode3.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

// Spoa.
#include "spoa/spoa.hpp"

// Seqan.
#include <seqan/align.h>

// Standard library.
#include "fstream.hpp"



// Assemble sequence for an AssemblyPath.
void AssemblyPath::assemble(const AssemblyGraph& assemblyGraph)
{
    // Assemble each segment on the path.
    assembleSegments(assemblyGraph);

    // Assemble links in this assembly path.
    assembleLinks(assemblyGraph);

    writeSegmentSequences();
    writeLinkSequences(assemblyGraph);

    assemble();
}



// Find the oriented reads to be used to assemble
// links between the segment at position0
// in the assembly path (which must be a reference segment)
// and the next reference segment in the path.
void AssemblyPath::findOrientedReadsForLinks(
    uint64_t position0,
    const AssemblyGraph& assemblyGraph,
    vector<OrientedReadId>& orientedReadsForLinks) const
{

    // Sanity checks.
    SHASTA_ASSERT(position0 < segments.size());
    const auto p0 = segments[position0];
    SHASTA_ASSERT(p0.second);   // Must be a reference segment.
    const uint64_t segmentId0 = p0.first;

    // Add the oriented reads in segmentId0, the reference segment
    // at position0.
    orientedReadsForLinks.clear();
    for(const auto& p: assemblyGraph.assemblyGraphJourneyInfos[segmentId0]) {
        orientedReadsForLinks.push_back(p.first);
    }

    // Look for the next reference segment in the path.
    uint64_t segmentId1 = invalid<uint64_t>;
    for(uint64_t position1=position0+1; position1<segments.size(); position1++) {
        const auto& q = segments[position1];
        if(q.second) {
            segmentId1 = q.first;
            break;
        }
    }
    SHASTA_ASSERT(segmentId1 != invalid<uint64_t>);

    // Add the oriented reads in segmentId1.
    for(const auto& p: assemblyGraph.assemblyGraphJourneyInfos[segmentId1]) {
        orientedReadsForLinks.push_back(p.first);
    }

    // Deduplicate and sort.
    deduplicate(orientedReadsForLinks);
}



// Assemble links in this assembly path.
void AssemblyPath::assembleLinks(const AssemblyGraph& assemblyGraph)
{
    const bool debug = false;
    SHASTA_ASSERT((assemblyGraph.k % 2) == 0);

    // Resize the vectors that describe link sequences.
    linksRleSequence.resize(segments.size()-1);
    linksRepeatCounts.resize(segments.size()-1);
    skipAtSegmentBegin.resize(segments.size());
    skipAtSegmentEnd.resize(segments.size());

    // Don't skip any bases at the beginning of the first
    // segment and at the end of the last segment.
    skipAtSegmentBegin.front() = 0;
    skipAtSegmentEnd.back() = 0;

    // The list of oriented reads that can be used to assemble links.
    // These are the oriented reads that appear in either the previous or next
    // reference segment in the path.
    // This list changes every time we encounter a new reference segment.
    // It is stored sorted.
    vector<OrientedReadId> orientedReadsForLinks;

    ofstream html("Msa.html");

    // Loop over pairs of consecutive segments in the path.
    for(uint64_t position0=0; position0<segments.size()-1 ; position0++) {
        const uint64_t position1 = position0 + 1;

        // Access the consecutive segments in this pair.
        // We will process the link between segmentId0 and segmentId1.
        const auto& p0 = segments[position0];
        const uint64_t segmentId0 = p0.first;
        const bool isReferenceSegment0 = p0.second;
        const AssembledSegment& assembledSegment0 = assembledSegments[position0];
        const auto& p1 = segments[position1];
        const uint64_t segmentId1 = p1.first;
        const AssembledSegment& assembledSegment1 = assembledSegments[position1];

        // The first segment in the path must be a reference segment.
        if(position0 == 0) {
            SHASTA_ASSERT(isReferenceSegment0);
        }

        // If this is a reference segment, update
        // the list of oriented reads that can be used to assemble links.
        // These are the oriented reads that appear in either the previous or next
        // reference segment in the path.
        if(isReferenceSegment0) {
            findOrientedReadsForLinks(position0, assemblyGraph, orientedReadsForLinks);
        }

        // If segmentId0 and segmentId1 are consecutive in the marker graph,
        // just skip the last k/2 RLE bases of segmentId0
        // and the first k/2 RLE bases of seegmentId1.
        // and take them from segmentId1 (they are identical because
        // they correspond to the same marker graph vertex).
        const auto markerGraphPath0 = assemblyGraph.markerGraphPaths[segmentId0];
        const auto markerGraphPath1 = assemblyGraph.markerGraphPaths[segmentId1];
        const MarkerGraph::Edge lastEdge0 = assemblyGraph.markerGraph.edges[markerGraphPath0.back()];
        const MarkerGraph::Edge firstEdge1 = assemblyGraph.markerGraph.edges[markerGraphPath1.front()];
        if(lastEdge0.target == firstEdge1.source) {
            skipAtSegmentEnd[position0] = assemblyGraph.k/2;
            skipAtSegmentBegin[position1] = assemblyGraph.k/2;

            // Leave empty the sequence for the link between these segments.
            SHASTA_ASSERT(linksRleSequence[position0].empty());
            SHASTA_ASSERT(linksRepeatCounts[position0].empty());

            // We are done.
            continue;
        }

        // If getting here, this segment and the next are not adjacent in the marker graph.
        // We need to assemble the link between them.
        // We can only use oriented reads that are in the link
        // and also in orientedReadsForLinks.

        // Locate the link between segmentId0 and segmentId1.
        uint64_t linkId01 = invalid<uint64_t>;
        for(uint64_t linkId: assemblyGraph.linksBySource[segmentId0]) {
            if(assemblyGraph.links[linkId].segmentId1 == segmentId1) {
                linkId01 = linkId;
                break;
            }
        }
        SHASTA_ASSERT(linkId01!= invalid<uint64_t>);
        if(debug) {
            cout << "Assembling link " << linkId01 << " " << segmentId0 << "->" << segmentId1 << endl;
        }

        // First, find:
        // - The position in segmentId0 of the leftmost transition.
        // - The position in segmentId1 of the rightmost transition.
        uint64_t minEdgePosition0 = markerGraphPath0.size();
        uint64_t maxEdgePosition1 = 0;
        for(const auto& p: assemblyGraph.transitions[linkId01]) {
            const OrientedReadId orientedReadId = p.first;

            // If not in orientedReadsForLinks, skip it.
            if(not std::binary_search(
                orientedReadsForLinks.begin(), orientedReadsForLinks.end(),
                orientedReadId)) {
                continue;
            }

            // Access the transition from segmentId0 to segmentId1 for this oriented read.
            const AssemblyGraph::Transition& transition = p.second;

            minEdgePosition0 = min(minEdgePosition0, uint64_t(transition[0].position));
            maxEdgePosition1 = max(maxEdgePosition1, uint64_t(transition[1].position));
        }

        // When getting here:
        // - minEdgePosition0 is the leftmost position of the transitions in path0.
        // - maxEdgePosition1 is the rightmost position of the transitions in path1.
        // These positions are edge positions in markerGraphPath0 and markerGraphPath1.
        // We will do a multiple sequence alignment of the oriented reads,
        // using the sequence of segmentId0 to extend to the left all reads to minEdgePosition0,
        // and using the sequence of segmentId1 to extend to the right all reads to maxEdgePosition1,

        // Get the corresponding vertex positions in segmentId0 and segmentId1.
        const uint64_t minVertexPosition0 = minEdgePosition0 + 1;
        const uint64_t maxVertexPosition1 = maxEdgePosition1;



        // Now extract the portion of each oriented read sequence that
        // will be used to assemble this link.
        vector<OrientedReadId> orientedReadIdsForAssembly;
        vector< vector<Base> > orientedReadsSequencesForAssembly;
        vector< vector<uint64_t> > orientedReadsRepeatCountsForAssembly;
        for(const auto& p: assemblyGraph.transitions[linkId01]) {
            const OrientedReadId orientedReadId = p.first;

            // If not in orientedReadsForLinks, skip it.
            if(not std::binary_search(
                orientedReadsForLinks.begin(), orientedReadsForLinks.end(),
                orientedReadId)) {
                continue;
            }

            // Access the transition from segmentId0 to segmentId1 for this oriented read.
            const AssemblyGraph::Transition& transition = p.second;

            // Get the ordinals of the last appearance of this oriented
            // read on segmentId0 and the first on segmentId1,
            // and the corresponding markers.
            const uint32_t ordinal0 = transition[0].ordinals[1];
            const uint32_t ordinal1 = transition[1].ordinals[0];
            const CompressedMarker& marker0 = assemblyGraph.markers[orientedReadId.getValue()][ordinal0];
            const CompressedMarker& marker1 = assemblyGraph.markers[orientedReadId.getValue()][ordinal1];

            // Get the positions of these markers on the oriented read.
            // If using RLE, these are RLE positions.
            const uint32_t position0 = marker0.position;
            const uint32_t position1 = marker1.position;

            // Extract the sequence between these markers (including the markers).
            vector<Base> orientedReadSequence;
            vector<uint8_t> orientedReadRepeatCounts;
            if(assemblyGraph.readRepresentation == 1) {
                // RLE.
                for(uint64_t position=position0; position<position1+assemblyGraph.k; position++) {
                    Base b;
                    uint8_t r;
                    tie(b, r) = assemblyGraph.reads.getOrientedReadBaseAndRepeatCount(orientedReadId, uint32_t(position));
                    orientedReadSequence.push_back(b);
                    orientedReadRepeatCounts.push_back(r);
                }
            } else {
                // Raw sequence.
                for(uint64_t position=position0; position<position1+assemblyGraph.k; position++) {
                    const Base b = assemblyGraph.reads.getOrientedReadBase(orientedReadId, uint32_t(position));
                    orientedReadSequence.push_back(b);
                    orientedReadRepeatCounts.push_back(uint8_t(1));
                }
            }

            // We need to extend the sequence of this read to the left,
            // using segmentId0 sequence, up to minVertexPosition0,
            // so the portions of all reads we will be using for the MSA
            // all begin in the same place.
            vector<Base> leftSequence;
            vector<uint32_t> leftRepeatCounts;
            const uint64_t vertexPosition0 = transition[0].position + 1;  // Add 1 to get vertex position.
            const uint64_t begin0 = assembledSegment0.vertexOffsets[minVertexPosition0];
            const uint64_t end0 = assembledSegment0.vertexOffsets[vertexPosition0];
            for(uint64_t position=begin0; position!=end0; position++) {
                leftSequence.push_back(assembledSegment0.runLengthSequence[position]);
                leftRepeatCounts.push_back(assembledSegment0.repeatCounts[position]);
            }

            vector<Base> rightSequence;
            vector<uint32_t> rightRepeatCounts;
            const uint64_t vertexPosition1 = transition[1].position;
            const uint64_t begin1 = assembledSegment1.vertexOffsets[vertexPosition1] + assemblyGraph.k;
            const uint64_t end1 = assembledSegment1.vertexOffsets[maxVertexPosition1] + assemblyGraph.k;
            for(uint64_t position=begin1; position!=end1; position++) {
                rightSequence.push_back(assembledSegment1.runLengthSequence[position]);
                rightRepeatCounts.push_back(assembledSegment1.repeatCounts[position]);
            }

            // Construct the extended sequence for this oriented read,
            // to be used in the MSA.
            vector<Base> orientedReadExtendedSequence;
            vector<uint64_t> orientedReadExtendedRepeatCounts;
            copy(leftSequence.begin(), leftSequence.end(), back_inserter(orientedReadExtendedSequence));
            copy(leftRepeatCounts.begin(), leftRepeatCounts.end(), back_inserter(orientedReadExtendedRepeatCounts));
            copy(orientedReadSequence.begin(), orientedReadSequence.end(), back_inserter(orientedReadExtendedSequence));
            copy(orientedReadRepeatCounts.begin(), orientedReadRepeatCounts.end(), back_inserter(orientedReadExtendedRepeatCounts));
            copy(rightSequence.begin(), rightSequence.end(), back_inserter(orientedReadExtendedSequence));
            copy(rightRepeatCounts.begin(), rightRepeatCounts.end(), back_inserter(orientedReadExtendedRepeatCounts));

            orientedReadIdsForAssembly.push_back(orientedReadId);
            orientedReadsSequencesForAssembly.push_back(orientedReadExtendedSequence);
            orientedReadsRepeatCountsForAssembly.push_back(orientedReadExtendedRepeatCounts);

            if(debug) {
                copy(orientedReadExtendedSequence.begin(), orientedReadExtendedSequence.end(), ostream_iterator<Base>(cout));
                cout << " " << orientedReadId << endl;
            }
        }

        // Compute the consensus sequence for the link.
        vector<Base> consensusRleSequence;
        vector<uint64_t> consensusRepeatCounts;
        html<< "<h2>Link " << linkId01 << "</h2>\n";
        computeLinkConsensusUsingSpoa(
            orientedReadIdsForAssembly,
            orientedReadsSequencesForAssembly,
            orientedReadsRepeatCountsForAssembly,
            assemblyGraph.readRepresentation,
            assemblyGraph.consensusCaller,
            debug,
            html,
            consensusRleSequence,
            consensusRepeatCounts
            );
        SHASTA_ASSERT(consensusRleSequence.size() == consensusRepeatCounts.size());

        if(debug) {
            cout << "Consensus RLE sequence length before trimming " << consensusRleSequence.size() << endl;
            cout << "Portion of segment on left involved in the MSA begins at position " <<
                assembledSegment0.vertexOffsets[minVertexPosition0] << endl;
            cout << "Portion of segment on right involved in the MSA ends at position " <<
                assembledSegment1.vertexOffsets[maxVertexPosition1] + assemblyGraph.k << endl;
        }

        // Count the number of identical (RLE) bases at the beginning of the
        // link consensus sequence and of the segmentId0 sequence portion
        // involved in assembling this link.
        uint64_t identicalOnLeft = 0;
        const uint64_t begin0 = assembledSegment0.vertexOffsets[minVertexPosition0];
        const uint64_t end0 = assembledSegment0.runLengthSequence.size();
        for(uint64_t i=begin0; (i!=end0 and (i-begin0)<consensusRleSequence.size()); i++) {
            if(consensusRleSequence[i-begin0] == assembledSegment0.runLengthSequence[i]) {
                // cout << "*** " << begin0 << " " << end0 << " " << i << endl;
                ++identicalOnLeft;
            } else {
                break;
            }
        }
        if(debug) {
            cout << "Identical on left: " << identicalOnLeft << endl;
        }

        // Count the number of identical (RLE) bases at the end of the
        // link consensus sequence and the beginning of segmentId1 .
        uint64_t identicalOnRight = 0;
        const uint64_t end1 = assembledSegment1.vertexOffsets[maxVertexPosition1] + assemblyGraph.k;
        for(uint64_t i=end1-1; ; i--) {
            const uint64_t j = consensusRleSequence.size() - (end1 - i);
            if(consensusRleSequence[j] == assembledSegment1.runLengthSequence[i]) {
                // cout << "*** " << i << " " << assembledSegment1.runLengthSequence[i] << " " <<
                //     j << " " << consensusRleSequence[j] << endl;
                ++identicalOnRight;
            } else {
                break;
            }
            if(i == 0) {
                break;
            }
            if(j == 0) {
                break;
            }
        }
        identicalOnRight = min(identicalOnRight, consensusRleSequence.size()-identicalOnLeft);
        if(debug) {
            cout << "Identical on right: " << identicalOnRight << endl;
        }

        // Trim these identical bases from the link consensus sequence.
        const uint64_t trimmedConsensusLength = consensusRleSequence.size() - identicalOnLeft - identicalOnRight;
        vector<Base> trimmedConsensusRleSequence(trimmedConsensusLength);
        vector<uint64_t> trimmedlconsensusRepeatCounts(trimmedConsensusLength);
        copy(
            consensusRleSequence.begin() + identicalOnLeft,
            consensusRleSequence.begin() + identicalOnLeft + trimmedConsensusLength,
            trimmedConsensusRleSequence.begin());
        copy(
            consensusRepeatCounts.begin() + identicalOnLeft,
            consensusRepeatCounts.begin() + identicalOnLeft + trimmedConsensusLength,
            trimmedlconsensusRepeatCounts.begin());

        // Store this as the link sequence.
        linksRleSequence[position0] = trimmedConsensusRleSequence;
        linksRepeatCounts[position0] = trimmedlconsensusRepeatCounts;

        // Compute and store the number of bases to be skipped at the end of segmentId0
        // and at the beginning of segmentId1.
        skipAtSegmentEnd[position0] =
            assembledSegment0.runLengthSequence.size() -
            assembledSegment0.vertexOffsets[minVertexPosition0] -
            identicalOnLeft;
        skipAtSegmentBegin[position1] =
            assembledSegment1.vertexOffsets[maxVertexPosition1] + assemblyGraph.k
            - identicalOnRight;
    }
}



// Assemble links in this assembly path.
void AssemblyPath::assembleLinksOld(
    const AssemblyGraph& assemblyGraph,
    bool debug)
{
    ofstream html;
    ofstream linksFasta;
    if(debug) {
        linksFasta.open("PathLinksSequence.fasta");
        html.open("Msa.html");
        writeHtmlBegin(html, "Mode 3 path assembly");
    }

    // The list of oriented reads that can be used to assemble links.
    // These are the oriented reads that appear in either the previous or next
    // reference segment in the path.
    // This list changes every time we encounter a new reference segment.
    // It is kept sorted for efficiency.
    vector<OrientedReadId> orientedReadsForLinks;


    // Assemble segment and link sequence into path sequence.
    string sequence;
    for(uint64_t i=0; ; i++) {
        const auto& p = segments[i];
        const uint64_t segmentId0 = p.first;
        if(debug) {
            cout << "Working on segment " << segmentId0 << " at position " << i <<
                " in assembly path." << endl;
        }
        const bool isReferenceSegment = p.second;
        AssembledSegment& assembledSegment0 = assembledSegments[i];

        // The first segment in the path must be a reference segment.
        if(i == 0) {
            SHASTA_ASSERT(isReferenceSegment);
        }



        // Add the sequence of this segment.
        for(const Base b: assembledSegment0.rawSequence) {
            sequence.push_back(b.character());
        }

        // If this is the last segment, we are done.
        if(i == segments.size() - 1) {
            break;
        }



        // If this is a reference segment, update
        // the list of oriented reads that can be used to assemble links.
        // These are the oriented reads that appear in either the previous or next
        // reference segment in the path.
        if(isReferenceSegment) {
            orientedReadsForLinks.clear();

            // Add the oriented reads in segmentId0, our new
            // reference segment.
            for(const auto& p: assemblyGraph.assemblyGraphJourneyInfos[segmentId0]) {
                orientedReadsForLinks.push_back(p.first);
            }

            // Look for the next reference segment in the path.
            uint64_t nextReferenceSegmentId = invalid<uint64_t>;
            for(uint64_t j=i+1; j<segments.size(); j++) {
                const auto& q = segments[i];
                if(q.second) {
                    nextReferenceSegmentId = q.first;
                    break;
                }
            }
            SHASTA_ASSERT(nextReferenceSegmentId != invalid<uint64_t>);
            for(const auto& p: assemblyGraph.assemblyGraphJourneyInfos[nextReferenceSegmentId]) {
                orientedReadsForLinks.push_back(p.first);
            }

            // Deduplicate and sort.
            deduplicate(orientedReadsForLinks);
        }

        // This is not the last segment. Access the next segment in the path.
        const uint64_t segmentId1 = segments[i+1].first;
        AssembledSegment& assembledSegment1 = assembledSegments[i+1];



        // If the segments are consecutive in the marker graph,
        // just take out the last k RLE bases from the sequence.
        // We will get them from the next segment.
        const auto path0 = assemblyGraph.markerGraphPaths[segmentId0];
        const auto path1 = assemblyGraph.markerGraphPaths[segmentId1];
        const MarkerGraph::Edge lastEdge0 = assemblyGraph.markerGraph.edges[path0.back()];
        const MarkerGraph::Edge firstEdge1 = assemblyGraph.markerGraph.edges[path1.front()];
        if(lastEdge0.target == firstEdge1.source) {

            // Compute the number of bases to be removed.
            uint64_t removeCount = 0;
            for(uint64_t j=assembledSegment0.repeatCounts.size() - assemblyGraph.k;
                j<assembledSegment0.repeatCounts.size(); j++) {
                removeCount += assembledSegment0.repeatCounts[j];
            }

            // Remove them from the sequence we already have.
            sequence.resize(sequence.size() - removeCount);

            // We are done.
            continue;
        }



        // If getting here, this segment and the next are not adjacent in the marker graph.
        // We need to assemble the link between them.
        // We can only use oriented reads that are in the link
        // and also in orientedReadsForLinks.

        // Locate the link between segmentId0 and segmentId1.
        uint64_t linkId01 = invalid<uint64_t>;
        for(uint64_t linkId: assemblyGraph.linksBySource[segmentId0]) {
            if(assemblyGraph.links[linkId].segmentId1 == segmentId1) {
                linkId01 = linkId;
                break;
            }
        }
        SHASTA_ASSERT(linkId01!= invalid<uint64_t>);

        // Access the oriented reads that transition over this link.
        const auto transitions01 = assemblyGraph.transitions[linkId01];

        if(debug) {
            cout << "Assembling link " << linkId01 << " " << segmentId0 << "->" << segmentId1 <<
                " with coverage " << transitions01.size() <<
                " and " << orientedReadsForLinks.size() <<
                " reference oriented reads." << endl;
            cout << segmentId0 << " path length " << path0.size() << endl;
            cout << segmentId1 << " path length " << path1.size() << endl;
        }



        // Consider all oriented reads that are both in transitions01 and in
        // orientedReadsForLinks. These are the oriented reads that we can
        // actually use to assemble the link.



        // First, find:
        // - The position in segmentId0 of the leftmost transition.
        // - The position in segmentId1 of the rightmost transition.
        uint64_t minEdgePosition0 = path0.size();
        uint64_t maxEdgePosition1 = 0;
        for(const auto& p: assemblyGraph.transitions[linkId01]) {
            const OrientedReadId orientedReadId = p.first;

            // If not in orientedReadsForLinks, skip it.
            if(not std::binary_search(
                orientedReadsForLinks.begin(), orientedReadsForLinks.end(),
                orientedReadId)) {
                continue;
            }

            // Access the transition from segmentId0 to segmentId1 for this oriented read.
            const AssemblyGraph::Transition& transition = p.second;

            minEdgePosition0 = min(minEdgePosition0, uint64_t(transition[0].position));
            maxEdgePosition1 = max(maxEdgePosition1, uint64_t(transition[1].position));
        }
        if(debug) {
            cout << "minEdgePosition0 " << minEdgePosition0 << endl;
            cout << "maxEdgePosition1 " << maxEdgePosition1 << endl;
        }

        // When getting here:
        // - minPosition0 is the leftmost position of the transitions in path0.
        // - minPosition1 is the rightmost position of the transitions in path1.
        // These positions are edge positions in path0 and path1.
        // We will do a multiple sequence alignment of the oriented reads,
        // using the sequence of segmentId0 to extend to the left all reads to minPosition0,
        // and using the sequence of segmentId1 to extend to the right all reads to minPosition1,

        // Get the corresponding vertex positions in segmentId0 and segmentId1.
        const uint64_t minVertexPosition0 = minEdgePosition0 + 1;
        const uint64_t maxVertexPosition1 = maxEdgePosition1;
        if(debug) {
            cout << "minVertexPosition0 " << minVertexPosition0 << endl;
            cout << "maxVertexPosition1 " << maxVertexPosition1 << endl;
        }



        // Now extract the portion of each oriented read sequence that
        // will be used to assemble this link.
        uint64_t actualLinkCoverage = 0;
        vector<OrientedReadId> orientedReadIdsForAssembly;
        vector< vector<Base> > orientedReadsSequencesForAssembly;
        vector< vector<uint64_t> > orientedReadsRepeatCountsForAssembly;
        for(const auto& p: assemblyGraph.transitions[linkId01]) {
            const OrientedReadId orientedReadId = p.first;

            // If not in orientedReadsForLinks, skip it.
            if(not std::binary_search(
                orientedReadsForLinks.begin(), orientedReadsForLinks.end(),
                orientedReadId)) {
                continue;
            }

            // Access the transition from segmentId0 to segmentId1 for this oriented read.
            const AssemblyGraph::Transition& transition = p.second;
            ++actualLinkCoverage;

            // Get the ordinals of the last appearance of this oriented
            // read on segmentId0 and the first on segmentId1,
            // and the corresponding markers.
            const uint32_t ordinal0 = transition[0].ordinals[1];
            const uint32_t ordinal1 = transition[1].ordinals[0];
            if(debug) {
                cout << orientedReadId << " " << ordinal0 << " " << ordinal1 << endl;
            }
            const CompressedMarker& marker0 = assemblyGraph.markers[orientedReadId.getValue()][ordinal0];
            const CompressedMarker& marker1 = assemblyGraph.markers[orientedReadId.getValue()][ordinal1];

            // Get the positions of these markers on the oriented read.
            // If using RLE, these are RLE positions.
            const uint32_t position0 = marker0.position;
            const uint32_t position1 = marker1.position;

            // Extract the sequence between these markers (including the markers).
            vector<Base> orientedReadSequence;
            vector<uint8_t> orientedReadRepeatCounts;
            if(assemblyGraph.readRepresentation == 1) {
                // RLE.
                for(uint64_t position=position0; position<position1+assemblyGraph.k; position++) {
                    Base b;
                    uint8_t r;
                    tie(b, r) = assemblyGraph.reads.getOrientedReadBaseAndRepeatCount(orientedReadId, uint32_t(position));
                    orientedReadSequence.push_back(b);
                    orientedReadRepeatCounts.push_back(r);
                }
            } else {
                // Raw sequence.
                for(uint64_t position=position0; position<position1+assemblyGraph.k; position++) {
                    const Base b = assemblyGraph.reads.getOrientedReadBase(orientedReadId, uint32_t(position));
                    orientedReadSequence.push_back(b);
                    orientedReadRepeatCounts.push_back(uint8_t(1));
                }
            }

            if(debug) {
                copy(orientedReadSequence.begin(), orientedReadSequence.end(), ostream_iterator<Base>(cout));
                cout << endl;
                for(uint8_t r: orientedReadRepeatCounts) {
                    if(r < 10) {
                        cout << int(r);
                    } else {
                        cout << "*";
                    }
                }
                cout << endl;
            }

            // We need to extend the sequence of this read to the left,
            // using segmentId0 sequence, up to minVertexPosition0,
            // so the portions of all reads we will be using for the MSA
            // all begin in the same place.
            vector<Base> leftSequence;
            vector<uint32_t> leftRepeatCounts;
            const uint64_t vertexPosition0 = transition[0].position + 1;  // Add 1 to get vertex position.
            const uint64_t begin0 = assembledSegment0.vertexOffsets[minVertexPosition0];
            const uint64_t end0 = assembledSegment0.vertexOffsets[vertexPosition0];
            for(uint64_t position=begin0; position!=end0; position++) {
                leftSequence.push_back(assembledSegment0.runLengthSequence[position]);
                leftRepeatCounts.push_back(assembledSegment0.repeatCounts[position]);
            }
            if(debug) {
                cout << "Left extend sequence:" << endl;
                copy(leftSequence.begin(), leftSequence.end(), ostream_iterator<Base>(cout));
                cout << endl;
                for(uint32_t r: leftRepeatCounts) {
                    if(r < 10) {
                        cout << int(r);
                    } else {
                        cout << "*";
                    }
                }
                cout << endl;
            }

            vector<Base> rightSequence;
            vector<uint32_t> rightRepeatCounts;
            const uint64_t vertexPosition1 = transition[1].position;
            const uint64_t begin1 = assembledSegment1.vertexOffsets[vertexPosition1] + assemblyGraph.k;
            const uint64_t end1 = assembledSegment1.vertexOffsets[maxVertexPosition1] + assemblyGraph.k;
            // cout << "*** begin1 " << begin1 << endl;
            // cout << "*** end1 " << end1 << endl;
            for(uint64_t position=begin1; position!=end1; position++) {
                // cout << "*** " << position << " " << assembledSegment1.runLengthSequence[position] <<
                //     " " << int(assembledSegment1.repeatCounts[position]) << endl;
                rightSequence.push_back(assembledSegment1.runLengthSequence[position]);
                rightRepeatCounts.push_back(assembledSegment1.repeatCounts[position]);
            }
            if(debug) {
                cout << "Right extend sequence:" << endl;
                copy(rightSequence.begin(), rightSequence.end(), ostream_iterator<Base>(cout));
                cout << endl;
                for(uint32_t r: rightRepeatCounts) {
                    if(r < 10) {
                        cout << int(r);
                    } else {
                        cout << "*";
                    }
                }
                cout << endl;
            }

            // Construct the extended sequence for this oriented read,
            // to be used in the MSA.
            vector<Base> orientedReadExtendedSequence;
            vector<uint64_t> orientedReadExtendedRepeatCounts;
            copy(leftSequence.begin(), leftSequence.end(), back_inserter(orientedReadExtendedSequence));
            copy(leftRepeatCounts.begin(), leftRepeatCounts.end(), back_inserter(orientedReadExtendedRepeatCounts));
            copy(orientedReadSequence.begin(), orientedReadSequence.end(), back_inserter(orientedReadExtendedSequence));
            copy(orientedReadRepeatCounts.begin(), orientedReadRepeatCounts.end(), back_inserter(orientedReadExtendedRepeatCounts));
            copy(rightSequence.begin(), rightSequence.end(), back_inserter(orientedReadExtendedSequence));
            copy(rightRepeatCounts.begin(), rightRepeatCounts.end(), back_inserter(orientedReadExtendedRepeatCounts));

            orientedReadIdsForAssembly.push_back(orientedReadId);
            orientedReadsSequencesForAssembly.push_back(orientedReadExtendedSequence);
            orientedReadsRepeatCountsForAssembly.push_back(orientedReadExtendedRepeatCounts);

            if(debug) {
                cout << "Sequence for this oriented read to be used for MSA:" << endl;
                copy(orientedReadExtendedSequence.begin(), orientedReadExtendedSequence.end(), ostream_iterator<Base>(cout));
                cout << endl;
                for(uint64_t r: orientedReadExtendedRepeatCounts) {
                    if(r < 10) {
                        cout << int(r);
                    } else {
                        cout << "*";
                    }
                }
                cout << endl;
            }
        }
        if(debug) {
            cout << "Link assembly will use " << actualLinkCoverage <<
                " oriented reads." << endl;
        }

        // Compute the consensus sequence for the link.
        vector<Base> consensusRleSequence;
        vector<uint64_t> consensusRepeatCounts;
        if(true) {
            html<< "<h2>Link " << linkId01 << "</h2>\n";
        }
        computeLinkConsensusUsingSpoa(
            orientedReadIdsForAssembly,
            orientedReadsSequencesForAssembly,
            orientedReadsRepeatCountsForAssembly,
            assemblyGraph.readRepresentation,
            assemblyGraph.consensusCaller,
            debug,
            html,
            consensusRleSequence,
            consensusRepeatCounts
            );
        SHASTA_ASSERT(consensusRleSequence.size() == consensusRepeatCounts.size());
        if(debug) {
            linksFasta << ">L" << i << " " << linkId01 << "\n";
            for(uint64_t position=0; position<consensusRleSequence.size(); position++) {
                const Base base = consensusRleSequence[position];
                const uint64_t repeatCount = consensusRepeatCounts[position];
                for(uint64_t j=0; j<repeatCount; j++) {
                    linksFasta << base;
                }
            }
            linksFasta << "\n";
        }



        // Align the link consensus sequence with the sequence of the adjacent segment
        // on the left.
        // We don't need to use the entire segment sequence to compute
        // the alignment - only the final portion that was used to construct read
        // extended sequences for this link.
        using SeqanSequence = seqan::String<char>;
        using SeqanStringSet = seqan::StringSet<SeqanSequence>;
        using SeqanDepStringSet = seqan::StringSet<SeqanSequence, seqan::Dependent<> >;
        using SeqanAlignGraph = seqan::Graph<seqan::Alignment<SeqanDepStringSet> >;
        {
            // The sequence on the left is the previous segment.
            const uint64_t begin0 = assembledSegment0.vertexOffsets[minVertexPosition0];
            const uint64_t end0 = assembledSegment0.runLengthSequence.size();
            SeqanSequence sequence0;
            for(uint64_t i=begin0; i!=end0; i++) {
                appendValue(sequence0, assembledSegment0.runLengthSequence[i].character());
            }

            // The sequence on the right is the link.
            SeqanSequence sequence1;
            for(uint64_t i=0; i!=consensusRleSequence.size(); i++) {
                appendValue(sequence1, consensusRleSequence[i].character());
            }

            SeqanStringSet sequences;
            appendValue(sequences, sequence0);
            appendValue(sequences, sequence1);
            SeqanAlignGraph graph(sequences);
            const int64_t alignmentScore =  seqan::globalAlignment(
                graph,
                seqan::Score<int64_t, seqan::Simple>(1, -1, -1),
                seqan::AlignConfig<false, false, true, true>(),
                seqan::LinearGaps());
            if(debug) {
                cout << "Alignment of link consensus to left segment has score " << alignmentScore << endl;
                cout << graph << endl;
            }



            // Align the link consensus sequence with the sequence of the adjacent segment
            // on the right.
            // We don't need to use the entire segment sequence to compute
            // the alignment - only the initial portion that was used to construct read
            // extended sequences for this link.
            {
                // The sequence on the left is the link.
                SeqanSequence sequence0;
                for(uint64_t i=0; i!=consensusRleSequence.size(); i++) {
                    appendValue(sequence0, consensusRleSequence[i].character());
                }

                // The sequence on the right is the next segment.
                SeqanSequence sequence1;
                const uint64_t begin1 = 0;
                const uint64_t end1 = assembledSegment1.vertexOffsets[maxVertexPosition1] + assemblyGraph.k;
                for(uint64_t i=begin1; i!=end1; i++) {
                    appendValue(sequence1, assembledSegment1.runLengthSequence[i].character());
                }

                SeqanStringSet sequences;
                appendValue(sequences, sequence0);
                appendValue(sequences, sequence1);
                SeqanAlignGraph graph(sequences);
                const int64_t alignmentScore =  seqan::globalAlignment(
                    graph,
                    seqan::Score<int64_t, seqan::Simple>(1, -1, -1),
                    seqan::AlignConfig<true, true, false, false>(),
                    seqan::LinearGaps());
                if(debug) {
                    cout << "Alignment of link consensus to right segment has score " << alignmentScore << endl;
                    cout << graph << endl;
                }
            }
        }
    }



    // Write out the path sequence.
    if(debug) {
        ofstream fasta("PathSequence.fasta");
        fasta << ">path\n" << sequence << "\n";
    }

    if(debug) {
        cout << timestamp << "AssemblyGraph::assemblePathSequence ends." << endl;
        writeHtmlEnd(html);
    }
}



void AssemblyPath::clear()
{
    segments.clear();
    assembledSegments.clear();
    linksRleSequence.clear();
    linksRepeatCounts.clear();
    skipAtSegmentBegin.clear();
    skipAtSegmentEnd.clear();
}



// Assemble each segment on the path.
void AssemblyPath::assembleSegments(const AssemblyGraph& assemblyGraph)
{
    assembledSegments.clear();
    assembledSegments.resize(segments.size());
    for(uint64_t i=0; i<segments.size(); i++) {
        const uint64_t segmentId = segments[i].first;
        AssembledSegment& assembledSegment = assembledSegments[i];
        assembleMarkerGraphPath(
            assemblyGraph.readRepresentation,
            assemblyGraph.k,
            assemblyGraph.markers,
            assemblyGraph.markerGraph,
            assemblyGraph.markerGraphPaths[segmentId], false, assembledSegment);
    }
}



void AssemblyPath::writeSegmentSequences()
{
    ofstream fasta("PathSegmentsSequence.fasta");
    ofstream txt("PathSegmentsRleSequence.txt");

    for(uint64_t i=0; i<segments.size(); i++) {
        const uint64_t segmentId = segments[i].first;
        AssembledSegment& assembledSegment = assembledSegments[i];

        if(skipAtSegmentBegin[i] + skipAtSegmentEnd[i] > assembledSegment.runLengthSequence.size()) {
            cout << "Segment " << segmentId <<
                " has overlapping skips on left/right." << endl;
            continue;
        }
        if(skipAtSegmentBegin[i] + skipAtSegmentEnd[i] == assembledSegment.runLengthSequence.size()) {
            cout << "Segment " << segmentId <<
                " is skipped entirely." << endl;
            continue;
        }

        const uint64_t rleLengthKept =
            assembledSegment.runLengthSequence.size() -
            skipAtSegmentBegin[i] -
            skipAtSegmentEnd[i];

        // Write out the RLE sequence.
        // Also compute the raw length kept.
        uint64_t rawLengthKept = 0;
        txt << "S" << i << " " << segmentId << "\n";
        copy(
            assembledSegment.runLengthSequence.begin() +  skipAtSegmentBegin[i],
            assembledSegment.runLengthSequence.end() - skipAtSegmentEnd[i],
            ostream_iterator<Base>(txt));
        txt << "\n";
        for(uint64_t j=skipAtSegmentBegin[i]; j<rleLengthKept; j++) {
            const uint32_t r = assembledSegment.repeatCounts[j];
            rawLengthKept += r;
            if(r < 10) {
                txt << r;
            } else {
                txt << "*";
            }
        }
        txt << "\n";

        fasta <<
            ">S" << i <<
            " segment " << segmentId <<
            ", length " << rawLengthKept << "\n";
        for(uint64_t j=skipAtSegmentBegin[i]; j<rleLengthKept; j++) {
            const Base b = assembledSegment.runLengthSequence[j];
            const uint32_t r = assembledSegment.repeatCounts[j];
            for(uint64_t k=0; k<r; k++) {
                fasta << b;
            }
        }
        fasta << "\n";
    }
}



void AssemblyPath::writeLinkSequences(const AssemblyGraph& assemblyGraph)
{
    ofstream fasta("PathLinksSequence.fasta");
    ofstream txt("PathLinksRleSequence.txt");

    for(uint64_t i=0; i<segments.size()-1; i++) {
        const uint64_t segmentId0 = segments[i].first;
        const uint64_t segmentId1 = segments[i+1].first;
        const uint64_t linkId = assemblyGraph.findLink(segmentId0, segmentId1);
        const vector<Base>& rleSequence = linksRleSequence[i];
        const vector<uint64_t>& repeatCounts = linksRepeatCounts[i];
        SHASTA_ASSERT(rleSequence.size() == repeatCounts.size());
        if(rleSequence.empty()) {
            continue;
        }

        fasta <<
            ">L" << i <<
            " link " << linkId << " " << segmentId0 << "->"<< segmentId1 << "\n";
        for(uint64_t j=0; j<rleSequence.size(); j++) {
            const Base b = rleSequence[j];
            const uint64_t repeatCount = repeatCounts[j];
            for(uint64_t k=0; k<repeatCount; k++) {
                fasta << b;
            }
        }
        fasta << "\n";

        txt << "L" << i <<
            " link " << linkId << " " << segmentId0 << "->"<< segmentId1 << "\n";
        copy(rleSequence.begin(), rleSequence.end(),
            ostream_iterator<Base>(txt));
        txt << "\n";
        for(const uint64_t r: repeatCounts) {
            if(r < 10) {
                txt << r;
            } else {
                txt << "*";
            }
        }
        txt << "\n";
    }
}


// Compute consensus sequence for Link, given sequences of
// the oriented reads, which must all be anchored on both sides.
void AssemblyPath::computeLinkConsensusUsingSpoa(
    const vector<OrientedReadId> orientedReadIds,
    const vector< vector<Base> > rleSequences,
    const vector< vector<uint64_t> > repeatCounts,
    uint64_t readRepresentation,
    const ConsensusCaller& consensusCaller,
    bool debug,
    ostream& html,
    vector<Base>& consensusRleSequence,
    vector<uint64_t>& consensusRepeatCounts
    ) const
{
    SHASTA_ASSERT(rleSequences.size() == orientedReadIds.size());
    SHASTA_ASSERT(repeatCounts.size() == orientedReadIds.size());

    // Create the spoa alignment engine and elignment graph.
    const spoa::AlignmentType alignmentType = spoa::AlignmentType::kNW;
    const int8_t match = 1;
    const int8_t mismatch = -1;
    const int8_t gap = -1;
    auto spoaAlignmentEngine = spoa::AlignmentEngine::Create(alignmentType, match, mismatch, gap);
    spoa::Graph spoaAlignmentGraph;

    // Add the oriented read sequences to the alignment.
    string sequenceString;
    for(const vector<Base>& sequence: rleSequences) {

        // Add it to the alignment.
        sequenceString.clear();
        for(const Base base: sequence) {
            sequenceString += base.character();
        }
        auto alignment = spoaAlignmentEngine->Align(sequenceString, spoaAlignmentGraph);
        spoaAlignmentGraph.AddAlignment(alignment, sequenceString);
    }

    // Compute the multiple sequence alignment.
    const vector<string> msa = spoaAlignmentGraph.GenerateMultipleSequenceAlignment();
    const string consensus = spoaAlignmentGraph.GenerateConsensus();
    const uint64_t msaLength = msa.front().size();
    if(debug) {
        cout << "Multiple sequence alignment has length " << msaLength << ":" << endl;
        for(const string& s: msa) {
            cout << s << endl;
        }
    }


    // Compute coverage for each base at each position of the MSA.
    // Use position 4 for gaps.
    vector<Coverage> coverage(msaLength);
    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const vector<Base>& rleSequence = rleSequences[i];
        const vector<uint64_t>& repeatCount = repeatCounts[i];
        const string& msaString = msa[i];

        // Here:
        // rPosition = position in rle sequence of oriented read.
        // aPosition = position in alignment
        uint64_t rPosition = 0;
        for(uint64_t aPosition=0; aPosition<msaLength; aPosition++) {
            const AlignedBase alignedBase = AlignedBase::fromCharacter(msaString[aPosition]);
            if(alignedBase.isGap()) {
                coverage[aPosition].addRead(alignedBase, orientedReadId.getStrand(), 0);
            } else {
                SHASTA_ASSERT(AlignedBase(rleSequence[rPosition]) == alignedBase);
                if(readRepresentation == 1) {
                    coverage[aPosition].addRead(
                        alignedBase,
                        orientedReadId.getStrand(),
                        repeatCount[rPosition]);
                } else {
                    coverage[aPosition].addRead(
                        alignedBase,
                        orientedReadId.getStrand(),
                        1);
                }
                ++rPosition;
            }
        }
        SHASTA_ASSERT(rPosition == rleSequence.size());
    }



    // Compute consensus base and repeat count at every position in the alignment.
    vector<AlignedBase> msaConsensusSequence(msaLength);
    vector<uint64_t> msaConsensusRepeatCount(msaLength);
    for(uint64_t aPosition=0; aPosition<msaLength; aPosition++) {
        const Consensus consensus = consensusCaller(coverage[aPosition]);
        msaConsensusSequence[aPosition] = consensus.base;
        msaConsensusRepeatCount[aPosition] = consensus.repeatCount;
    }



    // Fill in the output arguments.
    // These are the same as msaConsensusSequence and msaConsensusRepeatCount,
    // but with the gap bases removed.
    consensusRleSequence.clear();
    consensusRepeatCounts.clear();
    for(uint64_t aPosition=0; aPosition<msaLength; aPosition++) {
        const AlignedBase alignedBase = msaConsensusSequence[aPosition];
        if(not alignedBase.isGap()) {
            consensusRleSequence.push_back(Base(alignedBase));
            consensusRepeatCounts.push_back(msaConsensusRepeatCount[aPosition]);
        }
    }



    // Html output of the alignment.
    if(html.good()) {
        html << "Coverage " << rleSequences.size() << "<br>\n";
        html << "Alignment length " << msaLength << "<br>\n";
        html << "<div style='font-family:monospace;white-space:nowrap;'>\n";
        for(uint64_t i=0; i<orientedReadIds.size(); i++) {
            const OrientedReadId orientedReadId = orientedReadIds[i];
            const string& msaString = msa[i];

            for(const char c: msaString) {
                const AlignedBase alignedBase = AlignedBase::fromCharacter(c);
                if(alignedBase.isGap()) {
                    html << alignedBase;
                } else {
                    html << "<span style='background-color:" << alignedBase.htmlColor() <<
                        "'>" << alignedBase << "</span>";
                }
            }

            // If using the RLE representation, also write the
            // repeat count at each position.
            if(readRepresentation == 1) {
                const vector<Base>& rleSequence = rleSequences[i];
                const vector<uint64_t>& repeatCount = repeatCounts[i];

                // Here:
                // rPosition = position in RLE sequence of oriented read.
                // aPosition = position in alignment
                uint64_t rPosition = 0;
                html << "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
                for(uint64_t aPosition=0; aPosition<msaLength; aPosition++) {
                    const AlignedBase alignedBase = AlignedBase::fromCharacter(msaString[aPosition]);
                    if(alignedBase.isGap()) {
                        html << alignedBase;
                    } else {
                        SHASTA_ASSERT(AlignedBase(rleSequence[rPosition]) == alignedBase);
                        const uint64_t r = repeatCount[rPosition];
                        html << "<span style='background-color:" << alignedBase.htmlColor() <<
                            "'>";
                        if(r < 10) {
                            html << r;
                        } else {
                            html << "*";
                        }
                        html << "</span>";
                        ++rPosition;
                    }
                }
                SHASTA_ASSERT(rPosition == rleSequence.size());
            }

            html << " " << orientedReadId << "<br>\n";
        }



       // Also write the consensus.
        html << "<br>\n";
        for(uint64_t aPosition=0; aPosition<msaLength; aPosition++) {
            const AlignedBase alignedBase = msaConsensusSequence[aPosition];
            if(alignedBase.isGap()) {
                html << alignedBase;
            } else {
                html << "<span style='background-color:" << alignedBase.htmlColor() <<
                    "'>" << alignedBase << "</span>";
            }
        }
        if(readRepresentation == 1) {
            html << "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
            for(uint64_t aPosition=0; aPosition<msaLength; aPosition++) {
                const AlignedBase alignedBase = msaConsensusSequence[aPosition];
                if(alignedBase.isGap()) {
                    html << alignedBase;
                } else {
                    const uint64_t r = msaConsensusRepeatCount[aPosition];
                    html << "<span style='background-color:" << alignedBase.htmlColor() <<
                        "'>";
                    if(r < 10) {
                        html << r;
                    } else {
                        html << "*";
                    }
                    html << "</span>";
                }
            }
        }
        html << " Consensus<br>\n";
        html << "</div>\n";

        html << "<h3>Consensus</h3>";
        html << "<div style='font-family:monospace;white-space:nowrap;'>\n";
        for(const Base b: consensusRleSequence) {
            html << b;
        }
        html << "<br>\n";
        for(const uint64_t r: consensusRepeatCounts) {
            if(r < 10) {
                html << r;
            } else {
                html << "*";
            }
        }
        html << "<br>\n";
        html << "<br>\n";
        for(uint64_t i=0; i<consensusRleSequence.size(); i++) {
            const Base b = consensusRleSequence[i];
            const uint64_t r = consensusRepeatCounts[i];
            for(uint64_t j=0; j<r; j++) {
                html << b;
            }
        }
        html << "<br>\n";
        html << "</div>\n";
    }
}



// Final assembly of segments and links sequence into the path sequence.
void AssemblyPath::assemble()
{
    rleSequence.clear();
    repeatCounts.clear();
    rawSequence.clear();

    // Assemble RLE sequence.
    for(uint64_t i=0; i<segments.size(); i++) {
        const uint64_t segmentId = segments[i].first;
        AssembledSegment& assembledSegment = assembledSegments[i];

        if(skipAtSegmentBegin[i] + skipAtSegmentEnd[i] > assembledSegment.runLengthSequence.size()) {
            cout << "Segment " << segmentId <<
                " has overlapping skips on left/right." << endl;
            continue;
        }
        if(skipAtSegmentBegin[i] + skipAtSegmentEnd[i] == assembledSegment.runLengthSequence.size()) {
            cout << "Segment " << segmentId <<
                " is skipped entirely." << endl;
            continue;
        }

        // Add the sequence of this segment.
        copy(
            assembledSegment.runLengthSequence.begin() +  skipAtSegmentBegin[i],
            assembledSegment.runLengthSequence.end() - skipAtSegmentEnd[i],
            back_inserter(rleSequence));
        copy(
            assembledSegment.repeatCounts.begin() +  skipAtSegmentBegin[i],
            assembledSegment.repeatCounts.end() - skipAtSegmentEnd[i],
            back_inserter(repeatCounts));

        // Add the sequence of the link following this segment.
        if(i != segments.size() - 1) {
            copy(
                linksRleSequence[i].begin(),  linksRleSequence[i].end(),
                back_inserter(rleSequence));
            copy(
                linksRepeatCounts[i].begin(),  linksRepeatCounts[i].end(),
                back_inserter(repeatCounts));
        }
    }
    SHASTA_ASSERT(rleSequence.size() == repeatCounts.size());



    // Now we can compute the raw sequence.
    for(uint64_t i=0; i<rleSequence.size(); i++) {
        const Base b = rleSequence[i];
        const uint64_t r = repeatCounts[i];
        for(uint64_t k=0; k<r; k++) {
            rawSequence.push_back(b);
        }
    }

    ofstream fasta("PathSequence.fasta");
    fasta << ">Path" << endl;
    copy(rawSequence.begin(), rawSequence.end(), ostream_iterator<Base>(fasta));
    fasta << "\n";
}
