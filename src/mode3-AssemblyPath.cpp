// Shasta.
#include "mode3-AssemblyPath.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "ConsensusCaller.hpp"
#include "deduplicate.hpp"
#include "html.hpp"
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
    initializeLinks(assemblyGraph);
    assembleLinks(assemblyGraph);

    writeSegmentSequences();
    writeLinkSequences(assemblyGraph);

    assemble();
}

// Initialize the links.
// This only resizes the links vector and fills in the id and isTrivial
// fields of each link.
void AssemblyPath::initializeLinks(const AssemblyGraph& assemblyGraph)
{
    SHASTA_ASSERT(segments.size() > 1);
    links.resize(segments.size()-1);

    // Fill in the id and isTrivial fields of each link.
    for(uint64_t position0=0; position0<links.size(); position0++) {
        const uint64_t position1 = position0 + 1;

        // Access the source and target segments of this link.
        // We will process the link between segmentId0 and segmentId1.
        AssemblyPathSegment& segment0 = segments[position0];
        AssemblyPathSegment& segment1 = segments[position1];

        // Fill in the id and isTrivial fields.
        AssemblyPathLink& assemblyPathLink = links[position0];
        assemblyPathLink.id = assemblyGraph.findLink(segment0.id, segment1.id);
        const AssemblyGraph::Link& link = assemblyGraph.links[assemblyPathLink.id];
        assemblyPathLink.isTrivial = link.segmentsAreAdjacent;

        SHASTA_ASSERT(segment0.id == link.segmentId0);
        SHASTA_ASSERT(segment1.id == link.segmentId1);
    }


    // Fill in the previousPrimarySegmentId field of each link.
    SHASTA_ASSERT(segments.front().isPrimary);
    uint64_t lastPrimarySegmentSeen = invalid<uint64_t>;
    for(uint64_t position=0; position<links.size(); position++) {
        const AssemblyPathSegment& segment = segments[position];
        if(segment.isPrimary) {
            lastPrimarySegmentSeen = segment.id;
        }
        links[position].previousPrimarySegmentId = lastPrimarySegmentSeen;
    }



    // Fill in the nextPrimarySegmentId field of each link.
    SHASTA_ASSERT(segments.back().isPrimary);
    lastPrimarySegmentSeen = invalid<uint64_t>;
    for(uint64_t position = links.size() - 1; /* Check later */; position--) {
        const AssemblyPathSegment& segment = segments[position + 1];
        if(segment.isPrimary) {
            lastPrimarySegmentSeen = segment.id;
        }
        links[position].nextPrimarySegmentId = lastPrimarySegmentSeen;

        if(position == 0) {
            break;
        }
    }
}



// Assemble links in this assembly path.
void AssemblyPath::assembleLinks(const AssemblyGraph& assemblyGraph)
{
    SHASTA_ASSERT((assemblyGraph.k % 2) == 0);

    // Don't skip any bases at the beginning of the first
    // segment and at the end of the last segment.
    segments.front().leftTrim = 0;
    segments.back().rightTrim = 0;

    ofstream html("Msa.html");

    // Loop over links in the path.
    links.resize(segments.size()-1);
    for(uint64_t position0=0; position0<links.size(); position0++) {
        assembleLinkAtPosition(assemblyGraph, position0, html);
    }
}



void AssemblyPath::assembleLinkAtPosition(
    const AssemblyGraph& assemblyGraph,
    uint64_t position0,
    ostream& html)
{
    const bool debug = false;

    AssemblyPathLink& link = links[position0];
    const uint64_t position1 = position0 + 1;

    // Access the source and target segments of this link.
    // We will process the link between segmentId0 and segmentId1.
    AssemblyPathSegment& segment0 = segments[position0];
    AssemblyPathSegment& segment1 = segments[position1];

    if(debug) {
        cout << "Assembling link " << link.id << " " << segment0.id << "->" << segment1.id <<
            " at position " << position0 << " in the assembly path." << endl;
    }

    if(link.isTrivial) {

        // The two segments are consecutive in the marker graph.
        // This is a trivial link because the two segments share a terminal
        // marker graph vertex.
        // Just trim from the assembly the last k/2 RLE bases of segmentId0
        // and the first k/2 RLE bases of segmentId1.
        assembleTrivialLink(segment0, segment1, link, assemblyGraph.k);

    } else {

        assembleNonTrivialLink(
            assemblyGraph,
            segment0,
            segment1,
            link,
            html);
    }
}



void AssemblyPath::assembleNonTrivialLink(
    const AssemblyGraph& assemblyGraph,
    AssemblyPathSegment& segment0,
    AssemblyPathSegment& segment1,
    AssemblyPathLink& link,
    ostream& html)
{
    const bool debug = false;


    // First, find:
    // - The position in segmentId0 of the leftmost transition.
    // - The position in segmentId1 of the rightmost transition.
    uint64_t minEdgePosition0 = assemblyGraph.markerGraphPaths[segment0.id].size();
    uint64_t maxEdgePosition1 = 0;
    for(const auto& p: assemblyGraph.transitions[link.id]) {
        const OrientedReadId orientedReadId = p.first;

        // If not in previousPrimarySegmentId or nextPrimarySegmentId, skip it.
        if(not(
            assemblyGraph.segmentContainsOrientedRead(link.previousPrimarySegmentId, orientedReadId)
            or
            assemblyGraph.segmentContainsOrientedRead(link.nextPrimarySegmentId, orientedReadId)
            )) {
            continue;
        }

        // Access the transition from segmentId0 to segmentId1 for this oriented read.
        const Transition& transition = p.second;

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

    // To compute an MSA anchored at both sides,we will extend the
    // sequence of each read to the left/right using the sequence of
    // adjacent segments.
    const AssembledSegment& assembledSegment0 = segment0.assembledSegment;
    SHASTA_ASSERT(not assembledSegment0.runLengthSequence.empty());
    const AssembledSegment& assembledSegment1 = segment1.assembledSegment;
    SHASTA_ASSERT(not assembledSegment1.runLengthSequence.empty());


    // Now extract the portion of each oriented read sequence that
    // will be used to assemble this link.
    vector<OrientedReadId> orientedReadIdsForAssembly;
    vector< vector<Base> > orientedReadsSequencesForAssembly;
    vector< vector<uint32_t> > orientedReadsRepeatCountsForAssembly;
    for(const auto& p: assemblyGraph.transitions[link.id]) {
        const OrientedReadId orientedReadId = p.first;

        // If not in previousPrimarySegmentId or nextPrimarySegmentId, skip it.
        if(not(
            assemblyGraph.segmentContainsOrientedRead(link.previousPrimarySegmentId, orientedReadId)
            or
            assemblyGraph.segmentContainsOrientedRead(link.nextPrimarySegmentId, orientedReadId)
            )) {
            continue;
        }

        // Access the transition from segmentId0 to segmentId1 for this oriented read.
        const Transition& transition = p.second;

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
        vector<uint32_t> orientedReadExtendedRepeatCounts;
        const auto addToExtendedSequence = back_inserter(orientedReadExtendedSequence);
        copy(leftSequence, addToExtendedSequence);
        copy(orientedReadSequence, addToExtendedSequence);
        copy(rightSequence, addToExtendedSequence);
        const auto addToRepeatCounts = back_inserter(orientedReadExtendedRepeatCounts);
        copy(leftRepeatCounts, addToRepeatCounts);
        copy(orientedReadRepeatCounts, addToRepeatCounts);
        copy(rightRepeatCounts, addToRepeatCounts);

        orientedReadIdsForAssembly.push_back(orientedReadId);
        orientedReadsSequencesForAssembly.push_back(orientedReadExtendedSequence);
        orientedReadsRepeatCountsForAssembly.push_back(orientedReadExtendedRepeatCounts);

        if(debug) {
            copy(orientedReadExtendedSequence, ostream_iterator<Base>(cout));
            cout << " " << orientedReadId << endl;
        }
    }

    // Compute the consensus sequence for the link.
    html<< "<h2>Link " << link.id << "</h2>\n";
    computeLinkConsensusUsingSpoa(
        orientedReadIdsForAssembly,
        orientedReadsSequencesForAssembly,
        orientedReadsRepeatCountsForAssembly,
        assemblyGraph.readRepresentation,
        assemblyGraph.consensusCaller,
        debug,
        html,
        link.msaRleSequence,
        link.msaRepeatCounts
        );
    SHASTA_ASSERT(link.msaRleSequence.size() == link.msaRepeatCounts.size());

    if(debug) {
        cout << "Consensus RLE sequence length before trimming " << link.msaRleSequence.size() << endl;
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
    for(uint64_t i=begin0; (i!=end0 and (i-begin0)<link.msaRleSequence.size()); i++) {
        if(link.msaRleSequence[i-begin0] == assembledSegment0.runLengthSequence[i]) {
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
        const uint64_t j = link.msaRleSequence.size() - (end1 - i);
        if(link.msaRleSequence[j] == assembledSegment1.runLengthSequence[i]) {
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
    identicalOnRight = min(identicalOnRight, link.msaRleSequence.size()-identicalOnLeft);
    if(debug) {
        cout << "Identical on right: " << identicalOnRight << endl;
    }

    // Trim these identical bases from the link consensus sequence.
    link.leftTrim = identicalOnLeft;
    link.rightTrim = identicalOnRight;

    // Compute and store the number of bases to be trimmed at the end of segmentId0
    // and at the beginning of segmentId1.
    segment0.rightTrim =
        assembledSegment0.runLengthSequence.size() -
        assembledSegment0.vertexOffsets[minVertexPosition0] -
        identicalOnLeft;
    segment1.leftTrim =
        assembledSegment1.vertexOffsets[maxVertexPosition1] + assemblyGraph.k
        - identicalOnRight;
}



void AssemblyPath::assembleTrivialLink(
    AssemblyPathSegment& segment0,
    AssemblyPathSegment& segment1,
    AssemblyPathLink& link,
    uint64_t k)
{
    SHASTA_ASSERT(link.isTrivial);
    SHASTA_ASSERT(link.msaRleSequence.empty());
    SHASTA_ASSERT(link.msaRepeatCounts.empty());
    SHASTA_ASSERT(link.leftTrim == 0);
    SHASTA_ASSERT(link.rightTrim == 0);

    // Just trim k/2 bases from the adjacent segments,
    // because they are adjacent in the marker graph.
    segment0.rightTrim = k/2;
    segment1.leftTrim  = k/2;
}



void AssemblyPath::clear()
{
    segments.clear();
    links.clear();
}



// Assemble each segment on the path.
void AssemblyPath::assembleSegments(const AssemblyGraph& assemblyGraph)
{
    for(uint64_t i=0; i<segments.size(); i++) {
        AssemblyPathSegment& segment = segments[i];
        assembleMarkerGraphPath(
            assemblyGraph.readRepresentation,
            assemblyGraph.k,
            assemblyGraph.markers,
            assemblyGraph.markerGraph,
            assemblyGraph.markerGraphPaths[segment.id],
            false,
            segment.assembledSegment);
    }
}



void AssemblyPath::writeSegmentSequences()
{
    ofstream fasta("PathSegmentsSequence.fasta");
    ofstream txt("PathSegmentsRleSequence.txt");

    for(uint64_t i=0; i<segments.size(); i++) {
        const AssemblyPathSegment& segment = segments[i];
        const uint64_t segmentId = segment.id;
        const AssembledSegment& assembledSegment = segment.assembledSegment;

        if(segment.leftTrim + segment.rightTrim > assembledSegment.runLengthSequence.size()) {
            continue;
        }

        // Write the trimmed RLE sequence to txt.
        const auto trimmedRleSequence = segment.trimmedRleSequence();
        const auto trimmedRepeatCounts = segment.trimmedRepeatCounts();
        txt << "S" << i << " " << segmentId << "\n";
        copy( trimmedRleSequence, ostream_iterator<Base>(txt));
        txt << "\n";
        for(const uint32_t r: trimmedRepeatCounts) {
            txt << repeatCountCharacter(r);
        }
        txt << "\n";

        // Write the trimmed raw sequence to fasta.
        vector<Base> trimmedRawSequence;
        segment.getTrimmedRawSequence(trimmedRawSequence);
        fasta <<
            ">S" << i <<
            " segment " << segmentId <<
            ", length " << trimmedRawSequence.size() << "\n";
        copy(trimmedRawSequence, ostream_iterator<Base>(fasta));
        fasta << "\n";
    }
}



void AssemblyPath::writeLinkSequences(const AssemblyGraph& assemblyGraph)
{
    ofstream fasta("PathLinksSequence.fasta");
    ofstream txt("PathLinksRleSequence.txt");

    for(uint64_t i=0; i<segments.size()-1; i++) {
        const uint64_t segmentId0 = segments[i].id;
        const uint64_t segmentId1 = segments[i+1].id;
        const uint64_t linkId = assemblyGraph.findLink(segmentId0, segmentId1);
        const span<const Base> rleSequence = links[i].trimmedRleSequence();
        const span<const uint32_t> repeatCounts = links[i].trimmedRepeatCounts();
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
        copy(rleSequence, ostream_iterator<Base>(txt));
        txt << "\n";
        for(const uint32_t r: repeatCounts) {
            txt << repeatCountCharacter(r);
        }
        txt << "\n";
    }
}



// Compute consensus sequence for Link, given sequences of
// the oriented reads, which must all be anchored on both sides.// Lower level version.
void AssemblyPath::computeLinkConsensusUsingSpoa(
    const vector<OrientedReadId> orientedReadIds,
    const vector< vector<Base> > rleSequences,
    const vector< vector<uint32_t> > repeatCounts,
    uint64_t readRepresentation,
    const ConsensusCaller& consensusCaller,
    bool debug,
    ostream& html,
    vector<Base>& consensusRleSequence,
    vector<uint32_t>& consensusRepeatCounts
    )
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
        const vector<uint32_t>& repeatCount = repeatCounts[i];
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
    vector<uint32_t> msaConsensusRepeatCount(msaLength);
    vector<uint64_t> msaConsensusDiscordantCount(msaLength);
    for(uint64_t aPosition=0; aPosition<msaLength; aPosition++) {
        const Coverage& c = coverage[aPosition];
        const Consensus consensus = consensusCaller(c);
        msaConsensusSequence[aPosition] = consensus.base;
        msaConsensusRepeatCount[aPosition] = uint32_t(consensus.repeatCount);

        // Compute discordant count at this position of the alignment.
        msaConsensusDiscordantCount[aPosition] = 0;
        for(uint64_t b=0; b<5; b++) {
            if(b != consensus.base.value) {
                msaConsensusDiscordantCount[aPosition] += c.coverage(AlignedBase::fromInteger(b));
            }
        }
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
                const vector<uint32_t>& repeatCount = repeatCounts[i];

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

        // Write the discordant count.
        for(uint64_t aPosition=0; aPosition<msaLength; aPosition++) {
            const uint64_t d = msaConsensusDiscordantCount[aPosition];
            const double errorRate = double(d) / double(orientedReadIds.size());
            int hue;
            if(errorRate < .01) {
                hue = 120;    // Q>=20, green.
            } else {
                const double Q = -10. * log10(errorRate);
                hue = int(std::round(6. * Q));  // 60 at Q=10 (yellow), 120 at Q=20 (green).
            }
            const string color = "hsl(" + to_string(hue) + ",100%, 70%)";
            html << "<span style='background-color:" << color << "'>";
            html << repeatCountCharacter(uint32_t(d));
            html << "</span>";
        }
        html << "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
        for(uint64_t aPosition=0; aPosition<msaLength; aPosition++) {
            html << "&nbsp;";
        }
        html << " Discordant<br>\n";

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
        AssemblyPathSegment& segment = segments[i];
        const AssembledSegment& assembledSegment = segment.assembledSegment;
        segment.rlePosition = rleSequence.size();
        segment.rawPosition = rawSequence.size();


        if(segment.leftTrim + segment.rightTrim > assembledSegment.runLengthSequence.size()) {
            // The left and right trim of this segment overlap.
            // To handle this case, just take the excess number of bases out of the sequence
            // we already assembled.
            // This is not a great solution, but better than nothing.
            const uint64_t excessTrim =
               (segment.leftTrim + segment.rightTrim) - assembledSegment.runLengthSequence.size();
            SHASTA_ASSERT(excessTrim <= rleSequence.size());
            SHASTA_ASSERT(repeatCounts.size() == rleSequence.size());

            // Compute the excess trim in the raw sequence.
            uint64_t excessTrimRaw = 0;
            for(uint64_t i=0; i<excessTrim; i++) {
                excessTrimRaw += repeatCounts[repeatCounts.size() - 1 - i];
            }
            SHASTA_ASSERT(excessTrimRaw <= rawSequence.size());

            // Remove the excess trim from the sequence we already assembled.
            rleSequence.resize(rleSequence.size() - excessTrim);
            repeatCounts.resize(repeatCounts.size() - excessTrim);
            rawSequence.resize(rawSequence.size() - excessTrimRaw);
        } else {

            // This is the normal case.

            // Add the RLE sequence of this segment.
            const auto segmentTrimmedRleSequence = segment.trimmedRleSequence();
            const auto segmentTrimmedRepeatCounts = segment.trimmedRepeatCounts();
            copy(segmentTrimmedRleSequence, back_inserter(rleSequence));
            copy(segmentTrimmedRepeatCounts, back_inserter(repeatCounts));

            // Add the raw sequence of this segment.
            for(uint64_t i=0; i<segmentTrimmedRleSequence.size(); i++) {
                const Base b = segmentTrimmedRleSequence[i];
                const uint64_t r = segmentTrimmedRepeatCounts[i];
                for(uint64_t k=0; k<r; k++) {
                    rawSequence.push_back(b);
                }
            }
        }



        // Add the sequence of the link following this segment.
        if(i != segments.size() - 1) {
            AssemblyPathLink& link = links[i];
            link.rlePosition = rleSequence.size();
            link.rawPosition = rawSequence.size();

            // Add the RLE sequence of this link.
            const auto trimmedRleSequence = link.trimmedRleSequence();
            const auto trimmedRepeatCounts = link.trimmedRepeatCounts();
            copy(trimmedRleSequence, back_inserter(rleSequence));
            copy(trimmedRepeatCounts, back_inserter(repeatCounts));

            // Add the raw sequence of this link.
            for(uint64_t i=0; i<trimmedRleSequence.size(); i++) {
                const Base b = trimmedRleSequence[i];
                const uint64_t r = trimmedRepeatCounts[i];
                for(uint64_t k=0; k<r; k++) {
                    rawSequence.push_back(b);
                }
            }
        }
    }
    SHASTA_ASSERT(rleSequence.size() == repeatCounts.size());

    // Write it out.
    ofstream fasta("PathSequence.fasta");
    fasta << ">Path" << endl;
    copy(rawSequence, ostream_iterator<Base>(fasta));
    fasta << "\n";
}




AssemblyPathSegment::AssemblyPathSegment(
    uint64_t id,
    bool isPrimary) :
    id(id),
    isPrimary(isPrimary)
    {}



span<const Base> AssemblyPathSegment::trimmedRleSequence() const
{
    const auto begin = assembledSegment.runLengthSequence.begin() + leftTrim;
    const auto end = assembledSegment.runLengthSequence.end() - rightTrim;
    SHASTA_ASSERT(begin <= end);
    return span<const Base>(begin, end);
}



span<const uint32_t> AssemblyPathSegment::trimmedRepeatCounts() const
{
    const auto begin = assembledSegment.repeatCounts.begin() + leftTrim;
    const auto end = assembledSegment.repeatCounts.end() - rightTrim;
    SHASTA_ASSERT(begin <= end);
    return span<const uint32_t>(begin, end);
}



span<const Base> AssemblyPathLink::trimmedRleSequence() const
{
    const auto begin = msaRleSequence.begin() + leftTrim;
    const auto end = msaRleSequence.end() - rightTrim;
    SHASTA_ASSERT(begin <= end);
    return span<const Base>(begin, end);
}



span<const uint32_t> AssemblyPathLink::trimmedRepeatCounts() const
{
    const auto begin = msaRepeatCounts.begin() + leftTrim;
    const auto end = msaRepeatCounts.end() - rightTrim;
    SHASTA_ASSERT(begin <= end);
    return span<const uint32_t>(begin, end);
}



void AssemblyPathSegment::getTrimmedRawSequence(vector<Base>& trimmedRawSequence) const
{

    // Get the trimed RLE sequence and repeat counts.
    const span<const Base> trimmedRleSequenceSpan = trimmedRleSequence();
    const span<const uint32_t> trimmedRepeatCountsSpan = trimmedRepeatCounts();
    SHASTA_ASSERT(trimmedRleSequenceSpan.size() == trimmedRepeatCountsSpan.size());

    // Construct the raw sequence.
    trimmedRawSequence.clear();
    for(uint64_t i=0; i<trimmedRleSequenceSpan.size(); i++) {
        const Base b = trimmedRleSequenceSpan[i];
        const uint32_t r = trimmedRepeatCountsSpan[i];
        for(uint64_t k=0; k<r; k++) {
            trimmedRawSequence.push_back(b);
        }
    }
}



// Return a character to represent a repeat count
// when writing out RLE sequence.
char AssemblyPath::repeatCountCharacter(uint32_t r) {
    if(r < 10) {
        return '0' + char(r);
    } else if(r < 36) {
        return 'A' + char(r - 10);
    } else {
        return '*';
    }
}


void AssemblyPath::writeHtml(ostream& html) const
{
    SHASTA_ASSERT(segments.size() > 1);
    SHASTA_ASSERT(links.size() == segments.size() - 1);

    writeHtmlSummary(html);
    writeSequenceDialog(html);
    writeHtmlDetail(html);
}



void AssemblyPath::writeHtmlSummary(ostream& html) const
{
    html <<
        "<table>" <<
        "<tr><th class=left>First segment id<td class=centered>" << segments.front().id <<
        "<tr><th class=left>Last segment id<td class=centered>" << segments.back().id <<
        "<tr><th class=left>Number of segments<td class=centered>" << segments.size() <<
        "<tr><th class=left>Number of links<td class=centered>" << segments.size() - 1 <<
        "<tr><th class=left>Length of RLE sequence assembled<td class=centered>" << rleSequence.size() <<
        "<tr><th class=left>Length of raw sequence assembled<td class=centered>" << rawSequence.size() <<
        "</table>";
}



void AssemblyPath::writeHtmlDetail(ostream& html) const
{
    // Table legend.
    html <<
        "<p>Hover on table headers for information on the meaning of each column.";



    // Table header.
    html <<
        "<p>"
        "<table style='table-layout:fixed;font-family:monospace;font-size:9'>"
        "<tr>"

        "<th title ='S (segment) or L (link). "
        "Primary segments have a light blue background. "
        "Trivial links have a grey background.'>"
        "<span class=rotated>Type"

        "<th title='Segment or link id'>"
        "<span class=rotated>Id"

        "<th title='The id of the last primary segment preceding a link.'>"
        "<span class=rotated>Previous<br>primary<br>segment"

        "<th title='The id of the first primary segment following a link.'>"
        "<span class=rotated>Next<br>primary<br>segment"

        "<th title='The position of the trimmed raw sequence of this segment or link "
        "in the raw assembled sequence of the path.'>"
        "<span class=rotated>Raw<br>position"

        "<th title='The complete raw sequence for this segment or link. "
        "The red portion is trimmed out and not used for assembly.'>"
        "Raw sequence"

        "<th title='Assembly details for non-trivial links.'>"
        "<span class=rotated>Detail";



    // Main body of the table.
    // There is one row for each segment and one row for each link.
    for(uint64_t position=0; position<segments.size(); position++) {

        // Write a row for the segment at this position.
        const AssemblyPathSegment& segment = segments[position];
        const AssembledSegment& assembledSegment = segment.assembledSegment;
        html << "<tr";
        if(segment.isPrimary) {
            html << " style='background-color:LightCyan'";
        }
        html <<
            ">"
            "<td class=centered>S"
            "<td class=centered>" << segment.id <<
            "<td><td>"
            "<td class=centered>" << segment.rawPosition;



        // Raw sequence for this segment.
        html << "<td class=centered style='max-width:300px;word-wrap:break-word'>";
         if(segment.leftTrim + segment.rightTrim > assembledSegment.runLengthSequence.size()) {

            // Exceptional case where the left and right trim overlap.
             html << "<span style='background-color:LightCoral'>";
             for(uint64_t i=0; i<assembledSegment.runLengthSequence.size(); i++) {
                 const Base b = assembledSegment.runLengthSequence[i];
                 const uint32_t r = assembledSegment.repeatCounts[i];
                 if(i == assembledSegment.runLengthSequence.size() - segment.rightTrim) {
                     html << "</span><span style='background-color:Fuchsia'>";
                 }
                 for(uint32_t k=0; k<r; k++) {
                     html << b;
                 }
                 if(i == segment.leftTrim - 1) {
                     html << "</span><span style='background-color:LightCoral'>";
                 }
             }
             html << "</span><td>";

        } else {

            // Normal case.
            html << "<span style='background-color:LightCoral'>";
            for(uint64_t i=0; i<assembledSegment.runLengthSequence.size(); i++) {
                const Base b = assembledSegment.runLengthSequence[i];
                const uint32_t r = assembledSegment.repeatCounts[i];
                if(i == segment.leftTrim) {
                    html << "</span>";
                }
                for(uint32_t k=0; k<r; k++) {
                    html << b;
                }
                if(i == assembledSegment.runLengthSequence.size() -1 - segment.rightTrim) {
                    html << "<span style='background-color:LightCoral'>";
                }

            }
            html << "</span><td>";
        }



        // Write a row for the link.
        if(position == links.size()) {
            break;
        }
        const AssemblyPathLink& link = links[position];
        html << "<tr";
        if(link.isTrivial) {
            html << " style='background-color:LightGray'";
        }
        html <<
            "><td class=centered>L" <<
            "<td class=centered>" << link.id <<
            "<td class=centered>" << link.previousPrimarySegmentId <<
            "<td class=centered>" << link.nextPrimarySegmentId <<
            "<td class=centered>" << link.rawPosition;

        // Raw sequence for this link.
        html << "<td class=centered style='max-width:300px;word-wrap:break-word'>";
        html << "<span style='background-color:LightCoral'>";
        for(uint64_t i=0; i<link.msaRleSequence.size(); i++) {
            const Base b = link.msaRleSequence[i];
            const uint32_t r = link.msaRepeatCounts[i];
            if(i == link.leftTrim) {
                html << "</span>";
            }
            for(uint32_t k=0; k<r; k++) {
                html << b;
            }
            if(i == link.msaRleSequence.size() -1 - link.rightTrim) {
                html << "<span style='background-color:LightCoral'>";
            }
        }
        html << "</span>";

        html << "<td class=centered>";
        if(not link.isTrivial) {
            html << "<a href='exploreMode3LinkAssembly?linkId=" << link.id <<
                "&previousPrimarySegmentId=" << link.previousPrimarySegmentId <<
                "&nextPrimarySegmentId=" << link.nextPrimarySegmentId <<
                "'>Detail</a>";
        }
    }

    // End the table.
    html << "</table>";

}



// This writes out a dialog that permit displaying
// selected portions of the path assembled sequence.
void AssemblyPath::writeSequenceDialog(ostream& html) const
{
    html << "<script>var assembledSequence = '";
    copy(rawSequence, ostream_iterator<Base>(html));
    html << "';</script>";

    html << R"zzz(
<form onsubmit="displaySequence(); return false;">
<br><input type=submit value='Display assembled sequence'>
 in the position range <input type=text id=begin>
 to <input type=text id=end>
</form>
<script>
function displaySequence()
{
    var beginString = document.getElementById('begin').value;
    var endString = document.getElementById('end').value;
    var begin = parseInt(beginString);
    var end = parseInt(endString);
    if((end < begin) || (end > assembledSequence.length)) {
        document.getElementById("assembledSequence").innerText = "";
    } else {
        document.getElementById("assembledSequence").innerText = assembledSequence.substring(begin, end);
    }
}
</script>
<p id=assembledSequence style='font-family:monospace;font-size:9pt;word-wrap:break-word;'>
    )zzz";
}
