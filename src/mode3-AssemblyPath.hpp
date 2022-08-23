#ifndef SHASTA_MODE3_ASSEMBLY_PATH_HPP
#define SHASTA_MODE3_ASSEMBLY_PATH_HPP

// Shasta.
#include "AssembledSegment.hpp"

// Standard library.
#include "cstdint.hpp"
#include "span.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    namespace mode3 {
        class AssemblyPath;
        class AssemblyPathLink;
        class AssemblyPathSegment;

        class AssemblyGraph;
    }

    class Base;
    class ConsensusCaller;
    class OrientedReadId;
}



// A segment in an AssemblyPath.
class shasta::mode3::AssemblyPathSegment {
public:

    // The id of this segment, in the AssemblyGraph.
    uint64_t id;

    // Each primary segment in the path has high Jaccard similarity
    // with the previous primary segment.
    // The first and last segment are always primary segments.
    bool isPrimary;

    // The AssembledSegment contains the sequence for this segment
    // plus information on how the sequence was extracted from the
    // marker graph.
    // The sequence includes the first and last marker graph vertex
    // of this segment.
    AssembledSegment assembledSegment;

    // For assembly of the path sequence, we don't use the entire
    // sequence of the AssembledSegment.
    // We trim some bases at each end to avoid overlap
    // with adjacent segments and links.
    // When a segment is adjacent to a non-trivial link,
    // we give priority to link sequence over segment sequence.
    // The reason is that sequence assembled from links
    // is generally more accurate because it is assembled
    // using only a restricted set of
    // oriented reads that are believed to originate from the
    // sequence copy we are assembling.
    uint64_t leftTrim;
    uint64_t rightTrim;
    span<const Base> trimmedRleSequence() const
    {
        return span<const Base>(
            assembledSegment.runLengthSequence.begin() + leftTrim,
            assembledSegment.runLengthSequence.end() - rightTrim);
    }
    span<const uint32_t> trimmedRepeatCounts() const
    {
        return span<const uint32_t>(
            assembledSegment.repeatCounts.begin() + leftTrim,
            assembledSegment.repeatCounts.end() - rightTrim);
    }

    // Constructor.
    AssemblyPathSegment(uint64_t id, bool isPrimary);
};



// A link in an AssemblyPath.
class shasta::mode3::AssemblyPathLink {
public:

    // A link is trivial if the last marker graph vertex
    // of the source segment coincides with the first marker
    // graph vertex of the target segment.
    // In this case the link does not need to be assembled
    // and all the next fields are left empty.
    bool isTrivial;

    // The RLE sequence as computed by the MSA
    // of oriented reads in the link.
    // This overlaps with adjacent segments.
    vector<Base> msaRleSequence;
    vector<uint64_t> msaRepeatCounts;

    // The trimmed RLE sequence is obtained from
    // the MSA sequence by removing bases at the two ends
    // that are identical with the adjacent segments.
    vector<Base> trimmedRleSequence;
    vector<uint64_t> trimmedRepeatCounts;
};



// An assembly path in the mode3::AssemblyGraph
class shasta::mode3::AssemblyPath {
public:

    // The segments and links on the path.
    vector<AssemblyPathSegment> segments;
    vector<AssemblyPathLink> links;

    // Top level function to assemble sequence for this path.
    void assemble(const AssemblyGraph&);

    // Assemble the sequence of each segment.
    void assembleSegments(const AssemblyGraph&);
    void writeSegmentSequences();

    // Assemble links in this assembly path.
    void assembleLinks(const AssemblyGraph&);
    void writeLinkSequences(const AssemblyGraph&);

    // Final assembly of segments and links sequence into the path sequence.
    void assemble();
    vector<Base> rleSequence;
    vector<uint64_t> repeatCounts;
    vector<Base> rawSequence;

    void clear();

    // Use spoa to compute consensus sequence for a link, given sequences of
    // the oriented reads, which must all be anchored on both sides.
    void computeLinkConsensusUsingSpoa(
        const vector<OrientedReadId> orientedReadIds,
        const vector< vector<Base> > rleSequences,
        const vector< vector<uint64_t> > repeatCounts,
        uint64_t readRepresentation,
        const ConsensusCaller&,
        bool debug,
        ostream& html,
        vector<Base>& consensusRleSequence,
        vector<uint64_t>& consensusRepeatCounts
        ) const;

    // Find the oriented reads to be used to assemble
    // links between the segment at position0
    // in the assembly path (which must be a reference segment)
    // and the next reference segment in the path.
    // The oriented reads are returned sorted.
    void findOrientedReadsForLinks(
        uint64_t position0,
        const AssemblyGraph&,
        vector<OrientedReadId>&) const;

};


#endif

