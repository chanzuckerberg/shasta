#ifndef SHASTA_MODE3_ASSEMBLY_PATH_HPP
#define SHASTA_MODE3_ASSEMBLY_PATH_HPP

// Shasta.
#include "AssembledSegment.hpp"

// Standard library.
#include "cstdint.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {
    namespace mode3 {
        class AssemblyPath;
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
public:

    // The id of this segment, in the AssemblyGraph.
    uint64_t id;

    // Each primary segment in the path has high Jaccard similarity
    // with the previous primary segment.
    // The first and last segment are always primary segments.
    bool isPrimary;

    // Constructor.
    AssemblyPathSegment(uint64_t id, bool isPrimary);
};




// An assembly path in the mode3::AssemblyGraph
class shasta::mode3::AssemblyPath {
public:

    // The segments on the path.
     vector<AssemblyPathSegment > segments;

    // Top level function to assemble sequence for this path.
    void assemble(const AssemblyGraph&);

    // Each segment gets assembled and the result stored here.
    vector<AssembledSegment> assembledSegments;
    void assembleSegments(const AssemblyGraph&);
    void writeSegmentSequences();

    // Assemble links in this assembly path.
    void assembleLinksOld(const AssemblyGraph&, bool debug);
    void assembleLinks(const AssemblyGraph&);
    void writeLinkSequences(const AssemblyGraph&);
    vector< vector<Base> > linksRleSequence;
    vector< vector<uint64_t> > linksRepeatCounts;

    // When assembling path sequence, we give priority to
    // sequence assembled from links over sequence assembled
    // from segments. The reason is that sequence assembled from links
    // is generally more accurate because it is assembled
    // using only the "reference oriented reads" - that is,
    // the oriented reads that are believed to originate from the
    // sequence copy we are assembling.
    // As a result, only part of the sequence of each assembled segment
    // is used when assembling path sequence.
    // Here we store the number of (RLE) bases to be skipped at the
    // beginning and end of each assembled segment.
    // These are computed in assembleLinks.
    vector<uint64_t> skipAtSegmentBegin;
    vector<uint64_t> skipAtSegmentEnd;

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

