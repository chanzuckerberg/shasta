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

        class AssemblyGraph;
    }

    class Base;
    class ConsensusCaller;
    class OrientedReadId;
}



// An assembly path in the mode3::AssemblyGraph
class shasta::mode3::AssemblyPath {
public:

    // The segments on the path.
    // The bool is true for reference segments.
    // The first and last segment are always reference segments.
    // A reference segment is one that is believed to be exclusive
    // to the sequence copy described by this path (that is,
    // it does not appear in other copies or haplotypes).
    vector< pair<uint64_t, bool> > segments;

    // Top level function to assemble sequence for this path.
    void assemble(const AssemblyGraph&);

    // Each segment gets assembled and the result stored here.
    vector<AssembledSegment> assembledSegments;
    void assembleSegments(const AssemblyGraph&);
    void writeAssembledSegments();

    // Assemble links in this assembly path.
    void assembleLinks(const AssemblyGraph&, bool debug);
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

};


#endif

