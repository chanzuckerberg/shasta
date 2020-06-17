#ifndef SHASTA_SEQAN_HPP
#define SHASTA_SEQAN_HPP


// Wrappers to simplify SeqAn calls.

#include "SHASTA_ASSERT.hpp"

#include "utility.hpp"
#include "vector.hpp"

namespace shasta {

    // Align two integer sequences using SeqAn and
    // return the alignment score.
    // The alignment is returned as a vector of pairs of bools,
    // where false indicates a gap and true indicates a sequence element.
    // For example, the following alignment
    // ACTG-A
    // A-TGTA
    // would be returned as
    // {true, true},  A,A
    // {true, false}  C,-
    // {true, true},  T,T
    // {true, true}   G,G
    // {false, true}, -,T
    // {true, true}   A,A
    template<class Iterator>
        int64_t seqanAlign(
            Iterator begin0, Iterator end0,
            Iterator begin1, Iterator end1,
            int64_t matchScore,
            int64_t mismatchScore,
            int64_t gapScore,
            bool freeOnLeft,
            bool freeOnRight,
            vector< pair<bool, bool> >& alignment);

}



template<class Iterator>
    int64_t shasta::seqanAlign(
    Iterator begin0, Iterator end0,
    Iterator begin1, Iterator end1,
    int64_t matchScore,
    int64_t mismatchScore,
    int64_t gapScore,
    bool freeOnLeft,
    bool freeOnRight,
    vector< pair<bool, bool> >& alignment)
{
    // SeqAn does not handle empty sequences.
    SHASTA_ASSERT(begin0 != end0);
    SHASTA_ASSERT(begin1 != end1);

    // SeqAn types used below.
    using namespace seqan;
    using Int = typename Iterator::value_type;
    using Sequence = String<Int>;
    using StringSet = seqan::StringSet<Sequence>;
    using DepStringSet = seqan::StringSet<Sequence, Dependent<> >;
    using AlignGraph = Graph<seqan::Alignment<DepStringSet> >;

    // Fill in the sequences, adding 100 to all values
    // because SeqAn uses 45 to represent gaps.
    Sequence sequence0;
    for(Iterator it=begin0; it!=end0; ++it) {
        appendValue(sequence0, *it + 100);
    }
    Sequence sequence1;
    for(Iterator it=begin1; it!=end1; ++it) {
        appendValue(sequence1, *it + 100);
    }
    // Store them in a SeqAn string set.
    StringSet sequences;
    appendValue(sequences, sequence0);
    appendValue(sequences, sequence1);



    // Compute the alignment.
    // See https://docs.seqan.de/seqan/2.1.0/class_AlignConfig.html
    // for meaning of AlignConfig.
    AlignGraph graph(sequences);
    int64_t alignmentScore = 0;
    if(freeOnLeft) {
        if(freeOnRight) {
            // Free on both sides.
            alignmentScore = globalAlignment(
                graph,
                Score<int64_t, seqan::Simple>(matchScore, mismatchScore, gapScore),
                AlignConfig<true, true, true, true>(),
                LinearGaps());
        } else {
            // Free on left only.
            alignmentScore = globalAlignment(
                graph,
                Score<int64_t, seqan::Simple>(matchScore, mismatchScore, gapScore),
                AlignConfig<true, true, false, false>(),
                LinearGaps());
        }
    }else {
        if(freeOnRight) {
            // Free on right only.
            alignmentScore = globalAlignment(
                graph,
                Score<int64_t, seqan::Simple>(matchScore, mismatchScore, gapScore),
                AlignConfig<false, false, true, true>(),
                LinearGaps());
        } else {
            // Free on neither side.
            alignmentScore = globalAlignment(
                graph,
                Score<int64_t, seqan::Simple>(matchScore, mismatchScore, gapScore),
                AlignConfig<false, false, false, false>(),
                LinearGaps());
        }
    }



    // Extract the alignment from the graph.
    // This creates a single sequence consisting of the two rows
    // of the alignment, concatenated.
    Sequence align;
    convertAlignment(graph, align);
    const uint64_t totalAlignmentLength = seqan::length(align);
    SHASTA_ASSERT((totalAlignmentLength % 2) == 0);
    const uint64_t alignmentLength = totalAlignmentLength / 2;

    // Fill in the bool pairs representing the alignment.
    alignment.resize(alignmentLength);
    for(uint64_t i=0; i<alignmentLength; i++) {
        auto& p = alignment[i];
        p.first = not (align[i] == 45);
        p.second = not (align[i+alignmentLength] == 45);
    }


    return alignmentScore;
}

#endif

