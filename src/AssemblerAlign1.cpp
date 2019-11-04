// Alternative alignment functions with 1 suffix.
#include "Assembler.hpp"
using namespace shasta;

// Seqan.
#ifdef __linux__
#include <seqan/align.h>
#endif



#ifndef __linux__

// For macOS we don't have SeqAn, so we can't do any of this.
void Assembler::alignOrientedReads1(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1)
{
    throw runtime_error("alignOrientedReads1 is not available on macOS.");
}

#else



void Assembler::alignOrientedReads1(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    int matchScore,
    int mismatchScore,
    int gapScore)
{
    alignOrientedReads1(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1),
        matchScore, mismatchScore, gapScore);
}



void Assembler::alignOrientedReads1(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int matchScore,
    int mismatchScore,
    int gapScore)
{
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    alignOrientedReads1(
        orientedReadId0,
        orientedReadId1,
        matchScore,
        mismatchScore,
        gapScore,
        alignment,
        alignmentInfo);
}



void Assembler::alignOrientedReads1(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int matchScore,
    int mismatchScore,
    int gapScore,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo)
{

    // Use SeqAn to compute an alignment free at both ends.
    // https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html
    using namespace seqan;

    // Hide shasta::Alignment.
    using seqan::Alignment;

    // An oriented read is represented as a sequence of KmerId
    // (the KmerId's of its markers). We want to align a pair of
    // such sequences.
    using TSequence = String<KmerId>;

    // Other SeqAn types we need.
    using TStringSet = StringSet<TSequence>;
    using TDepStringSet = StringSet<TSequence, Dependent<> >;
    using TAlignGraph = Graph<Alignment<TDepStringSet> >;

    // Access the markers of our oriented reads.
    const MemoryAsContainer<CompressedMarker> markers0 =
        markers[orientedReadId0.getValue()];
    const MemoryAsContainer<CompressedMarker> markers1 =
        markers[orientedReadId1.getValue()];



    // Seqan uses the integer 45 to represent a gap
    // and I did not find a good way to control that.
    // So if KmerId 45 is a marker we replace it with the first KmerId
    // that does not represent a marker.
    // This is messy but I did not find a better solution.
    bool replacementIsNeeded = false;
    const KmerId seqanGapValue = 45;
    KmerId replacementValue = seqanGapValue;
    if(kmerTable[seqanGapValue].isMarker) {
        replacementIsNeeded = true;
        for(uint64_t i=0; i<kmerTable.size(); i++) {
            if(!kmerTable[i].isMarker) {
                replacementValue = KmerId(i);
                break;
            }
        }
        cout << "Replacement value " << replacementValue << endl;
        SHASTA_ASSERT(replacementValue != seqanGapValue);
    }



    // Construct the sequences of KmerId's we want to align.
    TSequence seq0;
    for(const CompressedMarker marker: markers0) {
        if(replacementIsNeeded && marker.kmerId == seqanGapValue) {
            appendValue(seq0, replacementValue);
        } else {
            appendValue(seq0, marker.kmerId);
        }
    }
    TSequence seq1;
    for(const CompressedMarker marker: markers1) {
        if(replacementIsNeeded && marker.kmerId == seqanGapValue) {
            appendValue(seq1, replacementValue);
        } else {
            appendValue(seq1, marker.kmerId);
        }
    }

    // Store them in a SeqAn string set.
    TStringSet sequences;
    appendValue(sequences, seq0);
    appendValue(sequences, seq1);

    // Compute the alignment.
    TAlignGraph graph(sequences);
    const int score = globalAlignment(
        graph,
        Score<int, Simple>(matchScore, mismatchScore, gapScore),
        AlignConfig<true, true, true, true>(),
        LinearGaps());
    cout << "Number of markers in these oriented reads: " <<
        markers0.size() << " " << markers1.size() << endl;
    cout << "Alignment score is " << score << endl;

    // Extract the alignment from the graph.
    // This creates a single sequence consisting of the two rows
    // of the alignment, concatenated.
    TSequence align;
    convertAlignment(graph, align);
    const int totalAlignmentLength = int(seqan::length(align));
    SHASTA_ASSERT((totalAlignmentLength % 2) == 0);    // Because we are aligning two sequences.
    const int alignmentLength = totalAlignmentLength / 2;
    cout << "Alignment length " << alignmentLength << endl;



    // Fill in the alignment.
    alignment.ordinals.clear();
    uint32_t ordinal0 = 0;
    uint32_t ordinal1 = 0;
    for(int i=0; i<alignmentLength; i++) {
        cout << i << " " << ordinal0 << " " << ordinal1 << endl;
        if( markers0[ordinal0].kmerId != seqanGapValue and
            markers0[ordinal0].kmerId == markers1[ordinal1].kmerId) {
            alignment.ordinals.push_back(array<uint32_t, 2>{ordinal0, ordinal1});
        }
        cout << "***A" << endl;
        if(align[i] != seqanGapValue) {
            ++ordinal0;
        }
        cout << "***B" << endl;
        if(align[i + alignmentLength] != seqanGapValue) {
            ++ordinal1;
        }
        cout << "***C" << endl;
    }
    SHASTA_ASSERT(ordinal0 == markers0.size());
    SHASTA_ASSERT(ordinal1 == markers1.size());

    // Store the alignment info.
    alignmentInfo.create(alignment, uint32_t(markers0.size()), uint32_t(markers1.size()));

#if 0
    // Extract the two rows of the alignment.
    array<vector<uint32_t>, 2> alignment;
    alignment[0].resize(alignmentLength);
    alignment[1].resize(alignmentLength);
    for(int i=0; i<alignmentLength; i++) {
        alignment[0][i] = align[i];
        alignment[1][i] = align[i + alignmentLength];
    }



    // Write out the alignment.
    for(int i=0; i<alignmentLength; i++) {
        cout << i << " ";
        if(alignment[0][i] == seqanGapValue) {
            cout << "-";
        } else {
            cout << alignment[0][i];
        }
        cout << " ";
        if(alignment[1][i] == seqanGapValue) {
            cout << "-";
        } else {
            cout << alignment[1][i];
        }
        if(
            alignment[0][i]!=seqanGapValue and
            alignment[1][i]!=seqanGapValue and
            alignment[0][i]==alignment[1][i]) {
            cout << "***";
        }
        cout << "\n";
    }
    cout << flush;
#endif

}

#endif


