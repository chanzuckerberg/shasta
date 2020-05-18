// PngImage.hpp must be included first because of png issues on Ubuntu 16.04.
#include "PngImage.hpp"

#include "Assembler.hpp"
using namespace shasta;


// Seqan.
#include <seqan/align.h>

#include <numeric>


// Align two oriented reads using SeqAn banded alignment.
// This id done in two steps:
// 1. Compute an alignment (unbanded) using downsampled marker
//    sequences for the two oriented reads.
// 2. Use the downsampled alignment to compute a band.
//    Then do a banded alignment using that band.

void Assembler::alignOrientedReads3(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int matchScore,
    int mismatchScore,
    int gapScore,
    double downsamplingFactor,  // The fraction of markers to keep in the first step.
    int bandExtend,             // How much to extend the band computed in the first step.
    Alignment& alignment,
    AlignmentInfo& alignmentInfo)
{
    const bool debug = false;
    if(debug) {
        cout << "Assembler::alignOrientedReads3 begins for " <<
            orientedReadId0 << " " << orientedReadId1 << endl;
    }

    // Hide shasta::Alignment.
    using namespace seqan;
    using seqan::Alignment;
    const uint32_t seqanGapValue = 45;

    // An oriented read is represented as a sequence of KmerId
    // (the KmerId's of its markers). We want to align a pair of
    // such sequences.
    using TSequence = String<KmerId>;

    // Other SeqAn types we need.
    using TStringSet = StringSet<TSequence>;
    using TDepStringSet = StringSet<TSequence, Dependent<> >;
    using TAlignGraph = Graph<Alignment<TDepStringSet> >;


    // Get the markers for the two oriented reads.
    array<span<CompressedMarker>, 2> allMarkers;
    allMarkers[0] = markers[orientedReadId0.getValue()];
    allMarkers[1] = markers[orientedReadId1.getValue()];

    // Vectors to contain downsampled markers.
    // For each of the two reads we store vectors of
    // (ordinal, KmerId).
    array< vector<pair<uint32_t, KmerId> >, 2> downsampledMarkers;
    array<TSequence, 2> downsampledSequences;

    // Fill in downsampled markers.
    // SeqAn uses 45 to represent gaps, so we add 45 to the KmerIds passed to SeqAn.
    // This means that we can't do k=16.
    const uint32_t hashThreshold =
        uint32_t(downsamplingFactor * double(std::numeric_limits<uint32_t>::max()));
    for(uint64_t i=0; i<2; i++) {
        for(uint32_t ordinal=0; ordinal<uint32_t(allMarkers[i].size()); ordinal++) {
            const KmerId kmerId = allMarkers[i][ordinal].kmerId;
             if(kmerTable[kmerId].hash < hashThreshold) {
                downsampledMarkers[i].push_back(make_pair(ordinal, kmerId));
                appendValue(downsampledSequences[i], kmerId + 100);
            }
        }
    }

    if(debug) {
        cout << "Aligning two oriented reads with " <<
            allMarkers[0].size() << " and " << allMarkers[1].size() << " markers." << endl;
        cout << "Downsampled markers for step 1 to " <<
            downsampledMarkers[0].size() << " and " <<
            downsampledMarkers[1].size() << " markers." << endl;


        for(uint64_t i=0; i<2; i++) {
            ofstream csv("OrientedReadDownsampled-" + to_string(i) + ".csv");
            for(const auto& p: downsampledMarkers[i]) {
                csv << p.first << "," << p.second << "\n";
            }
        }

    }



    // Use SeqAn to compute an alignment of the downsampled markers, free at both ends.
    // https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html

    // Store them in a SeqAn string set.
    TStringSet downsampledSequencesSet;
    appendValue(downsampledSequencesSet, downsampledSequences[0]);
    appendValue(downsampledSequencesSet, downsampledSequences[1]);

    // Compute the alignment.
    TAlignGraph downsampledGraph(downsampledSequencesSet);
    const int downsampledScore = globalAlignment(
        downsampledGraph,
        Score<int, Simple>(matchScore, mismatchScore, gapScore),
        AlignConfig<true, true, true, true>(),
        LinearGaps());
    if(debug) {
        cout << "Downsampled alignment score is " << downsampledScore << endl;
    }

    // Extract the alignment from the graph.
    // This creates a single sequence consisting of the two rows
    // of the alignment, concatenated.
    TSequence downsampledAlign;
    convertAlignment(downsampledGraph, downsampledAlign);
    const int totalDownsampledAlignmentLength = int(seqan::length(downsampledAlign));
    SHASTA_ASSERT((totalDownsampledAlignmentLength % 2) == 0);    // Because we are aligning two sequences.
    const int downsampledAlignmentLength = totalDownsampledAlignmentLength / 2;
    if(debug) {
        cout << "Downsampled alignment length " << downsampledAlignmentLength << endl;
    }



    // Write the downsampled alignment on its alignment matrix.
    if(debug) {
        PngImage image(int(downsampledMarkers[0].size()), int(downsampledMarkers[1].size()));
        for(uint64_t i0=0; i0<downsampledMarkers[0].size(); i0++) {
            const KmerId kmerId0 = downsampledMarkers[0][i0].second;
            for(uint64_t i1=0; i1<downsampledMarkers[1].size(); i1++) {
                const KmerId kmerId1 = downsampledMarkers[1][i1].second;
                if(kmerId1 == kmerId0) {
                    image.setPixel(int(i0), int(i1), 255, 0, 0);
                }
            }
        }

        uint32_t i0 = 0;
        uint32_t i1 = 0;
        for(int i=0;
            i<downsampledAlignmentLength and
            i0<downsampledMarkers[0].size() and
            i1<downsampledMarkers[1].size(); i++) {
            if( downsampledAlign[i] != seqanGapValue and
                downsampledAlign[i + downsampledAlignmentLength] != seqanGapValue) {
                if(downsampledMarkers[0][i0].second == downsampledMarkers[1][i1].second) {
                    image.setPixel(int(i0), int(i1), 0, 255, 0);
                } else {
                    image.setPixel(int(i0), int(i1), 80, 80, 0);
                }
            }
            if(downsampledAlign[i] != seqanGapValue) {
                ++i0;
            }
            if(downsampledAlign[i + downsampledAlignmentLength] != seqanGapValue) {
                ++i1;
            }
        }

        image.write("DownsampledAlignment.png");
    }



    // If the downsampled alignment is empty, just return an empty alignment.
    if(uint64_t(downsampledAlignmentLength) ==
        downsampledMarkers[0].size() + downsampledMarkers[1].size()) {
        alignment.clear();
        alignmentInfo.create(
            alignment, uint32_t(allMarkers[0].size()), uint32_t(allMarkers[1].size()));
        return;
    }



    // Use the downsampled alignment to compute the band to be used
    // for the full alignment.
    int32_t offsetMin = std::numeric_limits<int32_t>::max();
    int32_t offsetMax = std::numeric_limits<int32_t>::min();
    uint32_t i0 = 0;
    uint32_t i1 = 0;
    for(int i=0;
        i<downsampledAlignmentLength and
        i0<downsampledMarkers[0].size() and
        i1<downsampledMarkers[1].size(); i++) {
        if( downsampledAlign[i] != seqanGapValue and
            downsampledAlign[i + downsampledAlignmentLength] != seqanGapValue) {
            if(downsampledMarkers[0][i0].second == downsampledMarkers[1][i1].second) {
                const int32_t offset =
                    int32_t(downsampledMarkers[0][i0].first) -
                    int32_t(downsampledMarkers[1][i1].first);
                offsetMin = min(offsetMin, offset);
                offsetMax = max(offsetMax, offset);
            }
        }
        if(downsampledAlign[i] != seqanGapValue) {
            ++i0;
        }
        if(downsampledAlign[i + downsampledAlignmentLength] != seqanGapValue) {
            ++i1;
        }
    }
    const int32_t bandMin = offsetMin - bandExtend;
    const int32_t bandMax = offsetMax + bandExtend;
    if(debug) {
        cout << "Offset range " << offsetMin << " " << offsetMax << endl;
        cout << "Banded alignment will use band " << bandMin << " " << bandMax << endl;
    }



    // Now, do a alignment using this band and all markers.
    array<TSequence, 2> sequences;
    for(uint64_t i=0; i<2; i++) {
        for(uint32_t ordinal=0; ordinal<uint32_t(allMarkers[i].size()); ordinal++) {
            const KmerId kmerId = allMarkers[i][ordinal].kmerId;
            appendValue(sequences[i], kmerId + 100);
        }
    }
    TStringSet sequencesSet;
    appendValue(sequencesSet, sequences[0]);
    appendValue(sequencesSet, sequences[1]);
    TAlignGraph graph(sequencesSet);
    const int score = globalAlignment(
        graph,
        Score<int, Simple>(matchScore, mismatchScore, gapScore),
        AlignConfig<true, true, true, true>(),
        bandMin, bandMax,
        LinearGaps());
    if(debug) {
        cout << "Full alignment score is " << score << endl;
    }
    TSequence align;
    convertAlignment(graph, align);
    const int totalAlignmentLength = int(seqan::length(align));
    SHASTA_ASSERT((totalAlignmentLength % 2) == 0);    // Because we are aligning two sequences.
    const int alignmentLength = totalAlignmentLength / 2;
    if(debug) {
        cout << "Full alignment length " << alignmentLength << endl;
    }



    // Fill in the alignment.
    alignment.clear();
    uint32_t ordinal0 = 0;
    uint32_t ordinal1 = 0;
    for(int i=0;
        i<alignmentLength and ordinal0<allMarkers[0].size() and ordinal1<allMarkers[1].size(); i++) {
        if( align[i] != seqanGapValue and
            align[i + alignmentLength] != seqanGapValue and
            allMarkers[0][ordinal0].kmerId == allMarkers[1][ordinal1].kmerId) {
            alignment.ordinals.push_back(array<uint32_t, 2>{ordinal0, ordinal1});
        }
        if(align[i] != seqanGapValue) {
            ++ordinal0;
        }
        if(align[i + alignmentLength] != seqanGapValue) {
            ++ordinal1;
        }
    }

    // Check how close to the band we got.
    int32_t distanceToLowerDiagonal = std::numeric_limits<int32_t>::max();
    int32_t distanceToUpperDiagonal = std::numeric_limits<int32_t>::max();
    for(const auto& ordinals: alignment.ordinals) {
        const int32_t offset = int32_t(ordinals[0]) - int32_t(ordinals[1]);
        distanceToLowerDiagonal = min(distanceToLowerDiagonal, offset-bandMin);
        distanceToUpperDiagonal = min(distanceToUpperDiagonal, bandMax-offset);
    }
    if(debug) {
        cout << "Distance of alignment from band " <<
            distanceToLowerDiagonal << " " << distanceToUpperDiagonal << endl;
    }

    // Store the alignment info.
    alignmentInfo.create(alignment, uint32_t(allMarkers[0].size()), uint32_t(allMarkers[1].size()));

}

