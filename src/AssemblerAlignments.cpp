#include "Assembler.hpp"
#include "compressAlignment.hpp"
using namespace shasta;



// Analyze the stored alignments involving a given oriented read.
void Assembler::analyzeAlignments(ReadId readId0, Strand strand0) const
{
    const OrientedReadId orientedReadId0(readId0, strand0);
    cout << "Analyzing stored alignments for " << orientedReadId0 << endl;

    // Get the alignments involving this oriented read.
    // This returns a vector alignments with swaps and/or
    // reverse complementing already done, as necessary.
    vector< pair<OrientedReadId, Alignment> > alignments;
    getStoredAlignments(orientedReadId0, alignments);
    cout << "Found " << alignments.size() << " alignments." << endl;

    // Check that all alignments are strictly increasing.
    for(const auto& p: alignments) {
        p.second.checkStrictlyIncreasing();
    }



    // Create an ordinal table which contains, for each ordinal
    // of orientedReadId0, aligned ordinals for each of the aligned
    // oriented reads.
    const uint32_t markerCount0 = uint32_t(markers.size(orientedReadId0.getValue()));
    const uint32_t invalidOrdinal = std::numeric_limits<uint32_t>::max();
    vector< vector<uint32_t> > ordinalTable(
        markerCount0, vector<uint32_t>(alignments.size(), invalidOrdinal));
    for(uint64_t i=0; i<alignments.size(); i++) {
        const Alignment& alignment = alignments[i].second;
        for(const auto& o: alignment.ordinals) {
            const uint32_t ordinal0 = o[0];
            const uint32_t ordinal1 = o[1];
            SHASTA_ASSERT(ordinal0 < markerCount0);
            ordinalTable[ordinal0][i] = ordinal1;
        }
    }



    // Compute coverage for each marker and for each strand
    // (0=same strand, 1 = opposite strands).
    // Range coverage is the number of alignments whose range covers each ordinal.
    vector< array<uint32_t, 2> > coverage(markerCount0, {0, 0});
    vector< array<uint32_t, 2> > rangeCoverage(markerCount0, {0, 0});
    for(uint64_t i=0; i<alignments.size(); i++) {
        const OrientedReadId orientedReadId1 = alignments[i].first;
        const Alignment& alignment = alignments[i].second;
        const auto strandIndex = (orientedReadId0.getStrand() == orientedReadId1.getStrand()) ? 0 : 1;

        // Update coverage for this alignment.
        for(const auto& o: alignment.ordinals) {
            const uint32_t ordinal0 = o[0];
            SHASTA_ASSERT(ordinal0 < markerCount0);
            ++coverage[ordinal0][strandIndex];
        }

        // Update range coverage for this alignment.
        for(uint32_t ordinal0=alignment.ordinals.front()[0];
            ordinal0<=alignment.ordinals.back()[0]; ordinal0++) {
            ++rangeCoverage[ordinal0][strandIndex];
        }
    }



    // Create the csv file and write the header.
    ofstream csv("Alignments.csv");
    csv << "Ordinal0,Coverage,Same strand coverage,Opposite strand coverage,"
        "Range coverage,Same strand range coverage,Opposite strand range coverage,";
    for(const auto& p: alignments) {
        csv << p.first << ",";
    }
    csv << "\n";



    // Write the ordinal table to the csv file.
    for(uint32_t ordinal0=0; ordinal0<markerCount0; ordinal0++) {
        csv << ordinal0 << ",";
        csv << coverage[ordinal0][0]+coverage[ordinal0][1] << ",";
        csv << coverage[ordinal0][0] << ",";
        csv << coverage[ordinal0][1] << ",";
        csv << rangeCoverage[ordinal0][0]+rangeCoverage[ordinal0][1] << ",";
        csv << rangeCoverage[ordinal0][0] << ",";
        csv << rangeCoverage[ordinal0][1] << ",";
        for(const uint32_t ordinal1: ordinalTable[ordinal0]) {
            if(ordinal1 != invalidOrdinal) {
                csv << ordinal1;
            }
            csv << ",";
        }
        csv << "\n";
    }
}



// Get the stored compressed alignments involving a given oriented read.
// This performs swaps and reverse complementing as necessary,
// To return alignments in which the first oriented read is
// the one specified as the argument.
void Assembler::getStoredAlignments(
    OrientedReadId orientedReadId0,
    vector< pair<OrientedReadId, Alignment> > & alignments) const
{
    // Check that we have what we need.
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    SHASTA_ASSERT(compressedAlignments.isOpen());

    // Access the alignment table portion for this oriented read.
    // It contains indexes into alignmentData and compressedAlignments
    // for alignments involving this oriented read.
    const span<const uint32_t> alignmentIndexes = alignmentTable[orientedReadId0.getValue()];



    // Loop over alignments involving this oriented read.
    alignments.clear();
    for(const uint32_t alignmentIndex: alignmentIndexes) {

        // Access the stored information we have about this alignment.
        AlignmentData alignmentData = this->alignmentData[alignmentIndex];
        const span<const char> compressedAlignment = compressedAlignments[alignmentIndex];

        // The alignment is stored with its first read on strand 0.
        OrientedReadId alignmentOrientedReadId0(alignmentData.readIds[0], 0);
        OrientedReadId alignmentOrientedReadId1(alignmentData.readIds[1],
            alignmentData.isSameStrand ? 0 : 1);

        // Decompress the alignment.
        alignments.resize(alignments.size() + 1);
        Alignment& alignment = alignments.back().second;
        OrientedReadId& orientedReadId1 = alignments.back().first;
        decompress(compressedAlignment, alignment);
        SHASTA_ASSERT(alignment.ordinals.size() == alignmentData.info.markerCount);



        // Tweak the alignment to make sure its first oriented read is orientedReadId0.
        // This may require a swap and/or reverse complement.

        // Do a swap, if needed.
        if(alignmentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
            alignment.swap();
            swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
        }
        SHASTA_ASSERT(alignmentOrientedReadId0.getReadId() == orientedReadId0.getReadId());

        // Reverse complement, if needed.
        if(alignmentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
            alignment.reverseComplement(
                uint32_t(markers.size(alignmentOrientedReadId0.getValue())),
                uint32_t(markers.size(alignmentOrientedReadId1.getValue())));
            alignmentOrientedReadId0.flipStrand();
            alignmentOrientedReadId1.flipStrand();
        }
        SHASTA_ASSERT(alignmentOrientedReadId0 == orientedReadId0);
        orientedReadId1 = alignmentOrientedReadId1;
    }
}
