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
    vector<StoredAlignmentInformation> alignments;
    getStoredAlignments(orientedReadId0, alignments);
    cout << "Found " << alignments.size() << " alignments." << endl;

    // Check that all alignments are strictly increasing.
    for(const auto& p: alignments) {
        p.alignment.checkStrictlyIncreasing();
    }



    // Create an ordinal table which contains, for each ordinal
    // of orientedReadId0, aligned ordinals for each of the aligned
    // oriented reads.
    const uint32_t markerCount0 = uint32_t(markers.size(orientedReadId0.getValue()));
    const uint32_t invalidOrdinal = std::numeric_limits<uint32_t>::max();
    vector< vector<uint32_t> > ordinalTable(
        markerCount0, vector<uint32_t>(alignments.size(), invalidOrdinal));
    for(uint64_t i=0; i<alignments.size(); i++) {
        const Alignment& alignment = alignments[i].alignment;
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
        const OrientedReadId orientedReadId1 = alignments[i].orientedReadId;
        const Alignment& alignment = alignments[i].alignment;
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
        "Range coverage,Same strand range coverage,Opposite strand range coverage,"
        "Coverage ratio,Same strand coverage ratio,Opposite strand coverage ratio,";
    for(const auto& p: alignments) {
        csv << p.orientedReadId << ",";
    }
    csv << "\n";



    // Write the ordinal table to the csv file.
    for(uint32_t ordinal0=0; ordinal0<markerCount0; ordinal0++) {
        const uint64_t cSameStrand = coverage[ordinal0][0];
        const uint64_t cOppositeStrand = coverage[ordinal0][1];
        const uint64_t c = cSameStrand + cOppositeStrand;
        const uint64_t rcSameStrand = rangeCoverage[ordinal0][0];
        const uint64_t rcOppositeStrand = rangeCoverage[ordinal0][1];
        const uint64_t rc = rcSameStrand + rcOppositeStrand;
        const double rSameStrand  = double(cSameStrand) / double(rcSameStrand);
        const double rOppositeStrand  = double(cOppositeStrand) / double(rcOppositeStrand);
        const double r = double(c) / double(rc);

        csv << ordinal0 << ",";
        csv << c << ",";
        csv << cSameStrand << ",";
        csv << cOppositeStrand << ",";
        csv << rc << ",";
        csv << rcSameStrand << ",";
        csv << rcOppositeStrand << ",";
        csv << r << ",";
        csv << rSameStrand << ",";
        csv << rOppositeStrand << ",";
        for(uint64_t i=0; i<alignments.size(); i++) {
            const uint32_t ordinal1 = ordinalTable[ordinal0][i];
            if(ordinal1 != invalidOrdinal) {
                csv << ordinal1;
            } else {
                const Alignment& alignment = alignments[i].alignment;
                const uint32_t alignmentBegin0 = alignment.ordinals.front()[0];
                const uint32_t alignmentEnd0 = alignment.ordinals.back()[0];
                if((ordinal0 >= alignmentBegin0) and (ordinal0 <= alignmentEnd0)) {
                    csv << "No";
                }
            }
            csv << ",";
        }
        csv << "\n";
    }



    // Compute coverage histograms and write them out.
    // 0 = coverage
    // 1 = same strand coverage
    // 2 = opposite strand coverage.
    // 3 = range coverage
    // 4 = same strand range coverage
    // 5 = opposite strand range coverage.
    // Ratio histogram:
    // 0 = coverage ratio (binned).
    // 1 = same strand coverage ratio (binned).
    // 2 = opposite strand coverage ratio (binned).
    vector< array<uint64_t, 6> > histogram;
    const uint64_t binCount = 10;
    const double binSize = 1. / double(binCount);
    vector< array<uint64_t, 3> > ratioHistogram(binCount + 1, {0,0,0});
    for(uint32_t ordinal0=0; ordinal0<markerCount0; ordinal0++) {
        const uint64_t cSameStrand = coverage[ordinal0][0];
        const uint64_t cOppositeStrand = coverage[ordinal0][1];
        const uint64_t c = cSameStrand + cOppositeStrand;
        const uint64_t rcSameStrand = rangeCoverage[ordinal0][0];
        const uint64_t rcOppositeStrand = rangeCoverage[ordinal0][1];
        const uint64_t rc = rcSameStrand + rcOppositeStrand;
        const double rSameStrand  = (rcSameStrand==0 ? 0. : double(cSameStrand) / double(rcSameStrand));
        const double rOppositeStrand  = (rcOppositeStrand==0 ? 0. : double(cOppositeStrand) / double(rcOppositeStrand));
        const double r = (rc==0 ? 0. : double(c) / double(rc));
        const uint64_t irSameStrand = uint64_t(rSameStrand/binSize);
        const uint64_t irOppositeStrand = uint64_t(rOppositeStrand/binSize);
        const uint64_t ir = uint64_t(r/binSize);

        SHASTA_ASSERT(cSameStrand <= rcSameStrand);
        SHASTA_ASSERT(cOppositeStrand <= rcOppositeStrand);

        if(histogram.size() <= rc) {
            histogram.resize(rc + 1, {0,0,0,0,0,0,});
        }
        ++histogram[c][0];
        ++histogram[cSameStrand][1];
        ++histogram[cOppositeStrand][2];
        ++histogram[rc][3];
        ++histogram[rcSameStrand][4];
        ++histogram[rcOppositeStrand][5];
        ++ratioHistogram[ir][0];
        ++ratioHistogram[irSameStrand][1];
        ++ratioHistogram[irOppositeStrand][2];
    }
    ofstream csv2("AlignmentCoverageHistogram.csv");
    csv2 << "Coverage value,Total,Same strand,Opposite strand,"
        "Range total, Range same strand, Range opposite strand\n";
    for(uint64_t c=0; c<histogram.size(); c++) {
        csv2 << c << ",";
        for(uint64_t i=0; i<6; i++) {
            csv2 << histogram[c][i] << ",";
        }
        csv2 << "\n";
    }
    ofstream csv3("AlignmentCoverageRatioHistogram.csv");
    csv3 << "Coverage ratio,Total,Same strand,Opposite strand\n";
    for(uint64_t c=0; c<ratioHistogram.size(); c++) {
        csv3 << double(c)*binSize << ",";
        for(uint64_t i=0; i<3; i++) {
            csv3 << ratioHistogram[c][i] << ",";
        }
        csv3 << "\n";
    }
}



// Get the stored compressed alignments involving a given oriented read.
// This performs swaps and reverse complementing as necessary,
// To return alignments in which the first oriented read is
// the one specified as the argument.
void Assembler::getStoredAlignments(
    OrientedReadId orientedReadId0,
    vector<StoredAlignmentInformation> & alignments) const
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
        Alignment& alignment = alignments.back().alignment;
        OrientedReadId& orientedReadId1 = alignments.back().orientedReadId;
        alignments.back().alignmentId = alignmentIndex;
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
