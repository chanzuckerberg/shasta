// This file contains code for ReadGraph.creationMethod 2.

// Shasta.
#include "Assembler.hpp"
#include "orderPairs.hpp"
using namespace shasta;


void Assembler::createReadGraph2()
{
    // Parameters that control this function.
    // expose when code stabilizes.
    const uint64_t minSameStrandCoverage = 2;
    const uint64_t minOppositeStrandCoverage = 2;
    const uint64_t sampling = 100;
    const uint64_t maxAlignmentCount = 4;

    // Mark all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentData.size(), false);



    // Loop over reads. For each read, flag the alignments we want to keep.
    for(ReadId readId0=0; readId0<readCount(); readId0++) {

        // Work with this read on strand 0.
        const OrientedReadId orientedReadId0(readId0, 0);
        const uint32_t markerCount0 = uint32_t(markers.size(orientedReadId0.getValue()));

        // Get the alignments involving this oriented read.
        // This returns a vector of alignments with swaps and/or
        // reverse complementing already done, as necessary.
        vector<StoredAlignmentInformation> alignments;
        getStoredAlignments(orientedReadId0, alignments);

        // For each marker ordinal in orientedReadId0:
        // - Same strand alignment coverage is the number of alignments
        //   that marker is involved in, with other oriented reads
        //   on the same strand as orientedReadId0.
        // - Opposite strand alignment coverage is the number of alignments
        //   that marker is involved in, with other oriented reads
        //   on the opposite strand as orientedReadId0.
        // - Range coverage is the number of alignments such that
        //   that ordinal is internal to the alignment range.
        vector<uint64_t> sameStrandCoverage(markerCount0, 0);
        vector<uint64_t> oppositeStrandCoverage(markerCount0, 0);
        vector<uint64_t> rangeCoverage(markerCount0, 0);
        for(const StoredAlignmentInformation& s: alignments) {
            const OrientedReadId orientedReadId1 = s.orientedReadId;
            for(const auto& ordinals: s.alignment.ordinals) {
                const uint32_t ordinal0 = ordinals[0];
                SHASTA_ASSERT(ordinal0 < markerCount0);
                if(orientedReadId1.getStrand() == orientedReadId0.getStrand()) {
                    ++sameStrandCoverage[ordinal0];
                } else {
                    ++oppositeStrandCoverage[ordinal0];
                }
            }
            for(uint32_t ordinal0 = s.alignment.ordinals.front()[0];
                ordinal0 != s.alignment.ordinals.back()[0]; ordinal0++) {
                ++rangeCoverage[ordinal0];
            }
        }



        // Gather sampling marker ordinals such that:
        // - Same strand coverage >= minSameStrandCoverage.
        // - Opposite strand coverage >= minOppositeStrandCoverage.
        // The ordinals are stored in pairs with
        // coverage ratio (ratio of total coverage over range coverage).
        vector<pair<uint32_t, double> > samplingOrdinals;
        for(uint32_t ordinal0=0; ordinal0<markerCount0; ordinal0++) {
            if(sameStrandCoverage[ordinal0] < minSameStrandCoverage) {
                continue;
            }
            if(oppositeStrandCoverage[ordinal0] < minOppositeStrandCoverage) {
                continue;
            }
            const double coverageRatio =
                double(sameStrandCoverage[ordinal0] + oppositeStrandCoverage[ordinal0]) /
                double(rangeCoverage[ordinal0]);
            samplingOrdinals.push_back(make_pair(ordinal0, coverageRatio));
        }

        // Sort them in order of increasing coverage ratio.
        sort(samplingOrdinals.begin(), samplingOrdinals.end(), OrderPairsBySecondOnly<uint32_t, double>());

        // Keep only one every sampling markers.
        const uint64_t keepCount = markerCount0 / sampling;
        if(samplingOrdinals.size() > keepCount) {
            samplingOrdinals.resize(keepCount);
        }
        cout << readId0 << ": using " << samplingOrdinals.size() << " sampling ordinals." << endl;



        // Create an ordinal table which contains, for each ordinal
        // of orientedReadId0, aligned ordinals for each of the aligned
        // oriented reads.
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


        // For each aligned read, count the number of times each
        // of the sampling ordinals is:
        // - Part of the alignment.
        // - Inside the alignment range.
        vector< pair<uint64_t, uint64_t> > countTable;  // pair(i, alignedCount).
        for(uint64_t i=0; i<alignments.size(); i++) {
            const Alignment& alignment = alignments[i].alignment;
            const uint32_t firstAligned0 = alignment.ordinals.front()[0];
            const uint32_t lastAligned0 = alignment.ordinals.back()[0];
            uint64_t alignedCount = 0;
            uint64_t inRangeCount = 0;
            for(const auto& p: samplingOrdinals) {
                const uint32_t ordinal0 = p.first;
                if(ordinal0 < firstAligned0) {
                    continue;
                }
                if(ordinal0 > lastAligned0) {
                    continue;
                }
                ++inRangeCount;
                if(ordinalTable[ordinal0][i] != invalidOrdinal) {
                    ++alignedCount;
                }
            }
            countTable.push_back(make_pair(i, alignedCount));
        }

        // Only keep up to maxAlignmentCount.
        sort(countTable.begin(), countTable.end(),
            OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());
        if(countTable.size() > maxAlignmentCount) {
            countTable.resize(maxAlignmentCount);
        }

        cout << "*** " << readId0 << endl;
        for(const auto& p: countTable) {
            cout << p.first << " " << p.second << " " <<
                alignments[p.first].orientedReadId << " " << alignments[p.first].alignmentId << endl;
        }

        // Mark the alignments with the best aligned counts.
        for(const auto& p: countTable) {
            keepAlignment[alignments[p.first].alignmentId] = true;
        }
    }



    // Create the read graph using the alignments we selected.
    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;
    createReadGraphUsingSelectedAlignments(keepAlignment);

}
