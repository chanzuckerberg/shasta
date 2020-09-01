#include "Assembler.hpp"
using namespace shasta;



void Assembler::createReadGraph2(
    uint32_t maxAlignmentCount)
{
    // Percentile cutoffs for automated parameters
    double largePercentile = 0.12;
    double smallPercentile = 0.015;

    vector<bool> keepAlignment(alignmentData.size(), false);

    // Initialize histograms for measuring alignedFraction, markerCount, maxDrift, and maxSkip distributions
    Histogram2 alignedFractionHistogram(0, 1, 100);
    Histogram2 markerCountHistogram(0, 3000, 300);
    Histogram2 maxDriftHistogram(0, 100, 100);
    Histogram2 maxSkipHistogram(0, 100, 100);
    Histogram2 maxTrimHistogram(0, 100, 100);

    ofstream alignmentInfoCsv("alignmentInfo.csv");
    alignmentInfoCsv << "readId0" << ','
                     << "readId1" << ','
                     << "minAlignedFraction" << ','
                     << "markerCount" << ','
                     << "maxDrift" << ','
                     << "maxSkip" << ','
                     << "trim" << '\n';

    // Sample read pairs by discarding. Not ideal for large number of pairs AND small fractions
    for (size_t i=0; i<alignmentData.size(); i++){
        const auto info = alignmentData[i].info;
        const auto trims = info.computeTrim();
        const auto trim = max(trims.first, trims.second);

        alignedFractionHistogram.update(info.minAlignedFraction());
        markerCountHistogram.update(info.markerCount);
        maxDriftHistogram.update(info.maxDrift);
        maxSkipHistogram.update(info.maxSkip);
        maxTrimHistogram.update(trim);

        alignmentInfoCsv << alignmentData[i].readIds[0] << ','
                << alignmentData[i].readIds[1] << ','
                << info.minAlignedFraction() << ','
                << info.markerCount << ','
                << info.maxDrift << ','
                << info.maxSkip << ','
                << trim << '\n';
    }

    // Minimums
    double alignedFractionThreshold = alignedFractionHistogram.thresholdByCumulativeProportion(largePercentile);

    // MarkerCount is not at all gaussian, so it needs a different percentile.
    // This may be difficult to automate for varying length read sets.
    double markerCountThreshold = markerCountHistogram.thresholdByCumulativeProportion(smallPercentile);

    // Maximums use (1 - percentile)
    double maxDriftThreshold = maxDriftHistogram.thresholdByCumulativeProportion(1 - largePercentile);
    double maxSkipThreshold = maxSkipHistogram.thresholdByCumulativeProportion(1 - largePercentile);

    // MaxTrim also uses a small percentile because the default manual configuration was very permissive
    double maxTrimThreshold = maxTrimHistogram.thresholdByCumulativeProportion(1 - smallPercentile);

    cout << "Selected thresholds automatically for the following parameters:\n\t"
         << "alignedFraction:\t" << alignedFractionThreshold << "\n\t"
         << "markerCount:\t\t" << markerCountThreshold << "\n\t"
         << "maxDrift:\t\t" << maxDriftThreshold << "\n\t"
         << "maxSkip:\t\t" << maxSkipThreshold << "\n\t"
         << "maxTrim:\t\t" << maxTrimThreshold << "\n";

    // Find the number of reads and oriented reads.
    const ReadId orientedReadCount = uint32_t(markers.size());
    SHASTA_ASSERT((orientedReadCount % 2) == 0);
    const ReadId readCount = orientedReadCount / 2;

    // Vector to keep the alignments for each read,
    // with their number of markers.
    // Contains pairs(marker count, alignment id).
    vector< pair<uint32_t, uint32_t> > readAlignments;

    // Loop over reads.
    for(ReadId readId=0; readId<readCount; readId++) {

        // Gather the alignments for this read, each with its number of markers.
        readAlignments.clear();
        for(const uint32_t alignmentId: alignmentTable[OrientedReadId(readId, 0).getValue()]) {
            const AlignmentInfo& info = alignmentData[alignmentId].info;
            const auto trims = info.computeTrim();
            const auto trim = max(trims.first, trims.second);

            // If this alignment didnt pass the thresholding step, skip it
            if (info.minAlignedFraction() < alignedFractionThreshold){
                continue;
            }
            if (info.markerCount < markerCountThreshold){
                continue;
            }
            if (info.maxDrift > maxDriftThreshold){
                continue;
            }
            if (info.maxSkip > maxSkipThreshold){
                continue;
            }
            if (trim > maxTrimThreshold){
                continue;
            }

            // Otherwise add it to the list of candidate alignments for this read (to be ranked)
            readAlignments.push_back(make_pair(info.markerCount, alignmentId));
        }

        // Keep the best maxAlignmentCount.
        if(readAlignments.size() > maxAlignmentCount) {
            std::nth_element(
                    readAlignments.begin(),
                    readAlignments.begin() + maxAlignmentCount,
                    readAlignments.end(),
                    std::greater< pair<uint32_t, uint32_t> >());
            readAlignments.resize(maxAlignmentCount);
        }

        // Mark the surviving alignments as to be kept.
        for(const auto& p: readAlignments) {
            const uint32_t alignmentId = p.second;
            keepAlignment[alignmentId] = true;
        }
    }
    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

    ofstream alignedFractionHistogramCsv("AlignedFractionHistogram.csv");
    ofstream markerCountHistogramCsv("AlignmentMarkerCountHistogram.csv");
    ofstream maxDriftHistogramCsv("AlignmentDriftHistogram.csv");
    ofstream maxSkipHistogramCsv("AlignmentSkipHistogram.csv");
    ofstream maxTrimHistogramCsv("AlignmentTrimHistogram.csv");

    alignedFractionHistogram.writeToCsv(alignedFractionHistogramCsv, 3);
    markerCountHistogram.writeToCsv(markerCountHistogramCsv, 0);
    maxDriftHistogram.writeToCsv(maxDriftHistogramCsv, 0);
    maxSkipHistogram.writeToCsv(maxSkipHistogramCsv, 0);
    maxTrimHistogram.writeToCsv(maxTrimHistogramCsv, 0);

    createReadGraphUsingSelectedAlignments(keepAlignment);
}
