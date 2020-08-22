#include "Assembler.hpp"
using namespace shasta;



void Assembler::createReadGraph2(
    uint32_t maxAlignmentCount,
    uint32_t maxTrim)
{
    // This boilerplate code creates the read graph using
    // all available alignments.
    vector<bool> keepAlignment(alignmentData.size(), true);

    double candidateSampleFraction = 0.3;
    if (candidateSampleFraction > 1.0){
        throw runtime_error("ERROR: sample fraction must be <= 1.0");
    }

    // Initialize histograms for measuring alignedFraction, markerCount, maxDrift, and maxSkip distributions
    Histogram2 alignedFractionHistogram(0, 1, 50);
    Histogram2 markerCountHistogram(0, 3000, 240);
    Histogram2 maxDriftHistogram(0, 100, 100);
    Histogram2 maxSkipHistogram(0, 100, 100);

    // Sample read pairs by discarding. Not ideal for large number of pairs AND small fractions
    for (size_t i=0; i<alignmentData.size(); i++){
        if (i % 101 <= size_t(round(candidateSampleFraction * 101))){
            const auto info = alignmentData[i].info;

            alignedFractionHistogram.update(info.minAlignedFraction());
            markerCountHistogram.update(info.markerCount);
            maxDriftHistogram.update(info.maxDrift);
            maxSkipHistogram.update(info.maxSkip);
        }
    }

    // Evaluate thresholds by finding the 5th percentile
    double cumulativeProportion = 0.05;

    // Minimums
    double alignedFractionThreshold = alignedFractionHistogram.thresholdByCumulativeProportion(cumulativeProportion);
    double markerCountThreshold = markerCountHistogram.thresholdByCumulativeProportion(cumulativeProportion);

    // Maximums use (1 - percentile)
    double maxDriftThreshold = maxDriftHistogram.thresholdByCumulativeProportion(1 - cumulativeProportion);
    double maxSkipThreshold = maxSkipHistogram.thresholdByCumulativeProportion(1 - cumulativeProportion);

    cout << "Selected thresholds automatically for the following parameters:\n\t"
         << "alignedFraction:\t" << alignedFractionThreshold << "\n\t"
         << "markerCount:\t\t" << markerCountThreshold << "\n\t"
         << "maxDrift:\t\t" << maxDriftThreshold << "\n\t"
         << "maxSkip:\t\t" << maxSkipThreshold << "\n";

    // Flag failing alignments
    for (size_t i=0; i<alignmentData.size(); i++) {
        const auto info = alignmentData[i].info;

        if (info.minAlignedFraction() < alignedFractionThreshold){
            keepAlignment[i] = false;
        }
        if (info.markerCount < markerCountThreshold){
            keepAlignment[i] = false;
        }
        if (info.maxDrift > maxDriftThreshold){
            keepAlignment[i] = false;
        }
        if (info.maxSkip > maxSkipThreshold){
            keepAlignment[i] = false;
        }
    }

    createReadGraphUsingSelectedAlignments(keepAlignment);
}
