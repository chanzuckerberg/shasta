#include "Assembler.hpp"
using namespace shasta;



void Assembler::createReadGraph2(
    uint32_t maxAlignmentCount,
    uint32_t maxTrim)
{
    vector<bool> keepAlignment(alignmentData.size(), false);

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

    // Evaluate thresholds by finding the 12th percentile
    double cumulativeProportion = 0.12;

    // Minimums
    double alignedFractionThreshold = alignedFractionHistogram.thresholdByCumulativeProportion(cumulativeProportion);

    // MarkerCount is not at all gaussian, so it needs a different percentile.
    // This may be difficult to automate for varying length read sets.
    double markerCountThreshold = markerCountHistogram.thresholdByCumulativeProportion(0.015);

    // Maximums use (1 - percentile)
    double maxDriftThreshold = maxDriftHistogram.thresholdByCumulativeProportion(1 - cumulativeProportion);
    double maxSkipThreshold = maxSkipHistogram.thresholdByCumulativeProportion(1 - cumulativeProportion);

    cout << "Selected thresholds automatically for the following parameters:\n\t"
         << "alignedFraction:\t" << alignedFractionThreshold << "\n\t"
         << "markerCount:\t\t" << markerCountThreshold << "\n\t"
         << "maxDrift:\t\t" << maxDriftThreshold << "\n\t"
         << "maxSkip:\t\t" << maxSkipThreshold << "\n";

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

    createReadGraphUsingSelectedAlignments(keepAlignment);
}
