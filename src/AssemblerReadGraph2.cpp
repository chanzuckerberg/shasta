#include "Assembler.hpp"
#include "Histogram.hpp"
#include "Reads.hpp"
using namespace shasta;


/*******************************************************************************

	The original code had the following definitions and comments.
	These defaults are now set in AssemblerOptions.cpp and
	can be modified by the user using command line options
	in the [ReadGraph] section.

    double alignedFractionPercentile = 0.12;
    double maxDriftPercentile = 0.12;
    double maxSkipPercentile = 0.12;

    // MarkerCount is not at all gaussian, so it needs a different percentile.
    // This may be difficult to automate for varying length read sets.
    double markerCountPercentile = 0.015;

    // MaxTrim also uses a small percentile because the default manual configuration was very permissive
    double maxTrimPercentile = 0.015;

*******************************************************************************/


void Assembler::writeReadGraphEdges(bool useReadName) const{
    string path = "ReadGraphEdges.csv";
    path = filesystem::getAbsolutePath(path);
    ofstream file(path);

    if(not file.good()){
        throw runtime_error("ERROR: file could not be written: " + path);
    }

    if(not useReadName) {
        file << "ReadId0,ReadId1,SameStrand\n";

        // Loop over readgraph edges
        for (auto edge = readGraph.edges.begin(); edge != readGraph.edges.end(); std::advance(edge,2)) {
            bool isSameStrand = (edge->orientedReadIds[0].getStrand() == edge->orientedReadIds[1].getStrand());

            file << edge->orientedReadIds[0].getReadId() << ','
                 << edge->orientedReadIds[1].getReadId() << ','
                 << (isSameStrand ? "Yes" : "No") << '\n';
        }
    }
    else{
        file << "ReadName0,ReadName1,SameStrand\n";

        // Loop over readgraph edges
        for (auto edge = readGraph.edges.begin(); edge != readGraph.edges.end(); std::advance(edge,2)) {
            bool isSameStrand = (edge->orientedReadIds[0].getStrand() == edge->orientedReadIds[1].getStrand());

            file << reads->getReadName(edge->orientedReadIds[0].getReadId()) << ','
                 << reads->getReadName(edge->orientedReadIds[1].getReadId()) << ','
                 << (isSameStrand ? "Yes" : "No") << '\n';
        }
    }
}


bool Assembler::passesReadGraph2Criteria(const AlignmentInfo& info) const{
    const auto trims = info.computeTrim();
    const auto trim = max(trims.first, trims.second);

    // If this alignment doesn't pass the thresholds, skip it
    if (info.minAlignedFraction() < assemblerInfo->automatedAlignedFractionThreshold){
        return false;
    }
    if (info.markerCount < assemblerInfo->automatedMarkerCountThreshold){
        return false;
    }
    if (info.maxDrift > assemblerInfo->automatedMaxDriftThreshold){
        return false;
    }
    if (info.maxSkip > assemblerInfo->automatedMaxSkipThreshold){
        return false;
    }
    if (trim > assemblerInfo->automatedMaxTrimThreshold){
        return false;
    }

    return true;
}


void Assembler::setReadGraph2Criteria(
        double markerCountPercentile,
        double alignedFractionPercentile,
        double maxSkipPercentile,
        double maxDriftPercentile,
        double maxTrimPercentile
){
    const bool debug = false;

    // Initialize histograms for measuring alignedFraction, markerCount, maxDrift, and maxSkip distributions
    Histogram2 alignedFractionHistogram(0, 1, 100, false, false, true);
    Histogram2 markerCountHistogram(0, 3000, 300, false, false, true);
    Histogram2 maxDriftHistogram(0, 100, 100, false, false, true);
    Histogram2 maxSkipHistogram(0, 100, 100, false, false, true);
    Histogram2 maxTrimHistogram(0, 100, 100, false, false, true);

    ofstream alignmentInfoCsv;
    if (debug) {
        alignmentInfoCsv.open("AlignmentInfo.csv");
        alignmentInfoCsv
                << "readId0" << ','
                << "readId1" << ','
                << "minAlignedFraction" << ','
                << "markerCount" << ','
                << "maxDrift" << ','
                << "maxSkip" << ','
                << "trim" << '\n';
    }

    // Sample all available alignments that pass the initial permissive criteria
    for (size_t i=0; i<alignmentData.size(); i++){
        const auto info = alignmentData[i].info;
        const auto trims = info.computeTrim();
        const auto trim = max(trims.first, trims.second);

        alignedFractionHistogram.update(info.minAlignedFraction());
        markerCountHistogram.update(info.markerCount);
        maxDriftHistogram.update(info.maxDrift);
        maxSkipHistogram.update(info.maxSkip);
        maxTrimHistogram.update(trim);

        if (debug) {
            alignmentInfoCsv << alignmentData[i].readIds[0] << ','
                             << alignmentData[i].readIds[1] << ','
                             << info.minAlignedFraction() << ','
                             << info.markerCount << ','
                             << info.maxDrift << ','
                             << info.maxSkip << ','
                             << trim << '\n';
        }
    }

    // Minimums
    assemblerInfo->automatedAlignedFractionThreshold =
            alignedFractionHistogram.thresholdByCumulativeProportion(alignedFractionPercentile);
    assemblerInfo->automatedMarkerCountThreshold =
            markerCountHistogram.thresholdByCumulativeProportion(markerCountPercentile);

    // Maximums use (1 - percentile)
    assemblerInfo->automatedMaxDriftThreshold =
            maxDriftHistogram.thresholdByCumulativeProportion(1 - maxDriftPercentile);
    assemblerInfo->automatedMaxSkipThreshold =
            maxSkipHistogram.thresholdByCumulativeProportion(1 - maxSkipPercentile);
    assemblerInfo->automatedMaxTrimThreshold =
            maxTrimHistogram.thresholdByCumulativeProportion(1 - maxTrimPercentile);

    cout << "Selected thresholds automatically for the following parameters:\n\t"
         << "alignedFraction:\t" << assemblerInfo->automatedAlignedFractionThreshold << "\n\t"
         << "markerCount:\t\t" << assemblerInfo->automatedMarkerCountThreshold << "\n\t"
         << "maxDrift:\t\t" << assemblerInfo->automatedMaxDriftThreshold << "\n\t"
         << "maxSkip:\t\t" << assemblerInfo->automatedMaxSkipThreshold << "\n\t"
         << "maxTrim:\t\t" << assemblerInfo->automatedMaxTrimThreshold << "\n";

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

}


void Assembler::createReadGraph2(
    uint32_t maxAlignmentCount,
    double markerCountPercentile,
    double alignedFractionPercentile,
    double maxSkipPercentile,
    double maxDriftPercentile,
    double maxTrimPercentile)
{
    // First find thresholds based on the observed
    // distribution of alignment quality indicators
    setReadGraph2Criteria(
            markerCountPercentile,
            alignedFractionPercentile,
            maxSkipPercentile,
            maxDriftPercentile,
            maxTrimPercentile);

    vector<bool> keepAlignment(alignmentData.size(), false);

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

            // Discard each alignment if it does not pass the chosen thresholds
            if(not passesReadGraph2Criteria(info)){
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
