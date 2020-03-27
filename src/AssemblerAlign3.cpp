#include "Assembler.hpp"
using namespace shasta;

#include <numeric>


// Align two oriented reads using SeqAn banded alignment.
void Assembler::alignOrientedReads3(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo)
{
    const bool debug = true;

    // Get the markers for the two oriented reads.
    const span<CompressedMarker> markers0 = markers[orientedReadId0.getValue()];
    const span<CompressedMarker> markers1 = markers[orientedReadId1.getValue()];
    const uint64_t markerCount0 = markers0.size();
    const uint64_t markerCount1 = markers1.size();

    // Get the markers sorted by kmerId.
    vector<MarkerWithOrdinal> markersSortedByKmerId0;
    vector<MarkerWithOrdinal> markersSortedByKmerId1;
    getMarkersSortedByKmerId(orientedReadId0, markersSortedByKmerId0);
    getMarkersSortedByKmerId(orientedReadId1, markersSortedByKmerId1);
    SHASTA_ASSERT(markersSortedByKmerId0.size() == markerCount0);
    SHASTA_ASSERT(markersSortedByKmerId1.size() == markerCount1);

    // Some iterators we will need.
    using MarkerIterator = vector<MarkerWithOrdinal>::const_iterator;
    const MarkerIterator begin0 = markersSortedByKmerId0.begin();
    const MarkerIterator end0   = markersSortedByKmerId0.end();
    const MarkerIterator begin1 = markersSortedByKmerId1.begin();
    const MarkerIterator end1   = markersSortedByKmerId1.end();

    // For an element in the alignment matrix at coordinates (m0, m1)
    // we define the offset as m0-m1.
    // Compute the range for possible values of the offset.
    const int64_t offsetBegin = -int64_t(markerCount1-1);
    const int64_t offsetEnd = int64_t(markerCount0);



    // Compute a histogram of the offset of alignment matrix elements that are 1.
    // Use a joint loop over the markers sorted by KmerId, looking for common KmerId's.
    vector<uint64_t> histogram(offsetEnd-offsetBegin, 0);
    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->kmerId < it1->kmerId) {
            ++it0;
        } else if(it1->kmerId < it0->kmerId) {
            ++it1;
        } else {

            // We found a common k-mer id.
            const KmerId kmerId = it0->kmerId;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer
            // in markersSortedByKmerId0 and markersSortedByKmerId0.
            MarkerIterator it0Begin = it0;
            MarkerIterator it1Begin = it1;
            MarkerIterator it0End = it0Begin;
            MarkerIterator it1End = it1Begin;
            while(it0End!=end0 && it0End->kmerId==kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->kmerId==kmerId) {
                ++it1End;
            }

            // Loop over pairs in the two streaks.
            for(MarkerIterator jt0=it0Begin; jt0!=it0End; ++jt0) {
                for(MarkerIterator jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const int64_t offset = int64_t(jt0->ordinal) - int64_t(jt1->ordinal);
                    ++histogram[offset - offsetBegin];
                }
            }

            // Continue joint loop over k-mers.
            it0 = it0End;
            it1 = it1End;
        }
    }


    // Create a smoothed offset histogram.
    const double nDrift = 20.;  // Distance over which we expect unit diffusion relative drift.
    const int64_t halfWidth = int64_t(sqrt(double(min(markerCount0, markerCount1)) / nDrift));
    vector<uint64_t> smoothedHistogram(histogram.size(), 0);
    uint64_t sum = std::accumulate(histogram.begin(), histogram.begin() + (2*halfWidth+1), 0);
    for(uint64_t i=halfWidth; i<smoothedHistogram.size()-halfWidth-1; i++) {
        smoothedHistogram[i] = sum;
        sum -= histogram[i-halfWidth];
        sum += histogram[i+halfWidth];
    }



    // Write out the offset histogram.
    if(debug) {
        ofstream csv("OffsetHistogram.csv");
        csv << "Offset,Frequency,Smoothed frequency,Range,Smoothed frequency/Range\n";
        for(uint64_t i=0; i<histogram.size(); i++) {
            const uint64_t frequency = histogram[i];
            const uint64_t smoothedFrequency = smoothedHistogram[i];
            const int64_t offset = offsetBegin + int64_t(i);
            const int64_t range0 =
                min(int64_t(markerCount0), int64_t(markerCount1) + offset) -
                max(int64_t(0L), offset);
            const int64_t range1 =
                min(int64_t(markerCount0-offset), int64_t(markerCount1)) -
                max(-offset, int64_t(0L));
            SHASTA_ASSERT(range0 == range1);
            csv <<
                offset << "," <<
                frequency << "," <<
                smoothedFrequency << "," <<
                range0 << "," <<
                double(smoothedFrequency) / double(range0) << "\n";
        }
    }
}
