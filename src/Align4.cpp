#include "Align4.hpp"
#include "html.hpp"
using namespace shasta;

#include "algorithm.hpp"
#include "fstream.hpp"



// Align two arbitrary sequences  using alignment method 4.
// If debug is true, detailed output to html is produced.
// Otherwise, html is not used.
void shasta::align4(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const Align4Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug)
{
    switch(options.m) {
    case 1:
        align4<1>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    case 2:
        align4<2>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    case 3:
        align4<3>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    case 4:
        align4<4>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    default:
        SHASTA_ASSERT(0);
    }
}



// Version templated on m, the number of markers that define
// a "feature" used in the alignment.
template<uint64_t m> void shasta::align4(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const Align4Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug)
{
    Align4<m> graph(markers0, markers1,
        options, alignment, alignmentInfo,
        debug);
}



template<uint64_t m> shasta::Align4<m>::Align4(
    const Sequence& sequence0,
    const Sequence& sequence1,
    const Align4Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug)
{
    // Check that we are in the templated version consistent with
    /// the options.
    SHASTA_ASSERT(options.m == m);

    if(debug) {
        cout << "Computing a marker alignment of two sequences with " <<
            sequence0.size() << " and " << sequence1.size() << " markers." << endl;
    }

    // Create the FeatureMap for sequence0.
    const uint64_t inverseLoadFactor = 2;
    FeatureMap featureMap0(inverseLoadFactor * sequence0.size());
    fillFeatureMap(sequence0, featureMap0);

    // Create the alignment matrix.
    fillAlignmentMatrix(featureMap0, sequence1,
        uint32_t(sequence0.size()-(m-1)),
        uint32_t(sequence1.size()-(m-1)),
        uint32_t(options.deltaX), uint32_t(options.deltaY));

    if(debug) {
        writeMatrixCsv(alignmentMatrix, "Align4-Matrix.csv");
    }
}



template<uint64_t m> void shasta::Align4<m>::fillFeatureMap(
    const Sequence& sequence,
    FeatureMap& featureMap)
{
    SHASTA_ASSERT(featureMap.empty());

    // If the sequence is too short, it has no features.
    if(sequence.size() < m) {
        return;
    }

    // Start with the first m kmerIds.
    Feature feature;
    uint64_t nextOrdinal = 0;
    for(uint64_t i=0; i<m; i++) {
        feature[i] = sequence[nextOrdinal++].kmerId;
    }

    // Add the features.
    for(uint32_t ordinal=0; ; ordinal++) {
        featureMap.insert(make_pair(feature, ordinal));

        // Check if done.
        if(nextOrdinal == sequence.size()) {
            break;
        }

        // Shift by one.
        for(uint64_t i=1; i<m; i++) {
            feature[i-1] = feature[i];
        }
        feature.back() = sequence[nextOrdinal++].kmerId;
    }

    if(sequence.size() >= m) {
        SHASTA_ASSERT(featureMap.size() == sequence.size() + 1 - m);
    } else {
        SHASTA_ASSERT(featureMap.empty());
    }
}



template<uint64_t m> void shasta::Align4<m>::fillAlignmentMatrix(
    const FeatureMap& featureMap0,
    const Sequence& sequence1,
    int32_t nx,
    int32_t ny,
    int32_t deltaX,
    int32_t deltaY)
{
    // Start with the first m kmerIds.
    Feature feature;
    uint64_t nextOrdinal1 = 0;
    for(uint64_t i=0; i<m; i++) {
        feature[i] = sequence1[nextOrdinal1++].kmerId;
    }

    // Add the features.
    for(int32_t y=0; ; y++) {

        // Look up this feature in the feature map for sequence 0.
        typename FeatureMap::const_iterator begin, end;
        tie(begin, end) = featureMap0.equal_range(feature);

        // Loop over all the hits. Each hit generates an
        // entry in the alignment matrix.
        for(auto it=begin; it!=end; ++it) {
            const int32_t x = it->second;
            const int32_t X = x + y;
            const int32_t Y = y + nx - 1 - x;
            const int32_t iX = X / deltaX;
            const int32_t iY = Y / deltaY;
            AlignmentMatrixEntry alignmentMatrixEntry;
            alignmentMatrixEntry.xy = make_pair(x, y);
            alignmentMatrixEntry.XY = make_pair(X, Y);
            alignmentMatrixEntry.isNearLeft = (x < deltaX) ? 1 : 0;
            alignmentMatrixEntry.isNearTop  = (y < deltaX) ? 1 : 0;
            alignmentMatrixEntry.isNearRight = (nx-1-x < deltaX) ? 1 : 0;
            alignmentMatrixEntry.isNearBottom = (ny-1-y < deltaX) ? 1 : 0;
            alignmentMatrix.insert(make_pair(Coordinates(iX, iY), alignmentMatrixEntry));
        }

        // Check if done.
        if(nextOrdinal1 == sequence1.size()) {
            break;
        }

        // Shift by one.
        for(uint64_t i=1; i<m; i++) {
            feature[i-1] = feature[i];
        }
        feature.back() = sequence1[nextOrdinal1++].kmerId;
    }
}


template<uint64_t m> void shasta::Align4<m>::writeMatrixCsv(
    const AlignmentMatrix& alignmentMatrix,
    const string& fileName)
{
    ofstream csv(fileName);
    csv << "x,y,X,Y,iX,iY,Near left,Near top,Near right,Near bottom\n";

    for(const auto& p: alignmentMatrix) {
        const Coordinates& iXY = p.first;
        const AlignmentMatrixEntry& alignmentMatrixEntry = p.second;
        csv <<
            alignmentMatrixEntry.xy.first << "," <<
            alignmentMatrixEntry.xy.second << "," <<
            alignmentMatrixEntry.XY.first << "," <<
            alignmentMatrixEntry.XY.second << "," <<
            iXY.first << "," <<
            iXY.second << "," <<
            int(alignmentMatrixEntry.isNearLeft) << "," <<
            int(alignmentMatrixEntry.isNearTop) << "," <<
            int(alignmentMatrixEntry.isNearRight) << "," <<
            int(alignmentMatrixEntry.isNearBottom) << "\n";
    }
}



