#include "Align4.hpp"
#include "html.hpp"
using namespace shasta;

#include "algorithm.hpp"



// Align two arbitrary sequences  using alignment method 4.
// If debug is true, detailed output to html is produced.
// Otherwise, html is not used.
void shasta::align4(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const Align4Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug,
    ostream& html)
{
    switch(options.m) {
    case 1:
        align4<1>(markers0, markers1, options, alignment, alignmentInfo, debug, html);
        return;
    case 2:
        align4<2>(markers0, markers1, options, alignment, alignmentInfo, debug, html);
        return;
    case 3:
        align4<3>(markers0, markers1, options, alignment, alignmentInfo, debug, html);
        return;
    case 4:
        align4<4>(markers0, markers1, options, alignment, alignmentInfo, debug, html);
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
    bool debug,
    ostream& html)
{
    Align4<m> graph(markers0, markers1,
        options, alignment, alignmentInfo,
        debug, html);
}



template<uint64_t m> shasta::Align4<m>::Align4(
    const Sequence& sequence0,
    const Sequence& sequence1,
    const Align4Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug,
    ostream& html)
{
    // Check that we are in the templated version consistent with
    /// the options.
    SHASTA_ASSERT(options.m == m);

    if(debug) {
        html << "<p>Computing a marker alignment of two sequences with " <<
            sequence0.size() << " and " << sequence1.size() << " markers." << endl;
    }

    // Create the FeatureMap for sequence0.
    const uint64_t inverseLoadFactor = 2;
    FeatureMap featureMap0(inverseLoadFactor * sequence0.size());
    fillFeatureMap(sequence0, featureMap0);

    // Create the alignment matrix.
    AlignmentMatrix alignmentMatrix(2*max(sequence0.size(), sequence1.size()));
    alignmentMatrix.max_load_factor(4);
    const uint64_t cellSizeX = options.deltaX;
    const uint64_t cellSizeY = options.deltaY;
    fillAlignmentMatrix(featureMap0, sequence1,
        sequence0.size(), cellSizeX, cellSizeY, alignmentMatrix);

    // Find occupied cells.
    std::set<Cell> cells;
    for(const auto& p: alignmentMatrix) {
        cells.insert(p.first);
    }

    if(debug) {
        html << "<h2>Alignment matrix</h2>"
            "<p>The alignment matrix has " << alignmentMatrix.size() << " entries.\n" <<
            "<br>The number of occupied cells is " << cells.size() << ".\n";
        // write(alignmentMatrix, html);
        html << "<br>\n";
        writeSvg(cells, sequence0.size(), sequence1.size(), cellSizeX, cellSizeY, html);
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
    uint64_t nx,
    uint64_t cellSizeX,
    uint64_t cellSizeY,
    AlignmentMatrix& alignmentMatrix
)
{
    // Start with the first m kmerIds.
    Feature feature;
    uint64_t nextOrdinal1 = 0;
    for(uint64_t i=0; i<m; i++) {
        feature[i] = sequence1[nextOrdinal1++].kmerId;
    }

    // Add the features.
    for(uint32_t ordinal1=0; ; ordinal1++) {

        // Look up this feature in the feature map for sequence 0.
        typename FeatureMap::const_iterator begin, end;
        tie(begin, end) = featureMap0.equal_range(feature);
        for(auto it=begin; it!=end; ++it) {
            const uint32_t ordinal0 = it->second;
            const uint64_t X = (ordinal0 + ordinal1) / cellSizeX;
            const uint64_t Y = (ordinal1 + nx - 1 - ordinal0) / cellSizeY;
            alignmentMatrix.insert(make_pair(Cell(X, Y), OrdinalPair(ordinal0, ordinal1)));
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


template<uint64_t m> void shasta::Align4<m>::write(
    const AlignmentMatrix& alignmentMatrix, ostream& html)
{
    html << "<table>\n"
        "<tr><th>iX<th>iY<th>Ordinal0<th>Ordinal1\n";


    for(const auto& p: alignmentMatrix) {
        const Cell& cell = p.first;
        const OrdinalPair& ordinalPair = p.second;
        html <<
            "<tr>"
            "<td class=centered>" << cell.first <<
            "<td class=centered>" << cell.second <<
            "<td class=centered>" << ordinalPair.first <<
            "<td class=centered>" << ordinalPair.second << "\n";
    }

    html << "</table>\n";

}



template<uint64_t m> void shasta::Align4<m>::writeSvg(
    const std::set<Cell>& cells,
    uint64_t nx,
    uint64_t ny,
    uint64_t cellSizeX,
    uint64_t cellSizeY,
    ostream& html)
{
    const uint64_t svgSizePixels = 800;
    const int64_t squareSize = nx + ny -2;
    const int64_t borderSize = 2;

    html << "<svg width='" << svgSizePixels << "' height='" << svgSizePixels <<
        "' viewbox='" << -borderSize <<
        " " << -borderSize <<
        " " << squareSize + borderSize <<
        " " << squareSize + borderSize <<
        "' style='background-color:Beige'>\n";



    for(const Cell& cell: cells) {
        const uint64_t iX = cell.first;
        const uint64_t iY = cell.second;
        const uint64_t X = iX * cellSizeX;
        const uint64_t Y = iY * cellSizeY;
        const int64_t centerX = X + cellSizeX / 2;
        const int64_t centerY = Y + cellSizeY / 2;
        const int64_t centerx = (centerX - centerY + int64_t(nx) - 1) / 2;
        const int64_t centery = (centerX + centerY - int64_t(nx) + 1) / 2;
        html << "<rect x='" << X <<
            "' y='" << Y <<
            "' width='" << cellSizeX <<
            "' height='" << cellSizeY <<
            "' fill='Grey'><title>"
            "Cell (iX,iY)=(" << iX << "," << iY <<
            ") centered at (X,Y)=(" << centerX <<
            "," << centerY <<
            "), (x,y)=(" << centerx << "," << centery << ")"
            "</title></rect>\n";
    }



    html << "</svg>\n";
}

