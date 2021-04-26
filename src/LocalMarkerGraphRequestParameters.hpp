#ifndef SHASTA_LOCAL_MARKER_GRAPH_REQUEST_PARAMETERS_HPP
#define SHASTA_LOCAL_MARKER_GRAPH_REQUEST_PARAMETERS_HPP

#include "ReadId.hpp"
#include <map>

namespace shasta {
    class LocalMarkerGraphRequestParameters;
}


// Class describing the parameters in the form
// in the local marker graph page.
class shasta::LocalMarkerGraphRequestParameters {
public:

    uint64_t vertexId;
    bool vertexIdIsPresent;
    uint32_t maxDistance;
    bool maxDistanceIsPresent;
    string layoutMethod;    // dotLr, dotTb, or sfdp
    bool useWeakEdges;
    bool usePrunedEdges;
    bool useSuperBubbleEdges;
    bool useLowCoverageCrossEdges;

    // Vertex and edge label control:
    // 0 = no label
    // 1 = terse label
    // 2 = verbose tabel
    uint64_t vertexLabels;
    uint64_t edgeLabels;
    string vertexLabelsString() const;
    string edgeLabelsString() const;

    uint64_t minVertexCoverage;
    bool minVertexCoverageIsPresent;
    uint64_t minEdgeCoverage;
    bool minEdgeCoverageIsPresent;

    uint32_t sizePixels;
    bool sizePixelsIsPresent;
    double vertexScalingFactor;
    bool vertexScalingFactorIsPresent;
    string vertexScalingFactorString() const;
    double edgeThicknessScalingFactor;
    bool edgeThicknessScalingFactorIsPresent;
    string edgeThicknessScalingFactorString() const;
    double arrowScalingFactor;
    bool arrowScalingFactorIsPresent;
    string arrowScalingFactorString() const;
    int timeout;
    bool timeoutIsPresent;

    string vertexColoring;
    string edgeColoring;
    uint64_t vertexRedCoverage;
    uint64_t vertexGreenCoverage;
    uint64_t edgeRedCoverage;
    uint64_t edgeGreenCoverage;

    void writeForm(ostream&, uint64_t vertexCount) const;
    bool hasMissingRequiredParameters() const;

    string url() const;
    string urlForVertex(uint64_t newVertexId) const;



    // Data and  functions related to highlighted oriented reads.

    // The string describing highlighted oriented reads, as entered in the form.
    string highlightedOrientedReadsString;

    // The S and V values, in HSV color space, for colors used to highlight oriented reads.
    const double S = 0.7;
    const double V = 1.;

    // The highlighted oriented reads, and the corresponding Hue values.
    // Hue values are stored as defined in Graphviz - that is, with a range 0 to 1,
    // instead of the more common 0 to 360 degrees.
    std::map<OrientedReadId, double> highlightedOrientedReads;

    // This parses highlightedOrientedReadsString and creates
    // highlightedOrientedReads. Each oriented read is assigned a hue
    // via hashing of the OrientedReadId. This way, an oriented read
    // is always highlighted in the same color.
    void parseHighlightedOrientedReads();

};

#endif
