#ifndef SHASTA_LOCAL_MARKER_GRAPH_REQUEST_PARAMETERS_HPP
#define SHASTA_LOCAL_MARKER_GRAPH_REQUEST_PARAMETERS_HPP

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
};

#endif
