#include "AssembledSegment.hpp"
using namespace shasta;

#include "iterator.hpp"



void AssembledSegment::clear()
{
    assemblyGraphEdgeId = AssemblyGraph::invalidEdgeId;
    k = 0;
    vertexCount = 0;
    edgeCount = 0;

    vertexIds.clear();
    edgeIds.clear();
    vertexCoverage.clear();
    edgeCoverage.clear();

    vertexSequences.clear();
    vertexRepeatCounts.clear();

    edgeSequences.clear();
    edgeRepeatCounts.clear();
    edgeOverlappingBaseCounts.clear();

    vertexOffsets.clear();

    runLengthSequence.clear();
    repeatCounts.clear();
    rawSequence.clear();

    vertexRunLengthRange.clear();
    vertexRawRange.clear();
    edgeRunLengthRange.clear();
    edgeRawRange.clear();
}



void AssembledSegment::computeVertexOffsets()
{
    vertexOffsets.resize(vertexCount);
    vertexOffsets[0] = 0;

    for(size_t i=0; i<edgeCount; i++) {
        const uint8_t overlap = edgeOverlappingBaseCounts[i];
        if(overlap > 0) {
            SHASTA_ASSERT(edgeSequences[i].empty());
            SHASTA_ASSERT(edgeRepeatCounts[i].empty());
            vertexOffsets[i+1] = uint32_t(vertexOffsets[i] + k - overlap);
        } else {
            vertexOffsets[i+1] = uint32_t(vertexOffsets[i] + k + edgeSequences[i].size());
        }
    }
}



void AssembledSegment::computeVertexAssembledPortion()
{
    // Compute, for each vertex, the portion of vertex sequence that contributes
    // to the assembly. This is the portion that does not overlap a vertex with greater coverage.
    // (Break ties using vertex ids).
    // An edge with overlapping markers does not contribute to the assembly.
    // An edge with at least one intervening base contributes all of its bases
    // to the assembly.

    vertexAssembledPortion.resize(vertexCount);

    for(int i=0; i<int(vertexCount); i++) {

        // Check previous vertices.
        vertexAssembledPortion[i].first = 0;
        for(int j=i-1; j>=0; j--) {
            if(vertexOffsets[j]+k < vertexOffsets[i]) {
                break;
            }
            if(vertexCoverage[j]>vertexCoverage[i] ||
                (vertexCoverage[j]==vertexCoverage[i] && vertexIds[j]<vertexIds[i])) {
                vertexAssembledPortion[i].first =
                    vertexOffsets[j] + uint32_t(k) - vertexOffsets[i];
                break;
            }
        }

        // Check following vertices.
        vertexAssembledPortion[i].second = uint32_t(k);
        for(int j=i+1; j<int(vertexCount); j++) {
            if(vertexOffsets[i]+k < vertexOffsets[j]) {
                break;
            }
            if(vertexCoverage[j]>vertexCoverage[i] ||
                (vertexCoverage[j]==vertexCoverage[i] && vertexIds[j]<vertexIds[i])) {
                vertexAssembledPortion[i].second = vertexOffsets[j] - vertexOffsets[i];
                break;
            }
        }

        // Handle the case of a vertex that contributes nothing.
        if(vertexAssembledPortion[i].second <= vertexAssembledPortion[i].first) {
            vertexAssembledPortion[i].first = 0;
            vertexAssembledPortion[i].second = 0;
        }
        SHASTA_ASSERT(vertexAssembledPortion[i].second <= k);
    }
}



void AssembledSegment::assemble()
{
    // Figure out if we should store coverage data for assembled sequence.
    const bool storeCoverageData =
        vertexCoverageData.size() == vertexCount &&
        edgeCoverageData.size()   == edgeCount;

    vertexRunLengthRange.resize(vertexCount);
    vertexRawRange.resize(vertexCount);
    edgeRunLengthRange.resize(edgeCount);
    edgeRawRange.resize(edgeCount);

    for(size_t i=0; ; i++) {

        // Vertex.
        vertexRunLengthRange[i].first = uint32_t(runLengthSequence.size());
        vertexRawRange[i].first = uint32_t(rawSequence.size());
        for(uint32_t j=vertexAssembledPortion[i].first; j!=vertexAssembledPortion[i].second; j++) {
            const Base base = vertexSequences[i][j];
            const uint32_t repeatCount = vertexRepeatCounts[i][j];
            SHASTA_ASSERT(repeatCount > 0);
            runLengthSequence.push_back(base);
            repeatCounts.push_back(repeatCount);
            if(storeCoverageData) {
                assembledCoverageData.push_back(vertexCoverageData[i][j]);
            }
            for(uint32_t k=0; k!=repeatCount; k++) {
                rawSequence.push_back(base);
            }
        }
        vertexRunLengthRange[i].second = uint32_t(runLengthSequence.size());
        vertexRawRange[i].second = uint32_t(rawSequence.size());

        // This was the last vertex.
        if(i == edgeCount) {
            break;
        }

        // Edge.
        edgeRunLengthRange[i].first = uint32_t(runLengthSequence.size());
        edgeRawRange[i].first = uint32_t(rawSequence.size());
        if(edgeSequences[i].size() > 0) {
            for(uint32_t j=0; j!=uint32_t(edgeSequences[i].size()); j++) {
                const Base base = edgeSequences[i][j];
                const uint32_t repeatCount = edgeRepeatCounts[i][j];
                SHASTA_ASSERT(repeatCount > 0);
                runLengthSequence.push_back(base);
                repeatCounts.push_back(repeatCount);
                if(storeCoverageData) {
                    assembledCoverageData.push_back(edgeCoverageData[i][j]);
                }
                for(uint32_t k=0; k!=repeatCount; k++) {
                    rawSequence.push_back(base);
                }
            }
        }
        edgeRunLengthRange[i].second = uint32_t(runLengthSequence.size());
        edgeRawRange[i].second = uint32_t(rawSequence.size());
    }
}



// Write out details in html.
void AssembledSegment::writeHtml(ostream& html) const
{

    // Write a title.
    html <<
        "<h1>Assembly graph edge <a href="
        "'exploreAssemblyGraph?edgeId=" << assemblyGraphEdgeId <<
        "&maxDistance=6&detailed=on&sizePixels=1600&timeout=30'>" <<
        assemblyGraphEdgeId << "</a></h1>";

    // Assembled run-length sequence.
    html << "<p>Assembled run-length sequence (" << runLengthSequence.size() <<
        " bases):<br><span style='font-family:courier'>";
    copy(runLengthSequence.begin(), runLengthSequence.end(),
        ostream_iterator<Base>(html));
    html << "<br>";
    const uint32_t maxRepeatCount =
        *max_element(repeatCounts.begin(), repeatCounts.end());
    for(size_t j=0; j<repeatCounts.size(); j++) {
        const uint32_t repeatCount = repeatCounts[j];
        html << repeatCount%10;
    }
    if(maxRepeatCount >= 10) {
        html << "<br>";
        for(size_t j=0; j<repeatCounts.size(); j++) {
            const uint32_t repeatCount = repeatCounts[j];
            const uint32_t digit = (repeatCount/10) % 10;
            if(digit == 0) {
                html << "&nbsp;";
            } else {
                html << digit;
            }
        }
    }
    if(maxRepeatCount >= 100) {
        html << "<br>";
        for(size_t j=0; j<repeatCounts.size(); j++) {
            const uint32_t repeatCount = repeatCounts[j];
            if(repeatCount >= 1000) {
                html << "*";
            } else {
                const uint32_t digit = (repeatCount/100) % 10;
                if(digit == 0) {
                    html << "&nbsp;";
                } else {
                    html << digit;
                }
            }
        }
    }
    html << "</span>";

    // Assembled raw sequence.
    html << "<p>Assembled raw sequence (" << rawSequence.size() <<
        " bases):<br><span style='font-family:courier'>";
    copy(rawSequence.begin(), rawSequence.end(),
        ostream_iterator<Base>(html));
    html << "</span>";




    // Write a table with a row for each marker graph vertex or edge
    // in the marker graph chain.
    html <<
        "<p>This vertex of the assembly graph corresponds to a chain of " <<
        edgeIds.size() << " edges in the marker graph. "
        "The table below shows consensus sequences "
        "for the vertices and edges of this chain of the marker graph. "
        "All vertex and edge ids in the table refer to the marker graph."
        "<p><table><tr>"
        "<th rowspan=2>Vertex<br>or<br>edge"
        "<th rowspan=2>Id"
        "<th rowspan=2>Coverage"
        "<th colspan=4>Run-length"
        "<th colspan=3>Raw"
        "<tr>"
        "<th>Offset"
        "<th>Begin"
        "<th>End"
        "<th>Sequence"
        "<th>Begin"
        "<th>End"
        "<th>Sequence";
    const string urlPrefix = "exploreMarkerGraph?vertexId=";
    const string urlSuffix =
        "&maxDistance=5"
        "&detailed=on"
        "&minCoverage=3"
        "&minConsensus=3"
        "&sizePixels=600&timeout=10"
        "&useStoredConnectivity=on"
        "&showVertexId=on"
        "&showAssembledSequence=on";
    for(size_t i=0; ; i++) {

        // Vertex.
        const MarkerGraph::VertexId vertexId = vertexIds[i];
        const string url = urlPrefix + to_string(vertexId) + urlSuffix;
        const vector<Base>& vertexSequence = vertexSequences[i];
        const vector<uint32_t>& vertexRepeatCount = vertexRepeatCounts[i];
        // const uint32_t maxVertexRepeatCount =
        //     *std::max_element(vertexRepeatCount.begin(), vertexRepeatCount.end());
        html <<
             "<tr><td>Vertex" <<
            "<td class=centered><a href='" << url << "'>" << vertexId << "</a>"
            "<td class=centered>" << vertexCoverage[i] <<
            "<td class=centered>" << vertexOffsets[i] <<
            "<td class=centered>" << vertexRunLengthRange[i].first <<
            "<td class=centered>" << vertexRunLengthRange[i].second <<
            "<td style='font-family:courier'>";
        for(size_t j=0; j<vertexSequence.size(); j++) {
            if(j==vertexAssembledPortion[i].first &&
                vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                html << "<span style='background-color:LightGreen'>";
            }
            html << vertexSequence[j];
            if(j==vertexAssembledPortion[i].second-1  &&
                vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                html << "</span>";
            }
        }
        html << "<br>";
        for(size_t j=0; j<vertexSequence.size(); j++) {
            const uint32_t repeatCount = vertexRepeatCount[j];
            if(j==vertexAssembledPortion[i].first &&
                vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                html << "<span style='background-color:LightGreen'>";
            }
            if(repeatCount < 10) {
                html << repeatCount % 10;
            } else {
                html << "*";
            }
            if(j==vertexAssembledPortion[i].second-1 &&
                vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                html << "</span>";
            }
        }
        html <<
            "<td class=centered>" << vertexRawRange[i].first <<
            "<td class=centered>" << vertexRawRange[i].second <<
            "<td style='font-family:courier'>";
        for(size_t j=0; j<vertexSequence.size(); j++) {
            if(j==vertexAssembledPortion[i].first &&
                vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                html << "<span style='background-color:LightGreen'>";
            }
            const Base b = vertexSequence[j];
            const uint32_t repeatCount = vertexRepeatCount[j];
            for(uint32_t k=0; k<repeatCount; k++) {
                html << b;
            }
            if(j==vertexAssembledPortion[i].second-1 &&
                vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                html << "</span>";
            }
        }

        // This was the last vertex.
        if(i == edgeCount) {
            break;
        }



        // Edge.
        const MarkerGraph::EdgeId edgeId = edgeIds[i];
        // const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
        const string sourceUrl = urlPrefix + to_string(vertexIds[i]) + urlSuffix;
        const string targetUrl = urlPrefix + to_string(vertexIds[i+1]) + urlSuffix;
        const vector<Base>& edgeSequence = edgeSequences[i];
        const vector<uint32_t>& edgeRepeatCount = edgeRepeatCounts[i];
        const size_t edgeSequenceLength = edgeSequence.size();
        SHASTA_ASSERT(edgeRepeatCount.size() == edgeSequenceLength);
        // const uint32_t maxEdgeRepeatCount =
        //    *std::max_element(edgeRepeatCount.begin(), edgeRepeatCount.end());
        // SHASTA_ASSERT(maxEdgeRepeatCount < 10);  // For now. Add additional code when this fails.
        html <<
            "<tr><td>Edge<td class=centered>" << edgeId <<
            "<td class=centered>" << edgeCoverage[i] <<
            "<td class=centered>" <<
            "<td class=centered>" << edgeRunLengthRange[i].first <<
            "<td class=centered>" << edgeRunLengthRange[i].second <<
            "<td style='font-family:courier'>";
        for(size_t j=0; j<edgeSequenceLength; j++) {
            if(edgeSequenceLength>2*k && j==k) {
                html << "<span style='background-color:LightGreen'>";
            }
            html << edgeSequence[j];
            if(edgeSequenceLength>2*k && j==edgeSequenceLength-k-1) {
                html << "</span>";
            }
        }
        html << "<br>";
        for(size_t j=0; j<edgeSequenceLength; j++) {
            const uint32_t repeatCount = edgeRepeatCount[j];
            if(edgeSequenceLength>2*k && j==k) {
                html << "<span style='background-color:LightGreen'>";
            }
            if(repeatCount < 10) {
                html << repeatCount % 10;
            } else {
                html << "*";
            }
            if(edgeSequenceLength>2*k && j==edgeSequenceLength-k-1) {
                html << "</span>";
            }
        }
        html <<
            "<td class=centered>" << edgeRawRange[i].first <<
            "<td class=centered>" << edgeRawRange[i].second <<
            "<td style='font-family:courier'>";
        for(size_t j=0; j<edgeSequenceLength; j++) {
            const Base b = edgeSequence[j];
            const uint32_t repeatCount = edgeRepeatCount[j];
            if(edgeSequenceLength>2*k && j==k) {
                html << "<span style='background-color:LightGreen'>";
            }
            for(uint32_t k=0; k<repeatCount; k++) {
                html << b;
            }
            if(edgeSequenceLength>2*k && j==edgeSequenceLength-k-1) {
                html << "</span>";
            }
        }
     }
}
