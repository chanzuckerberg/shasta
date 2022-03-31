#ifdef SHASTA_HTTP_SERVER

// Shasta.
#include "LocalMarkerGraph.hpp"
#include "ConsensusCaller.hpp"
#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "orderPairs.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard libraries.
#include "fstream.hpp"



// Write the graph in Graphviz format.
void LocalMarkerGraph::write(
    const string& fileName,
    const LocalMarkerGraphRequestParameters& localMarkerGraphRequestParameters) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, localMarkerGraphRequestParameters);
}
void LocalMarkerGraph::write(
    ostream& s,
    const LocalMarkerGraphRequestParameters& localMarkerGraphRequestParameters) const
{
    Writer writer(*this, localMarkerGraphRequestParameters);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&LocalMarkerGraphVertex::vertexId, *this));
}

LocalMarkerGraph::Writer::Writer(
    const LocalMarkerGraph& graph,
    const LocalMarkerGraphRequestParameters& parameters) :
    LocalMarkerGraphRequestParameters(parameters),
    graph(graph)
{
}



// Vertex and edge colors.
const string LocalMarkerGraph::Writer::vertexColorZeroDistance                          = "#6666ff";
const string LocalMarkerGraph::Writer::vertexColorIntermediateDistance                  = "#00ccff";
const string LocalMarkerGraph::Writer::vertexColorMaxDistance                           = "#66ffff";
const string LocalMarkerGraph::Writer::edgeArrowColorRemovedDuringTransitiveReduction   = "#ff0000";
const string LocalMarkerGraph::Writer::edgeArrowColorRemovedDuringPruning               = "#ff00ff";
const string LocalMarkerGraph::Writer::edgeArrowColorRemovedDuringSuperBubbleRemoval    = "#009900";
const string LocalMarkerGraph::Writer::edgeArrowColorRemovedAsLowCoverageCrossEdge      = "#c0c000";
const string LocalMarkerGraph::Writer::edgeArrowColorRemovedWhileSplittingSecondaryEdges      = "#ff0000";
const string LocalMarkerGraph::Writer::edgeArrowColorNotRemovedNotAssembled             = "#fcba03";
const string LocalMarkerGraph::Writer::edgeArrowColorNotRemovedAssembled                = "#000000";
const string LocalMarkerGraph::Writer::edgeLabelColorRemovedDuringTransitiveReduction   = "#ff9999";
const string LocalMarkerGraph::Writer::edgeLabelColorRemovedDuringPruning               = "#c03280";
const string LocalMarkerGraph::Writer::edgeLabelColorRemovedDuringSuperBubbleRemoval    = "#99ff99";
const string LocalMarkerGraph::Writer::edgeLabelColorRemovedAsLowCoverageCrossEdge      = "#e0e000";
const string LocalMarkerGraph::Writer::edgeLabelColorNotRemovedNotAssembled             = "#996600";
const string LocalMarkerGraph::Writer::edgeLabelColorNotRemovedAssembled                = "#999999";



string LocalMarkerGraph::Writer::vertexColor(const LocalMarkerGraphVertex& vertex) const
{
    if(vertexColoring == "none") {
        return "black";

    } else if(vertexColoring == "byCoverage") {

        const uint64_t coverage = vertex.markerInfos.size();
        const array<uint64_t, 2> strandCoverage =  vertex.strandCoverage();
        const uint64_t minStrandCoverage = min(strandCoverage[0], strandCoverage[1]);
        double hTotal =
            (double(coverage) - double(vertexRedCoverage)) /
            (double(vertexGreenCoverage) - double(vertexRedCoverage));
        double hStrand =
            (double(minStrandCoverage) - double(vertexRedCoveragePerStrand)) /
            (double(vertexGreenCoveragePerStrand) - double(vertexRedCoveragePerStrand));
        // Color by the worst of the three.
        double h = min(hTotal, hStrand);
        h = max(h, 0.);
        h = min(h, 1.);
        return to_string(h/3.) + ",1.,0.9";

    } else if(vertexColoring == "byDistance") {

        if(vertex.distance == 0) {
            return vertexColorZeroDistance;
        } else if(vertex.distance == maxDistance) {
            return vertexColorMaxDistance;
        } else {
            return vertexColorIntermediateDistance;
        }

    } else {
        throw runtime_error("Invalid vertex coloring " + vertexColoring);
    }
}



string LocalMarkerGraph::Writer::edgeArrowColor(const LocalMarkerGraphEdge& edge) const
{

    if(edgeColoring == "none") {
        return "black";

    } else if(edgeColoring == "byCoverage") {

        const uint64_t coverage = edge.coverage();
        const array<uint64_t, 2> strandCoverage =  edge.strandCoverage();
        const uint64_t minStrandCoverage = min(strandCoverage[0], strandCoverage[1]);
        double hTotal =
            (double(coverage) - double(edgeRedCoverage)) /
            (double(edgeGreenCoverage) - double(edgeRedCoverage));
        double hStrand =
            (double(minStrandCoverage) - double(edgeRedCoveragePerStrand)) /
            (double(edgeGreenCoveragePerStrand) - double(edgeRedCoveragePerStrand));
        // Color by the worst of the three.
        double h = min(hTotal, hStrand);
        h = max(h, 0.);
        h = min(h, 1.);
        return to_string(h/3.) + ",1.,0.9";

    } else if(edgeColoring == "byFlags") {

        if(graph.assemblyMode == 2) {
            if(edge.wasRemovedWhileSplittingSecondaryEdges) {
                return edgeArrowColorRemovedWhileSplittingSecondaryEdges;
            } else {
                if(edge.wasAssembled) {
                    return edgeArrowColorNotRemovedAssembled;
                } else {
                    return edgeArrowColorNotRemovedNotAssembled;
                }
            }

        } else {

            if(edge.wasRemovedByTransitiveReduction) {
                return edgeArrowColorRemovedDuringTransitiveReduction;
            } else if(edge.wasPruned) {
                return edgeArrowColorRemovedDuringPruning;
            } else if (edge.isSuperBubbleEdge) {
                return edgeArrowColorRemovedDuringSuperBubbleRemoval;
            } else if (edge.isLowCoverageCrossEdge) {
                return edgeArrowColorRemovedAsLowCoverageCrossEdge;
            } else {
                if(edge.wasAssembled) {
                    return edgeArrowColorNotRemovedAssembled;
                } else {
                    return edgeArrowColorNotRemovedNotAssembled;
                }
            }
        }
    } else {
        throw runtime_error("Invalid edge coloring " + edgeColoring);
    }
}



string LocalMarkerGraph::Writer::edgeLabelColor(const LocalMarkerGraphEdge& edge) const
{
    if(edgeColoring == "none") {
        return "white";

    } else if(edgeColoring == "byCoverage") {

        const uint64_t coverage = edge.coverage();
        const array<uint64_t, 2> strandCoverage =  edge.strandCoverage();
        const uint64_t minStrandCoverage = min(strandCoverage[0], strandCoverage[1]);
        double hTotal =
            (double(coverage) - double(edgeRedCoverage)) /
            (double(edgeGreenCoverage) - double(edgeRedCoverage));
        double hStrand =
            (double(minStrandCoverage) - double(edgeRedCoveragePerStrand)) /
            (double(edgeGreenCoveragePerStrand) - double(edgeRedCoveragePerStrand));
        // Color by the worst of the three.
        double h = min(hTotal, hStrand);
        h = max(h, 0.);
        h = min(h, 1.);
        return to_string(h/3.) + ",1.,0.9";

    } else if(edgeColoring == "byFlags") {

            if(edge.wasRemovedByTransitiveReduction) {
            return edgeLabelColorRemovedDuringTransitiveReduction;
        } else if(edge.wasPruned) {
            return edgeLabelColorRemovedDuringPruning;
        } else if (edge.isSuperBubbleEdge) {
            return edgeLabelColorRemovedDuringSuperBubbleRemoval;
        } else if (edge.isLowCoverageCrossEdge) {
            return edgeLabelColorRemovedAsLowCoverageCrossEdge;
        } else {
            if(edge.wasAssembled) {
                return edgeLabelColorNotRemovedAssembled;
            } else {
                return edgeLabelColorNotRemovedNotAssembled;
            }
        }

    } else {
        throw runtime_error("Invalid edge coloring " + edgeColoring);
    }
}



void LocalMarkerGraph::writeColorLegendVerticesByDistance(ostream& html)
{
    html <<
        "<table>"
        "<tr><td>Zero distance<td style='width:50px;background-color:" <<
        Writer::vertexColorZeroDistance << "'>"
        "<tr><td>Intermediate distances<td style='width:50px;background-color:" <<
        Writer::vertexColorIntermediateDistance << "'>"
        "<tr><td>Maximum distance<td style='width:50px;background-color:" <<
        Writer::vertexColorMaxDistance << "'>"
        "</table>";
}



void LocalMarkerGraph::writeColorLegendEdgeArrowsByFlags(ostream& html)
{
    if(assemblyMode == 2) {
        html <<
            "<table>"
            "<tr><td>Removed while splitting secondary edges<td style='width:50px;background-color:" <<
            Writer::edgeArrowColorRemovedWhileSplittingSecondaryEdges << "'>"
            "<tr><td>Not removed, not assembled<td style='width:50px;background-color:" <<
            Writer::edgeArrowColorNotRemovedNotAssembled << "'>"
            "<tr><td>Not removed, assembled<td style='width:50px;background-color:" <<
            Writer::edgeArrowColorNotRemovedAssembled << "'>"
            "</table>";

    } else {
        html <<
            "<table>"
            "<tr><td>Removed during transitive reduction<td style='width:50px;background-color:" <<
            Writer::edgeArrowColorRemovedDuringTransitiveReduction << "'>"
            "<tr><td>Removed during pruning<td style='width:50px;background-color:" <<
            Writer::edgeArrowColorRemovedDuringPruning << "'>"
            "<tr><td>Removed during bubble/superbubble removal<td style='width:50px;background-color:" <<
            Writer::edgeArrowColorRemovedDuringSuperBubbleRemoval << "'>"
            "<tr><td>Removed as low coverage cross edge<td style='width:50px;background-color:" <<
            Writer::edgeArrowColorRemovedAsLowCoverageCrossEdge << "'>"
            "<tr><td>Not removed, opposite strand assembled<td style='width:50px;background-color:" <<
            Writer::edgeArrowColorNotRemovedNotAssembled << "'>"
            "<tr><td>Not removed, assembled<td style='width:50px;background-color:" <<
            Writer::edgeArrowColorNotRemovedAssembled << "'>"
            "</table>";
    }

}



void LocalMarkerGraph::writeColorLegendEdgeLabelsByFlags(ostream& html)
{
    html <<
        "<table>"
        "<tr><td>Removed during transitive reduction<td style='width:50px;background-color:" <<
        Writer::edgeLabelColorRemovedDuringTransitiveReduction << "'>"
        "<tr><td>Removed during pruning<td style='width:50px;background-color:" <<
        Writer::edgeLabelColorRemovedDuringPruning << "'>"
        "<tr><td>Removed during bubble/superbubble removal<td style='width:50px;background-color:" <<
        Writer::edgeLabelColorRemovedDuringSuperBubbleRemoval << "'>"
        "<tr><td>Removed as low coverage cross edge<td style='width:50px;background-color:" <<
        Writer::edgeLabelColorRemovedAsLowCoverageCrossEdge << "'>"
        "<tr><td>Not removed, opposite strand assembled<td style='width:50px;background-color:" <<
        Writer::edgeLabelColorNotRemovedNotAssembled << "'>"
        "<tr><td>Not removed, assembled<td style='width:50px;background-color:" <<
        Writer::edgeLabelColorNotRemovedAssembled << "'>"
        "</table>";

}



void LocalMarkerGraph::Writer::operator()(std::ostream& s) const
{
    // This turns off the tooltip on the graph and the edges.
    s << "tooltip = \" \";\n";

    if((vertexLabels > 0) or (edgeLabels > 0)) {
        s << "overlap = false;\n";
    }
    if(vertexLabels > 0) {
        s << "node [fontname = \"Courier New\" shape=rectangle];\n";
    } else {
        s << "node [shape=point];\n";
    }
    if(edgeLabels > 0) {
        s << "edge [fontname = \"Courier New\" shape=rectangle];\n";
    }

    if(layoutMethod == "dotLr") {
        s << "layout=dot;\n";
        s << "rankdir=LR;\n";
    } else if(layoutMethod == "dotTb") {
        s << "layout=dot;\n";
        s << "rankdir=TB;\n";
    } else if(layoutMethod == "sfdp") {
        s << "layout=sfdp;\n";
        s << "smoothing=triangle;\n";
    } else {
        throw runtime_error("Invalid layout method " + layoutMethod);
    }
}



void LocalMarkerGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const LocalMarkerGraphVertex& vertex = graph[v];
    const auto coverage = vertex.markerInfos.size();
    const string color = vertexColor(vertex);
    SHASTA_ASSERT(coverage > 0);

    // Begin vertex attributes.
    s << "[";

    // Id, so we can use JavaScript code to manipulate the vertex.
    s << "id=vertex" << vertex.vertexId;

    // Tooltip.
    s << " tooltip=\"";
    s << "Vertex " << vertex.vertexId << ", coverage ";
    s << coverage << ", distance " << vertex.distance;
    s << ", Ctrl-click to recenter graph here, right click for detail\"";



    if(vertexLabels == 0) {



        // Vertex area is proportional to coverage.
        s << " width=\"";
        const auto oldPrecision = s.precision(4);
        s << vertexScalingFactor * 0.05 * sqrt(double(coverage));
        s.precision(oldPrecision);
        s << "\"";

        // Color.
        s << " fillcolor=\"" << color << "\" color=\"" << color << "\"";

    } else if(vertexLabels == 1) {



        // Vertex terse label.

        // Color.
        if(vertexColoring != "none") {
            s << " style=filled";
            s << " fillcolor=\"" << color << "\"";
        }

        // Label.
        s << " label=\"";
        s << "Vertex " << vertex.vertexId << "\\n";
        s << "Coverage " << coverage << "\\n";
        s << "Distance " << vertex.distance << "\\n";

        // Marker sequence (run-length).
        const size_t k = graph.k;
        const KmerId kmerId = graph.getKmerId(v);
        const Kmer kmer(kmerId, k);
        kmer.write(s, k);
        s << "\\n";

        // Consensus, if available.
        if(vertex.storedConsensusRepeatCounts.size() == k) {
            // Consensus repeat counts.
            for(size_t i=0; i<k; i++) {
                s << int(vertex.storedConsensusRepeatCounts[i]);
            }
            s << "\\n";

            // Consensus sequence (raw).
            for(size_t i=0; i<k; i++) {
                const Base base = kmer[i];
                const int repeatCount = int(vertex.storedConsensusRepeatCounts[i]);
                for(int l=0; l<repeatCount; l++) {
                    s << base;
                }
            }
            s << "\\n";
        }

        // End the label.
        s << "\"";

    } else {

        // Vertex verbose label.

        // Color.
        if(vertexColoring != "none") {
            s << " style=filled";
            s << " fillcolor=\"" << color << "\"";
        }

        // Label.
        s << " label=<<table border='0' cellborder='0' cellspacing='2'>";

        // Vertex id.
        s << "<tr><td align='left'>Vertex</td><td>" << vertex.vertexId << "</td></tr>";

        // Coverage.
        s << "<tr><td align='left'>Coverage</td><td>" << coverage << "</td></tr><hr/>";

        // Marker information.
        for(const auto& markerInfo: vertex.markerInfos) {
            s << "<tr><td align='left'>" << markerInfo.orientedReadId <<
                "</td><td>" << markerInfo.ordinal << "</td></tr>";
        }
        s << "<hr/>";

        // Marker sequence (run-length).
        const size_t k = graph.k;
        const KmerId kmerId = graph.getKmerId(v);
        const Kmer kmer(kmerId, k);
        s << "<tr><td align='left'>Marker</td><td>";
        kmer.write(s, k);
        s << "</td></tr>";

        // Consensus, if available.
        if(vertex.storedConsensusRepeatCounts.size() == k) {
            // Consensus repeat counts.
            s << "<tr><td align='left'>Consensus counts</td><td>";
            for(size_t i=0; i<k; i++) {
                s << int(vertex.storedConsensusRepeatCounts[i]);
            }
            s << "</td></tr>";

            // Consensus sequence (raw).
            s << "<tr><td align='left'>Consensus sequence</td><td>";
            for(size_t i=0; i<k; i++) {
                const Base base = kmer[i];
                const int repeatCount = int(vertex.storedConsensusRepeatCounts[i]);
                for(int l=0; l<repeatCount; l++) {
                    s << base;
                }
            }
            s << "</td></tr>";
        }

        s<< "</table>>";
    }



    // End vertex attributes.
    s << "]";
}



void LocalMarkerGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{

    const LocalMarkerGraphEdge& edge = graph[e];
    const size_t coverage = edge.coverage();
    const string arrowColor = edgeArrowColor(edge);
    const string labelColor = edgeLabelColor(edge);
    SHASTA_ASSERT(coverage > 0);

    // Begin edge attributes.
    s << "[";

    // Id, so we can use JavaScript code to manipulate the edge.
    s << "id=edge" << edge.edgeId;

    // Tooltip.
    const string tooltipText =
        "Edge " + to_string(edge.edgeId) +
        ", coverage " + to_string(coverage) +
        ", Ctrl-click to recenter graph here, right click for detail";
    s << " tooltip=\"" << tooltipText << "\"";
    s << " labeltooltip=\"" << tooltipText << "\"";

    // Hack to convince graphviz to not ignore the labeltooltip.
    // It seems that this is nolonger necessary, and having it makes it
    // impossible to select text in the edge label.
    // s << " URL=\"#a\"";

    if(edgeThickness == "constant") {
        const double thickness = edgeThicknessScalingFactor;
        s << " penwidth=\"";
        const auto oldPrecision = s.precision(4);
        s <<  thickness;
        s.precision(oldPrecision);
        s << "\"";
    } else {
        // Thickness and weight are determined by coverage.
        double thickness = edgeThicknessScalingFactor;
        if(highlightedOrientedReads.empty()) {
            thickness *= 0.2 * double(max(coverage, size_t(1)));
        }
        s << " penwidth=\"";
        const auto oldPrecision = s.precision(4);
        s <<  thickness;
        s.precision(oldPrecision);
        s << "\" weight=" << coverage;
    }

    // Arrow size.
    s << " arrowsize=\"" << arrowScalingFactor << "\"";

    // Secondary edges in assembly mode 1 are drawn dashed.
    if(edge.isSecondary) {
        s << " style=dashed";
    }

    // Color.
    if(highlightedOrientedReads.empty()) {
        s << " fillcolor=\"" << arrowColor << "\"";
        s << " color=\"" << arrowColor << "\"";
    } else {

        // Highlighting of oriented reads overrides the other color options.

        // Gather the oriented read ids.
        vector<OrientedReadId> orientedReadIds;
        for(const auto& info: edge.infos) {
            const auto& intervals = info.second;
            for(const auto& interval: intervals) {
                orientedReadIds.push_back(interval.orientedReadId);
            }
        }
        sort(orientedReadIds.begin(), orientedReadIds.end());

        vector<string> colors;
        for(const OrientedReadId orientedReadId: orientedReadIds) {
            auto it = highlightedOrientedReads.find(orientedReadId);
            if(it == highlightedOrientedReads.end()) {
                colors.push_back("black");
            } else {
                const double H = it->second;
                colors.push_back(to_string(H) + " " + to_string(S) + " " + to_string(V));
            }
        }


        if(colors.empty()) {
            s << " color=\"black\"";
        } else {

            s << "color=\"";
            for(uint64_t i=0; i<colors.size(); i++) {
                if(i > 0) {
                    s << ":";
                }
                s << colors[i];
            }
            s << "\"";
        }
    }



    // If the edge was not marked as a DAG edge during approximate topological sort,
    // tell graphviz not to use it in constraint assignment.
    // This results in better graph layouts when using dot,
    // because back-edges tend to be low coverage edges.
    if((layoutMethod=="dotLr" or layoutMethod=="dotTb") and (not edge.isDagEdge)) {
        s << " constraint=false";
    }

    // Label.
    if(edgeLabels > 0) {

        s << " label=<<font color=\"black\">";
        s << "<table";
        s << " color=\"black\"";
        s << " bgcolor=\"" << labelColor << "\"";
        s << " border=\"1\"";
        s << " cellborder=\"0\"";
        s << " cellspacing=\"2\"";
        s << ">";

        // Edge id.
        SHASTA_ASSERT(edge.edgeId != MarkerGraph::invalidEdgeId);
        s << "<tr><td align='left'>Edge</td><td>" << edge.edgeId << "</td></tr>";

        // Assembly graph locations.
        for(const auto& p: edge.assemblyGraphLocations) {
            const AssemblyGraph::EdgeId edgeId = p.first;
            const uint32_t position = p.second;
            s << "<tr><td align='left'>Assembly</td><td>" << edgeId << "-" <<
                position << "</td></tr>";
        }

        // Coverage.
        s << "<tr><td align='left'>Coverage</td><td>" << coverage << "</td></tr>";



        // Verbose labels include the detail of all oriented read ids on this edge.
        if(edgeLabels == 2) {

            vector< pair<OrientedReadId, LocalMarkerGraphEdge::Sequence> > table;
            for(const auto& info: edge.infos) {
                const auto& sequence = info.first;
                const auto& intervals = info.second;
                for(const auto& interval: intervals) {
                    table.push_back(make_pair(interval.orientedReadId, sequence));
                }
            }
            sort(table.begin(), table.end(),
                OrderPairsByFirstOnly<OrientedReadId, LocalMarkerGraphEdge::Sequence>());

            s << "<hr/>";
            for(const auto& p: table) {
                s << "<tr><td align='left'>" << p.first << "</td><td>";
                const auto& sequence = p.second;
                if(sequence.overlappingBaseCount > 0) {
                    s << int(sequence.overlappingBaseCount);
                } else {
                    for(const auto base: sequence.sequence) {
                        s << base;
                    }
                }
                s << "</td></tr>";
            }
            s << "<hr/>";
        }



        // Consensus.
        if(edge.consensusOverlappingBaseCount > 0) {
            s << "<tr><td align='left'>Consensus overlap</td><td>";
            s << int(edge.consensusOverlappingBaseCount) << "</td></tr>";
        } else {
            // Consensus sequence (run-length).
            s << "<tr><td align='left'>Consensus (RLE)</td><td>";
            for(const Base base: edge.consensusSequence) {
                s << base;
            }
            s << "</td></tr>";
            // Consensus repeat counts.
            s << "<tr><td align='left'>Consensus counts</td><td>";
            for(const int repeatCount: edge.consensusRepeatCounts) {
                s << repeatCount;
            }
            s << "</td></tr>";
            // Consensus sequence (raw).
            s << "<tr><td align='left'>Consensus (raw)</td><td>";
            for(size_t i=0; i<edge.consensusSequence.size(); i++) {
                const Base base = edge.consensusSequence[i];
                const int repeatCount = edge.consensusRepeatCounts[i];
                for(int l=0; l<repeatCount; l++) {
                    s << base;
                }
            }
            s << "</td></tr>";
        }

        // End the label.
        s << "</table></font>> decorate=true";

    }

    // End edge attributes.
    s << "]";

}
#endif
