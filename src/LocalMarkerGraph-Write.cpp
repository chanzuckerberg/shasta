#ifndef SHASTA_STATIC_EXECUTABLE

// Shasta.
#include "LocalMarkerGraph.hpp"
#include "ConsensusCaller.hpp"
#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard libraries.
#include "fstream.hpp"



// Write the graph in Graphviz format.
void LocalMarkerGraph::write(
    const string& fileName,
    int maxDistance,
    bool detailed) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, maxDistance, detailed);
}
void LocalMarkerGraph::write(
    ostream& s,
    int maxDistance,
    bool detailed) const
{
    Writer writer(*this, maxDistance, detailed);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&LocalMarkerGraphVertex::vertexId, *this));
}

LocalMarkerGraph::Writer::Writer(
    const LocalMarkerGraph& graph,
    int maxDistance,
    bool detailed) :
    graph(graph),
    maxDistance(maxDistance),
    detailed(detailed)
{
}



// Vertex and edge colors.
const string LocalMarkerGraph::Writer::vertexColorZeroDistance                          = "#6666ff";
const string LocalMarkerGraph::Writer::vertexColorIntermediateDistance                  = "#00ccff";
const string LocalMarkerGraph::Writer::vertexColorMaxDistance                           = "#66ffff";
const string LocalMarkerGraph::Writer::edgeArrowColorRemovedDuringTransitiveReduction   = "#ff0000";
const string LocalMarkerGraph::Writer::edgeArrowColorRemovedDuringPruning               = "#ff00ff";
const string LocalMarkerGraph::Writer::edgeArrowColorRemovedDuringSuperBubbleRemoval    = "#009900";
const string LocalMarkerGraph::Writer::edgeArrowColorNotRemovedNotAssembled             = "#663300";
const string LocalMarkerGraph::Writer::edgeArrowColorNotRemovedAssembled                = "#000000";
const string LocalMarkerGraph::Writer::edgeLabelColorRemovedDuringTransitiveReduction   = "#ff9999";
const string LocalMarkerGraph::Writer::edgeLabelColorRemovedDuringPruning               = "#c03280";
const string LocalMarkerGraph::Writer::edgeLabelColorRemovedDuringSuperBubbleRemoval    = "#99ff99";
const string LocalMarkerGraph::Writer::edgeLabelColorNotRemovedNotAssembled             = "#996600";
const string LocalMarkerGraph::Writer::edgeLabelColorNotRemovedAssembled                = "#999999";
const string& LocalMarkerGraph::Writer::vertexColor(const LocalMarkerGraphVertex& vertex) const
{
    if(vertex.distance == 0) {
        return vertexColorZeroDistance;
    } else if(vertex.distance == maxDistance) {
        return vertexColorMaxDistance;
    } else {
        return vertexColorIntermediateDistance;
    }
}
const string& LocalMarkerGraph::Writer::edgeArrowColor(const LocalMarkerGraphEdge& edge) const
{
    if(edge.wasRemovedByTransitiveReduction) {
        return edgeArrowColorRemovedDuringTransitiveReduction;
    } else if(edge.wasPruned) {
        return edgeArrowColorRemovedDuringPruning;
    } else if (edge.isSuperBubbleEdge) {
        return edgeArrowColorRemovedDuringSuperBubbleRemoval;
    } else {
        if(edge.wasAssembled) {
            return edgeArrowColorNotRemovedAssembled;
        } else {
            return edgeArrowColorNotRemovedNotAssembled;
        }
    }
}
const string& LocalMarkerGraph::Writer::edgeLabelColor(const LocalMarkerGraphEdge& edge) const
{
    if(edge.wasRemovedByTransitiveReduction) {
        return edgeLabelColorRemovedDuringTransitiveReduction;
    } else if(edge.wasPruned) {
        return edgeLabelColorRemovedDuringPruning;
    } else if (edge.isSuperBubbleEdge) {
        return edgeLabelColorRemovedDuringSuperBubbleRemoval;
    } else {
        if(edge.wasAssembled) {
            return edgeLabelColorNotRemovedAssembled;
        } else {
            return edgeLabelColorNotRemovedNotAssembled;
        }
    }
}



void LocalMarkerGraph::writeColorLegend(ostream& html)
{
    html <<
        "<table style='font-size:10px'>"
        "<tr><th class=centered colspan=3>Marker graph color legend"
        "<tr><td rowspan=4>Vertices"
        "<tr><td>Zero distance<td style='width:50px;background-color:" <<
        Writer::vertexColorZeroDistance << "'>"
        "<tr><td>Intermediate distances<td style='width:50px;background-color:" <<
        Writer::vertexColorIntermediateDistance << "'>"
        "<tr><td>Maximum distance<td style='width:50px;background-color:" <<
        Writer::vertexColorMaxDistance << "'>"
        "<tr><td rowspan=6>Edge<br>arrows"
        "<tr><td>Removed during transitive reduction<td style='width:50px;background-color:" <<
        Writer::edgeArrowColorRemovedDuringTransitiveReduction << "'>"
        "<tr><td>Removed during pruning<td style='width:50px;background-color:" <<
        Writer::edgeArrowColorRemovedDuringPruning << "'>"
        "<tr><td>Removed during bubble/superbubble removal<td style='width:50px;background-color:" <<
        Writer::edgeArrowColorRemovedDuringSuperBubbleRemoval << "'>"
        "<tr><td>Not removed, not assembled<td style='width:50px;background-color:" <<
        Writer::edgeArrowColorNotRemovedNotAssembled << "'>"
        "<tr><td>Not removed, assembled<td style='width:50px;background-color:" <<
        Writer::edgeArrowColorNotRemovedAssembled << "'>"
        "<tr><td rowspan=6>Edge<br>labels"
        "<tr><td>Removed during transitive reduction<td style='width:50px;background-color:" <<
        Writer::edgeLabelColorRemovedDuringTransitiveReduction << "'>"
        "<tr><td>Removed during pruning<td style='width:50px;background-color:" <<
        Writer::edgeLabelColorRemovedDuringPruning << "'>"
        "<tr><td>Removed during bubble/superbubble removal<td style='width:50px;background-color:" <<
        Writer::edgeLabelColorRemovedDuringSuperBubbleRemoval << "'>"
        "<tr><td>Not removed, not assembled<td style='width:50px;background-color:" <<
        Writer::edgeLabelColorNotRemovedNotAssembled << "'>"
        "<tr><td>Not removed, assembled<td style='width:50px;background-color:" <<
        Writer::edgeLabelColorNotRemovedAssembled << "'>"
       "</table>";
}



void LocalMarkerGraph::Writer::operator()(std::ostream& s) const
{
    // This turns off the tooltip on the graph and the edges.
    s << "tooltip = \" \";\n";

    if(detailed) {
        s << "layout=dot;\n";
        s << "rankdir=LR;\n";
        s << "ratio=expand;\n";
        s << "node [fontname = \"Courier New\" shape=rectangle];\n";
        s << "edge [fontname = \"Courier New\"];\n";
    } else {
        s << "layout=sfdp;\n";
        s << "smoothing=triangle;\n";
        s << "ratio=expand;\n";
        s << "node [shape=point];\n";
    }
}



void LocalMarkerGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const LocalMarkerGraphVertex& vertex = graph[v];
    const auto coverage = vertex.markerInfos.size();
    const string& color = vertexColor(vertex);
    CZI_ASSERT(coverage > 0);

    // Begin vertex attributes.
    s << "[";

    // Id, so we can use JavaScript code to manipulate the vertex.
    s << "id=vertex" << vertex.vertexId;

    // Tooltip.
    s << " tooltip=\"";
    s << "Vertex " << vertex.vertexId << ", coverage ";
    s << coverage << ", distance " << vertex.distance;
    s << ", click to recenter graph here, right click for detail\"";



    if(!detailed) {

        // Compact output: point, no label.

        // Vertex area is proportional to coverage.
        s << " width=\"";
        const auto oldPrecision = s.precision(4);
        s << 0.05 * sqrt(double(coverage));
        s.precision(oldPrecision);
        s << "\"";

        // Color.
        s << " fillcolor=\"" << color << "\" color=\"" << color << "\"";


    } else {

        // Detailed output.

        // Color.
        s << " style=filled";
        s << " fillcolor=\"" << color << "\"";

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

        // End the label.
        s << "\"";
    }

    // End vertex attributes.
    s << "]";
}



void LocalMarkerGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{

    const LocalMarkerGraphEdge& edge = graph[e];
    const size_t coverage = edge.coverage();
    const string& arrowColor = edgeArrowColor(edge);
    const string& labelColor = edgeLabelColor(edge);
    CZI_ASSERT(coverage > 0);

    if(!detailed) {

        // Compact output.

        // Begin edge attributes.
        s << "[";

        // Tooltip.
        s << "tooltip=\"Edge " << edge.edgeId << ", coverage " << coverage << "\"";

        // Color.
        s << " fillcolor=\"" << arrowColor << "\"";
        s << " color=\"" << arrowColor << "\"";

        // Id, so we can use JavaScript code to manipulate the edge.
        s << "id=edge" << edge.edgeId;

        // Thickness is determined by coverage.
        const double thickness = 0.2 * double(coverage==0 ? 1 : coverage);
        s << " penwidth=";
        const auto oldPrecision = s.precision(4);
        s <<  thickness;
        s.precision(oldPrecision);

        // Weight;
        s << " weight=" << coverage;

        // If the edge was not marked as a DAG edge during approximate topological sort,
        // tell graphviz not to use it in constraint assignment.
        if(!edge.isDagEdge) {
            s << " constraint=false";
        }

        // End edge attributes.
        s << "]";

    } else {

        // Detailed output.

        // Begin edge attributes.
        s << "[";

        const string tooltipText = "Edge " + to_string(edge.edgeId) + ", coverage " + to_string(coverage);
        s << " tooltip=\"" << tooltipText << "\"";
        s << " labeltooltip=\"" << tooltipText << "\"";
        // s << " URL=\"#abcdef\"";   // Hack to convince graphviz to not ignore the labeltooltip.

        // Thickness is determined by coverage.
        const double thickness = 0.5 * double(coverage==0 ? 1 : coverage);
        s << " penwidth=";
        const auto oldPrecision = s.precision(4);
        s <<  thickness;
        s.precision(oldPrecision);

        // Color.
        s << " fillcolor=\"" << arrowColor << "\"";
        s << " color=\"" << arrowColor << "\"";

        // Id, so we can use JavaScript code to manipulate the edge.
        s << "id=edge" << edge.edgeId;

        // Weight;
        s << " weight=" << coverage;

        // Label.
        s << " label=<<font color=\"black\">";
        s << "<table";
        s << " color=\"black\"";
        s << " bgcolor=\"" << labelColor << "\"";
        s << " border=\"0\"";
        s << " cellborder=\"1\"";
        s << " cellspacing=\"1\"";
        s << ">";

        // Edge id.
        const int columnCount = 5;
        CZI_ASSERT(edge.edgeId != MarkerGraph::invalidEdgeId);
        s << "<tr><td colspan=\"" << columnCount << "\"><b>Edge " << edge.edgeId << "</b></td></tr>";

        // Assembly vertex id.
        if((edge.assemblyEdgeId != std::numeric_limits<AssemblyGraph::VertexId>::max())) {
            s << "<tr><td colspan=\"" << columnCount << "\"><b>Position " << edge.positionInAssemblyEdge <<
                " in assembly graph edge " << edge.assemblyEdgeId << "</b></td></tr>";
        }

        // Coverage.
        s << "<tr><td colspan=\"" << columnCount << "\"><b>Coverage " << coverage << "</b></td></tr>";

        // Header row.
        s <<
            "<tr>"
            "<td align=\"center\"><b>Read</b></td>"
            "<td align=\"center\"><b>Ord0</b></td>"
            "<td align=\"center\"><b>Ord1</b></td>"
            "<td align=\"center\"><b>Seq</b></td>";
        s << "<td align=\"center\"><b>Repeat</b></td>";
        s << "</tr>";

        // Loop over the infos table for this edge.
        for(const auto& p: edge.infos) {
            const auto& sequence = p.first;
            const auto& infos = p.second;

            // Construct the string representing this sequence.
            string sequenceString;
            if(sequence.sequence.empty()) {
                sequenceString = to_string(sequence.overlappingBaseCount);
            } else {
                for(const Base base: sequence.sequence) {
                    sequenceString.push_back(base.character());
                }
            }



            for(auto it=infos.begin(); it!=infos.end(); ++it) {
                const auto& info = *it;
                s << "<tr><td align=\"right\"";
                s << " href=\"exploreRead?readId&amp;" << info.orientedReadId.getReadId();
                s << "&amp;strand=" << info.orientedReadId.getStrand() << "\"";
                s << "><font color=\"blue\"><b><u>" << info.orientedReadId << "</u></b></font></td>";

                s << "<td align=\"right\"";
                s << " href=\"exploreRead?readId&amp;" << info.orientedReadId.getReadId();
                s << "&amp;strand=" << info.orientedReadId.getStrand();
                s << "&amp;highlightMarker=" << info.ordinals[0];
                s << "&amp;highlightMarker=" << info.ordinals[1];
                s << "\"";
                s << "><font color=\"blue\"><b><u>" << info.ordinals[0] << "</u></b></font></td>";

                s << "<td align=\"right\"";
                s << " href=\"exploreRead?readId&amp;" << info.orientedReadId.getReadId();
                s << "&amp;strand=" << info.orientedReadId.getStrand();
                s << "&amp;highlightMarker=" << info.ordinals[0];
                s << "&amp;highlightMarker=" << info.ordinals[1];
                s << "\"";
                s << "><font color=\"blue\"><b><u>" << info.ordinals[1] << "</u></b></font></td>";

                s << "<td align=\"center\"><b>";
                if(it == infos.begin()) {
                    if(sequenceString.size() > 100) {
                        s << "Too long";
                    } else {
                        s << sequenceString;
                    }
                } else {
                    s << "=";
                }
                s << "</b></td>";

                // Write out the repeat counts, if necessary.
                if(!info.repeatCounts.empty()) {
                    s << "<td align=\"center\"><b>";
                    if(sequenceString.size() > 100) {
                        s << "Too long";
                    } else {
                        for(const uint8_t repeatCount: info.repeatCounts) {
                            if(repeatCount < 10) {
                                s << int(repeatCount);
                            } else {
                                s << "*";
                            }
                        }
                    }
                    s << "</b></td>";
                }

                s << "</tr>";
            }
        }



        // End the label.
        s << "</table></font>> decorate=true";

        // If the edge was not marked as a DAG edge during approximate topological sort,
        // tell graphviz not to use it in constraint assignment.
        if(!edge.isDagEdge) {
            s << " constraint=false";
        }

        // End edge attributes.
        s << "]";
    }

}
#endif
