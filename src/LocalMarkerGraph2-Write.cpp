// Shasta.
#include "LocalMarkerGraph2.hpp"
#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard libraries.
#include "fstream.hpp"



// Write the graph in Graphviz format.
void LocalMarkerGraph2::write(
    const string& fileName,
    size_t minCoverage,
    int maxDistance,
    bool detailed,
    bool showVertexId) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, minCoverage, maxDistance, detailed, showVertexId);
}
void LocalMarkerGraph2::write(
    ostream& s,
    size_t minCoverage,
    int maxDistance,
    bool detailed,
    bool showVertexId) const
{
    Writer writer(*this, minCoverage, maxDistance, detailed, showVertexId);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&LocalMarkerGraph2Vertex::vertexId, *this));
}

LocalMarkerGraph2::Writer::Writer(
    const LocalMarkerGraph2& graph,
    size_t minCoverage,
    int maxDistance,
    bool detailed,
    bool showVertexId) :
    graph(graph),
    minCoverage(minCoverage),
    maxDistance(maxDistance),
    detailed(detailed),
    showVertexId(showVertexId)
{
}



void LocalMarkerGraph2::Writer::operator()(std::ostream& s) const
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



void LocalMarkerGraph2::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const LocalMarkerGraph2Vertex& vertex = graph[v];
    const auto coverage = vertex.markerInfos.size();
    CZI_ASSERT(coverage > 0);


    // For compact output, the node shape is already defaulted to point,
    // and we don't write a label. The tooltip contains the vertex id,
    // which can be used to create a local subgraph to be looked at
    // in detailed format (use scripts/CreateLocalSubgraph.py).
    if(!detailed) {

        // Compact output.

        // Begin vertex attributes.
        s << "[";

        // Id, so we can use JavaScript code to manipulate the vertex.
        s << "id=vertex" << vertex.vertexId;

        // Tooltip.
        s << " tooltip=\"";
        if(showVertexId) {
            s << "Vertex " << vertex.vertexId << ", coverage ";
        } else {
            s << "Coverage ";
        }
        s << coverage << ", distance " << vertex.distance << ", rank " << vertex.rank;
        s << ", click to recenter graph here, right click for detail\"";

        // Vertex size.
        s << " width=\"";
        const auto oldPrecision = s.precision(4);
        s << 0.05 * sqrt(double(coverage));
        s.precision(oldPrecision);
        s << "\"";

        // Color.
        string color;
        if(vertex.distance == maxDistance) {
            color = "cyan";
        } else if(vertex.distance == 0) {
            color = "#90ee90";
        } else  if(coverage >= minCoverage) {
            color = "black";
        } else {
            color = "red";
        }
        s << " fillcolor=\"" << color << "\" color=\"" << color << "\"";

        // End vertex attributes.
        s << "]";

    } else {

        // Detailed output.
        const size_t k = graph.k;
        const KmerId kmerId = graph.getKmerId(v);
        const Kmer kmer(kmerId, k);

        // Begin vertex attributes.
        s << "[";

        // Color.
        string color;
        if(vertex.distance == maxDistance) {
            color = "cyan";
        } else if(vertex.distance == 0) {
            color = "#90ee90";
        } else if(coverage >= minCoverage) {
            color = "green";
        } else {
            color = "red";
        }
        s << " style=filled";
        s << " fillcolor=\"" << color << "\"";

        // Id, so we can use JavaScript code to manipulate the vertex.
        s << " id=vertex" << vertex.vertexId;

        // Tooltip.
        s << " tooltip=\"";
        if(showVertexId) {
            s << "Vertex " << vertex.vertexId << ", coverage ";
        } else {
            s << "Coverage ";
        }
        s << coverage << ", distance " << vertex.distance << ", rank "  << vertex.rank << "\"";

        // Write the label using Graphviz html-like functionality.
        s << " label=<<font><table border=\"0\">";
        const int columnCount = graph.useRunLengthReads ? 4 : 3;

        // Vertex id.
        if(showVertexId) {
            s << "<tr><td colspan=\"" << columnCount << "\"><b>";
            s << "Vertex " << vertex.vertexId;
            s << "</b></td></tr>";
        }

        // Kmer.
        s << "<tr><td colspan=\"" << columnCount << "\"><b>";
        kmer.write(s, k);
        s << "</b></td></tr>";

        // Coverage.
        s << "<tr><td colspan=\"" << columnCount << "\"><b>";
        s << "Coverage " << coverage;
        s << "</b></td></tr>";

        // Distance.
        s << "<tr><td colspan=\"" << columnCount << "\" ";
        s << " href=\"\"";  // Necessary to activate tooltip.
        s << " id=\"vertexDistance" << vertex.vertexId << "\" tooltip=\"Click to recenter graph here\">";
        s << "<font color=\"blue\"><b><u>Distance " << vertex.distance;
        s << "</u></b></font></td></tr>";

        // Rank.
        s << "<tr><td colspan=\"" << columnCount << "\">";
        s << "<b>Rank " << vertex.rank;
        s << "</b></td></tr>";

        // Column headers.
        s << "<tr><td><b>Read</b></td><td><b>Ord</b></td><td><b>Pos</b></td>";
        if(graph.useRunLengthReads) {
            s << "<td><b>Repeat</b></td>";
        }
        s << "</tr>";

        // A row for each marker of this vertex.
        for(const auto& markerInfo: vertex.markerInfos) {
            const CompressedMarker& marker = graph.markers.begin()[markerInfo.markerId];

            // OrientedReadId
            s << "<tr><td align=\"right\"";
            s << " href=\"exploreRead?readId&amp;" << markerInfo.orientedReadId.getReadId();
            s << "&amp;strand=" << markerInfo.orientedReadId.getStrand() << "\"";
            s << "><font color=\"blue\"><b><u>" << markerInfo.orientedReadId << "</u></b></font></td>";

            // Ordinal.
            s << "<td align=\"right\"";
            s << " href=\"exploreRead?readId=" << markerInfo.orientedReadId.getReadId();
            s << "&amp;strand=" << markerInfo.orientedReadId.getStrand();
            s << "&amp;highlightMarker=" << markerInfo.ordinal;
            s << "\"";
            s << "><font color=\"blue\"><b><u>" << markerInfo.ordinal << "</u></b></font></td>";

            // Position.
            s << "<td align=\"right\"><b>" << marker.position << "</b></td>";

            // Repeat counts.
            if(graph.useRunLengthReads) {
                const vector<uint8_t> counts = graph.getRepeatCounts(markerInfo);
                s << "<td><b>";
                for(size_t i=0; i<k; i++) {
                    if(counts[i] < 10) {
                        s << int(counts[i]);
                    } else {
                        s << "*";
                    }
                }
                s << "</b></td>";
            }

            s << "</tr>";
        }



        // Repeat count consensus.
        if(graph.useRunLengthReads) {
            s << "<tr><td colspan=\"3\" align=\"left\"><b>Repeat consensus</b></td>";
            s << "<td><b>";
            for(size_t position=0; position<graph.k; position++) {
                const size_t bestRepeatCount =
                    vertex.coverages[position].bestRepeatCount(kmer[position]);
                if(bestRepeatCount < 10) {
                    s << bestRepeatCount;
                } else {
                    s << "*";
                }
            }
            s << "</b></td></tr>";

            // Coverage for each repeat count at each position.
            const std::set<size_t> repeatCounts =
                Coverage::findRepeatCounts(vertex.coverages);
            for(const size_t repeatCount: repeatCounts) {
                s << "<tr>";
                s << "<td colspan=\"3\" align=\"left\"><b>Coverage for repeat ";
                s << repeatCount << "</b></td>";
                s << "<td><b>";
                for(size_t position=0; position<graph.k; position++) {
                    const Coverage& consensusInfo = vertex.coverages[position];
                    const Base bestBase = Base(consensusInfo.bestBase());
                    s << vertex.coverages[position].coverageCharacter(bestBase, repeatCount);
                }
                s << "</b></td></tr>";
            }

            // Coverage for the best repeat count at each position.
            s << "<tr><td colspan=\"3\" align=\"left\"><b>Consensus coverage</b></td>";
            s << "<td><b>";
            for(size_t position=0; position<graph.k; position++) {
                const Coverage& consensusInfo = vertex.coverages[position];
                const Base bestBase = Base(consensusInfo.bestBase());
                const size_t bestRepeatCount = consensusInfo.bestBaseBestRepeatCount();
                s << vertex.coverages[position].coverageCharacter(bestBase, bestRepeatCount);
            }
            s << "</b></td></tr>";

            // The raw sequence, based on the best repeat counts.
            s << "<tr><td colspan=\"3\" align=\"left\"><b>Raw consensus</b></td>";
            s << "<td align=\"left\"><b>";
            for(size_t position=0; position<graph.k; position++) {
                const Coverage& coverage = vertex.coverages[position];
                const Base bestBase = Base(coverage.bestBase());
                const size_t bestRepeatCount = coverage.bestBaseBestRepeatCount();
                for(size_t k=0; k<bestRepeatCount; k++) {
                    s << bestBase;
            }
            }
            s << "</b></td></tr>";
        }



        // End the table.
        s << "</table></font>>";

        // End vertex attributes.
        s << "]";
    }
}



void LocalMarkerGraph2::Writer::operator()(std::ostream& s, edge_descriptor e) const
{

    const LocalMarkerGraph2Edge& edge = graph[e];
    const size_t coverage = edge.coverage();
    CZI_ASSERT(coverage > 0);
    const size_t consensus = edge.consensus();

    if(!detailed) {

        // Compact output.

        // Begin edge attributes.
        s << "[";

        // Tooltip.
        s << "tooltip=\"Coverage " << coverage << ", consensus " << consensus << "\"";

        // Color.
        string color;
        if(edge.isSpanningTreeEdge) {
            color = "violet";
        } else if(coverage >= minCoverage) {
            color = "black";
        } else {
            color = "red";
        }
        s << " fillcolor=\"" << color << "\"";
        s << " color=\"" << color << "\"";

        // Thickness is determined by coverage.
        const double thickness = 0.2 * double(coverage==0 ? 1 : coverage);
        s << " penwidth=";
        const auto oldPrecision = s.precision(4);
        s <<  thickness;
        s.precision(oldPrecision);

        // Style.
        if(edge.isSpanningTreeEdge && !edge.isSpanningTreeBestPathEdge) {
            s << " style=dashed";
        }

        // Weight;
        s << " weight=" << coverage;

        // End edge attributes.
        s << "]";

    } else {

        // Detailed output.

        // Begin edge attributes.
        s << "[";

        const string tooltipText = "Coverage " + to_string(coverage) + ", consensus " +to_string(consensus);
        s << " tooltip=\"" << tooltipText << "\"";
        s << " labeltooltip=\"" << tooltipText << "\"";
        // s << " URL=\"#abcdef\"";   // Hack to convince graphviz to not ignore the labeltooltip.

        // Thickness is determined by coverage.
        const double thickness = 0.5 * double(coverage==0 ? 1 : coverage);
        s << " penwidth=";
        const auto oldPrecision = s.precision(4);
        s <<  thickness;
        s.precision(oldPrecision);

        // Style.
        if(edge.isSpanningTreeEdge && !edge.isSpanningTreeBestPathEdge) {
            s << " style=dashed";
        }

        // Color.
        string color;
        if(edge.isSpanningTreeEdge) {
            color = "violet";
        } else if(coverage >= minCoverage) {
            color = "black";
        } else {
            color = "red";
        }
        s << " fillcolor=\"" << color << "\"";
        s << " color=\"" << color << "\"";

        // Label color (used below).
        string labelColor;
        if(color == "black") {
            labelColor = "green";
        } else {
            labelColor = color;
        }


        // Weight;
        s << " weight=" << coverage;

        // If the edge was not marked as a DAG edge during approximate topological sort,
        // tell graphviz not to use it in constraint assignment.
        if(!edge.isDagEdge) {
            s << " constraint=false";
        }

        // Label.
        s << " label=<<font color=\"black\">";
        s << "<table";
        s << " color=\"black\"";
        s << " bgcolor=\"" << labelColor << "\"";
        s << " border=\"0\"";
        s << " cellborder=\"1\"";
        s << " cellspacing=\"1\"";
        s << ">";

        // Consensus and coverage.
        const int columnCount = graph.useRunLengthReads ? 5 : 4;
        s << "<tr><td colspan=\"" << columnCount << "\"><b>Coverage " << coverage << "</b></td></tr>";
        s << "<tr><td colspan=\"" << columnCount << "\"><b>Consensus " << consensus << "</b></td></tr>";

        // Header row.
        s <<
            "<tr>"
            "<td align=\"center\"><b>Read</b></td>"
            "<td align=\"center\"><b>Ord0</b></td>"
            "<td align=\"center\"><b>Ord1</b></td>"
            "<td align=\"center\"><b>Seq</b></td>";
        if(graph.useRunLengthReads) {
            s << "<td align=\"center\"><b>Repeat</b></td>";
        }
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
                if(graph.useRunLengthReads && !info.repeatCounts.empty()) {
                    s << "<td align=\"center\"><b>";
                    for(const uint8_t repeatCount: info.repeatCounts) {
                        if(repeatCount < 10) {
                            s << int(repeatCount);
                        } else {
                            s << "*";
                        }
                    }
                    s << "</b></td>";
                }

                s << "</tr>";
            }
        }



        // If the SeqAn alignment was computed, also write it to the table.
        if(edge.seqanAlignmentWasComputed) {
            s << "<tr><td colspan=\"" << columnCount << "\"><b>SeqAn alignment</b></td></tr>";

            // Add one row to the table for each read.
            for(size_t i=0; i<edge.alignmentInfos.size(); i++) {
                const auto& alignmentInfo = edge.alignmentInfos[i];

                // Begin a new row of the table.
                s << "<tr>";

                // Read id and ordinals.
                s << "<td align=\"right\"";
                s << " href=\"exploreRead?readId&amp;" << alignmentInfo.orientedReadId.getReadId();
                s << "&amp;strand=" << alignmentInfo.orientedReadId.getStrand() << "\"";
                s << "><font color=\"blue\"><b><u>" << alignmentInfo.orientedReadId << "</u></b></font></td>";

                s << "<td align=\"right\"";
                s << " href=\"exploreRead?readId&amp;" << alignmentInfo.orientedReadId.getReadId();
                s << "&amp;strand=" << alignmentInfo.orientedReadId.getStrand();
                s << "&amp;highlightMarker=" << alignmentInfo.ordinals[0];
                s << "&amp;highlightMarker=" << alignmentInfo.ordinals[1];
                s << "\"";
                s << "><font color=\"blue\"><b><u>" << alignmentInfo.ordinals[0] << "</u></b></font></td>";

                s << "<td align=\"right\"";
                s << " href=\"exploreRead?readId&amp;" << alignmentInfo.orientedReadId.getReadId();
                s << "&amp;strand=" << alignmentInfo.orientedReadId.getStrand();
                s << "&amp;highlightMarker=" << alignmentInfo.ordinals[0];
                s << "&amp;highlightMarker=" << alignmentInfo.ordinals[1];
                s << "\"";
                s << "><font color=\"blue\"><b><u>" << alignmentInfo.ordinals[1] << "</u></b></font></td>";

                // SeqAn alignment.
                const seqan::Gaps< seqan::String<seqan::Dna> >& alignmentRow =
                    seqan::row(edge.seqanAlignment, i);
                s << "<td><b>" << alignmentRow << "</b></td>";

                // Repeat counts on SeqAn alignment.
                s << "<td><b>";
                const auto n = seqan::length(alignmentRow);
                size_t position = 0;
                for(size_t j=0; j<n; j++) {
                    if(seqan::isGap(alignmentRow, j)) {
                        s << "-";
                    } else {
                        const int repeatCount = edge.alignmentInfos[i].repeatCounts[position++];
                        if(repeatCount < 10) {
                            s << repeatCount;
                        } else {
                            s << "*";
                        }
                    }
                }
                s << "</b></td>";

                // End this row of the table.
                s << "</tr>";
            }



            // Seqan consensus (run-length sequence).
            s << "<tr><td colspan=\"3\" align=\"left\"><b>Consensus (run-length)</b></td>";
            s << "<td><b>";
            for(const Coverage& coverage: edge.coverages) {
                s << coverage.bestBase();
            }
            s << "</b></td>";
            s << "<td><b>";
            for(const Coverage& coverage: edge.coverages) {
                if(coverage.bestBase().isGap()) {
                    s << "-";
                } else {
                    if(coverage.bestBaseBestRepeatCount() < 10) {
                        s << coverage.bestBaseBestRepeatCount();
                    } else {
                        s << "*";
                    }
                }
            }
            s << "</b></td>";
            s << "</tr>";



            // Consensus coverage for each base.
            for(uint8_t b=0; b<=4; b++) {
                const AlignedBase base = AlignedBase::fromInteger(b);
                s << "<tr><td colspan=\"3\" align=\"left\"><b>Coverage for " << base;
                s << "</b></td>";
                s << "<td><b>";
                for(const Coverage& coverage: edge.coverages) {
                    s << coverage.coverageCharacter(base);
                }
                s << "</b></td>";
                s << "</tr>";
            }

            // Consensus coverage for the best base.
            s << "<tr><td colspan=\"3\" align=\"left\"><b>Consensus coverage</b></td>";
            s << "<td><b>";
            for(const Coverage& coverage: edge.coverages) {
                s << coverage.bestBaseCoverageCharacter();
            }
            s << "</b></td>";
            s << "</tr>";



            // Find the repeat counts that have non-zero coverage on the best base
            // at any position.
            const std::set<size_t> repeatCounts = Coverage::findRepeatCounts(edge.coverages);

            // Coverage for the best base at each position, broken down
            // by repeat count.
            for(size_t repeatCount: repeatCounts) {
                s << "<tr><td colspan=\"4\" align=\"left\"><b>Coverage for repeat count ";
                s << repeatCount << "</b></td>";
                s << "<td><b>";
                for(const Coverage& coverage: edge.coverages) {
                    const AlignedBase bestBase = coverage.bestBase();
                    if(bestBase.isGap()) {
                        s << "-";
                        continue;
                    }
                    s << coverage.coverageCharacter(Base(bestBase), repeatCount);
                }
                s << "</b></td>";
                s << "</tr>";
            }
            s << "<tr><td colspan=\"4\" align=\"left\"><b>Coverage for best repeat count</b></td>";
            s << "<td><b>";
            for(const Coverage& coverage: edge.coverages) {
                if(coverage.bestBase().isGap()) {
                    s << "-";
                    continue;
                }
                s << coverage.coverageCharacter(
                    Base(coverage.bestBase()), coverage.bestBaseBestRepeatCount());
            }
            s << "</b></td>";
            s << "</tr>";



            // Seqan consensus (raw sequence) and its coverage.
            s << "<tr><td colspan=\"3\" align=\"left\"><b>Consensus (raw)</b></td>";
            s << "<td colspan=\"2\"><b>";
            for(const Coverage& coverage: edge.coverages) {
                if(!coverage.bestBase().isGap()) {
                    for(size_t k=0; k<coverage.bestBaseBestRepeatCount(); k++) {
                        s << coverage.bestBase();
                    }
                }
            }
            s << "</b></td>";
            s << "</tr>";
            s << "<tr><td colspan=\"3\" align=\"left\"><b>Consensus (raw) coverage</b></td>";
            s << "<td colspan=\"2\"><b>";
            for(const Coverage& coverage: edge.coverages) {
                const AlignedBase bestBase = coverage.bestBase();
                if(!bestBase.isGap()) {
                    const size_t bestRepeatCount = coverage.bestBaseBestRepeatCount();
                    const char coverageCharacter = coverage.coverageCharacter(Base(bestBase), bestRepeatCount);
                    for(size_t k=0; k<bestRepeatCount; k++) {
                        s << coverageCharacter;
                    }
                }
            }
            s << "</b></td>";
            s << "</tr>";
        }


        // End the label.
        s << "</table></font>> decorate=true";


        // End edge attributes.
        s << "]";
    }

}
