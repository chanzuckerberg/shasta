// shasta
#include "LocalMarkerGraph2.hpp"
#include "CZI_ASSERT.hpp"
#include "findMarkerId.hpp"
#include "LongBaseSequence.hpp"
#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "orderPairs.hpp"
#include "ReadId.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard libraries.
#include "fstream.hpp"
#include "stdexcept.hpp"
#include "tuple.hpp"
#include "utility.hpp"


LocalMarkerGraph2::LocalMarkerGraph2(
    uint32_t k,
    LongBaseSequences& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MemoryMapped::Vector<CompressedGlobalMarkerGraphVertexId>& globalMarkerGraphVertex
    ) :
    k(k),
    reads(reads),
    markers(markers),
    globalMarkerGraphVertex(globalMarkerGraphVertex)
{

}


// Find out if a vertex with the given GlobalMarkerGraphVertexId exists.
// If it exists, return make_pair(true, v).
// Otherwise, return make_pair(false, null_vertex());
std::pair<bool, LocalMarkerGraph2::vertex_descriptor>
    LocalMarkerGraph2::findVertex(GlobalMarkerGraphVertexId vertexId) const
{
    const auto it = vertexMap.find(vertexId);
    if(it == vertexMap.end()) {
        return make_pair(false, null_vertex());
    } else {
        const vertex_descriptor v = it->second;
        return make_pair(true, v);
    }
}


// Add a vertex with the given GlobalMarkerGraphVertexId
// and return its vertex descriptor.
// A vertex with this GlobalMarkerGraphVertexId must not exist.
LocalMarkerGraph2::vertex_descriptor
    LocalMarkerGraph2::addVertex(
    GlobalMarkerGraphVertexId vertexId,
    int distance,
    MemoryAsContainer<MarkerId> vertexMarkers)
{
    // Check that the vertex does not already exist.
    CZI_ASSERT(vertexMap.find(vertexId) == vertexMap.end());

    // Add the vertex and store it in the vertex map.
    const vertex_descriptor v = add_vertex(LocalMarkerGraph2Vertex(vertexId, distance), *this);
    vertexMap.insert(make_pair(vertexId, v));

    // Fill in the marker information for this vertex.
    LocalMarkerGraph2Vertex& vertex = (*this)[v];
    vertex.markerInfos.reserve(vertexMarkers.size());
    for(const MarkerId markerId: vertexMarkers) {
        LocalMarkerGraph2Vertex::MarkerInfo markerInfo;
        markerInfo.markerId = markerId;
        tie(markerInfo.orientedReadId, markerInfo.ordinal) =
            findMarkerId(markerId, markers);
        vertex.markerInfos.push_back(markerInfo);
    }

    return v;
}



// Get the KmerId for a vertex.
KmerId LocalMarkerGraph2::getKmerId(vertex_descriptor v) const
{
    const LocalMarkerGraph2Vertex& vertex = (*this)[v];
    CZI_ASSERT(!vertex.markerInfos.empty());
    const MarkerId firstMarkerId = vertex.markerInfos.front().markerId;
    const CompressedMarker& firstMarker = markers.begin()[firstMarkerId];
    const KmerId kmerId = firstMarker.kmerId;

    // Sanity check that all markers have the same kmerId.
    // At some point this can be removed.
    for(const auto& markerInfo: vertex.markerInfos){
        const CompressedMarker& marker = markers.begin()[markerInfo.markerId];
        CZI_ASSERT(marker.kmerId == kmerId);
    }

    return kmerId;
}



void LocalMarkerGraph2::storeEdgeInfo(edge_descriptor e)
{
    // Access the graph, the edge, and its vertices.
    LocalMarkerGraph2& graph = *this;
    LocalMarkerGraph2Edge& edge = graph[e];
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);
    const LocalMarkerGraph2Vertex& vertex0 = graph[v0];
    const LocalMarkerGraph2Vertex& vertex1 = graph[v1];

    // Map to store the oriented read ids and ordinals, grouped by sequence.
    std::map<LocalMarkerGraph2Edge::Sequence, vector<LocalMarkerGraph2Edge::Info> > sequenceTable;



    // Find pairs of markers for the same oriented read in the two vertices.
    // We exploit the fact that the markers in each
    // of the vertices are sorted.
    auto it0 = vertex0.markerInfos.begin();
    auto it1 = vertex1.markerInfos.begin();
    const auto end0 = vertex0.markerInfos.end();
    const auto end1 = vertex1.markerInfos.end();
    while(it0!=end0 && it1!=end1) {
        const OrientedReadId orientedReadId0 = it0->orientedReadId;
        const OrientedReadId orientedReadId1 = it1->orientedReadId;
        if(orientedReadId0 < orientedReadId1) {
            ++it0;
            continue;
        }
        if(orientedReadId1 < orientedReadId0) {
            ++it1;
            continue;
        }

        // If getting here, the two oriented read ids are the same.
        CZI_ASSERT(orientedReadId0 == orientedReadId1);
        const OrientedReadId orientedReadId = orientedReadId0;


        // Find the streaks of markers for the same oriented readId.
        auto it0StreakEnd = it0;
        while(it0StreakEnd!=end0 && it0StreakEnd->orientedReadId == orientedReadId) {
            ++it0StreakEnd;
        }
        auto it1StreakEnd = it1;
        while(it1StreakEnd!=end1 && it1StreakEnd->orientedReadId == orientedReadId) {
            ++it1StreakEnd;
        }


        // Only do it if both streaks contain one marker,
        // the ordinal for the source vertex
        // is less than the ordinal for the target vertex,
        // and there are no intervening markers that also belong to a
        // vertex of the marker graph.
        if(it0StreakEnd-it0==1 && it1StreakEnd-it1==1 && it0->ordinal<it1->ordinal) {
            const MarkerId markerId0 = it0->markerId;
            const MarkerId markerId1 = it1->markerId;

            // Check that there are no intervening markers that also belong to a
            // vertex of the marker graph.
            bool interveningVertexFound = false;
            for(MarkerId markerId=markerId0+1; markerId!=markerId1; markerId++) {
                if(globalMarkerGraphVertex[markerId] != invalidCompressedGlobalMarkerGraphVertexId) {
                    interveningVertexFound = true;
                    break;
                }

            }
            if(!interveningVertexFound) {

                const CompressedMarker& marker0 = markers.begin()[markerId0];
                const CompressedMarker& marker1 = markers.begin()[markerId1];

                // Fill in the sequence information.
                LocalMarkerGraph2Edge::Sequence sequence;
                if(marker1.position <= marker0.position + k) {
                    sequence.overlappingBaseCount = uint8_t(marker0.position + k - marker1.position);
                } else {
                    sequence.overlappingBaseCount = 0;
                    const auto read = reads[orientedReadId.getReadId()];
                    const uint32_t readLength = uint32_t(read.baseCount);
                    for(uint32_t position=marker0.position+k;  position!=marker1.position; position++) {
                        shasta::Base base;
                        if(orientedReadId.getStrand() == 0) {
                            base = read.get(position);
                        } else {
                            base = read.get(readLength - 1 - position);
                            base.complementInPlace();
                        }
                        sequence.sequence.push_back(base);
                    }

                }

                // Store it.
                sequenceTable[sequence].push_back(LocalMarkerGraph2Edge::Info(
                    orientedReadId,
                    it0->ordinal,
                    it1->ordinal));
            }
        }

        // Update the iterators to point to the end of the streaks.
        it0 = it0StreakEnd;
        it1 = it1StreakEnd;

#if 0
        if(it1->orientedReadId == orientedReadId) {

            // Fill in the sequence information.
            LocalMarkerGraph2Edge::Sequence sequence;
            const CompressedMarker& marker0 = markers.begin()[markerId0];
            const CompressedMarker& marker1 = markers.begin()[markerId1];

            if(marker1.position <= marker0.position + k) {
                sequence.overlappingBaseCount = uint8_t(marker0.position + k - marker1.position);
            } else {
                sequence.overlappingBaseCount = 0;
                const auto read = reads[orientedReadId.getReadId()];
                const uint32_t readLength = uint32_t(read.baseCount);
                for(uint32_t position=marker0.position+k;  position!=marker1.position; position++) {
                    shasta::Base base;
                    if(orientedReadId.getStrand() == 0) {
                        base = read.get(position);
                    } else {
                        base = read.get(readLength - 1 - position);
                        base.complementInPlace();
                    }
                    sequence.sequence.push_back(base);
                }

            }

            sequenceTable[sequence].push_back(LocalMarkerGraph2Edge::Info(orientedReadId, it0->ordinal));
        }
        ++it0;
        ++it1;
#endif
    }

    // Copy to the edge infos data structure.
    edge.infos.clear();
    copy(sequenceTable.begin(), sequenceTable.end(), back_inserter(edge.infos));

    // Sort by decreasing size of the infos vector.
    sort(edge.infos.begin(), edge.infos.end(),
        OrderPairsBySizeOfSecondGreater<
        LocalMarkerGraph2Edge::Sequence,
        vector<LocalMarkerGraph2Edge::Info> >());
}



// Write the graph in Graphviz format.
void LocalMarkerGraph2::write(
    const string& fileName,
    size_t minCoverage,
    size_t minConsensus,
    int maxDistance,
    bool detailed) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, minCoverage, minConsensus, maxDistance, detailed);
}
void LocalMarkerGraph2::write(
    ostream& s,
    size_t minCoverage,
    size_t minConsensus,
    int maxDistance,
    bool detailed) const
{
    Writer writer(*this, minCoverage, minConsensus, maxDistance, detailed);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&LocalMarkerGraph2Vertex::vertexId, *this));
}

LocalMarkerGraph2::Writer::Writer(
    const LocalMarkerGraph2& graph,
    size_t minCoverage,
    size_t minConsensus,
    int maxDistance,
    bool detailed) :
    graph(graph),
    minCoverage(minCoverage),
    minConsensus(minConsensus),
    maxDistance(maxDistance),
    detailed(detailed)
{
}



void LocalMarkerGraph2::Writer::operator()(std::ostream& s) const
{
    // This turns off the tooltip on the graph.
    s << "tooltip = \" \";\n";

    if(detailed) {
        s << "layout=dot;\n";
        s << "ratio=expand;\n";
        s << "node [fontname = \"Courier New\" shape=rectangle];\n";
        s << "edge [fontname = \"Courier New\"];\n";
    } else {
        s << "layout=sfdp;\n";
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
        s << " tooltip=\"Coverage " << coverage << ", distance " << vertex.distance;
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
            color = "lightGreen";
        } else  if(coverage >= minCoverage) {
            color = "black";
        } else if(coverage == 1) {
            color = "#ff000080";  // Red, half way transparent
        } else if(coverage == 2) {
            color = "#ff800080";  // Orange, half way transparent
        } else {
            color = "#ff80ff80";  // Purple, half way transparent
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
            color = "#7cfc00";  // Bright green
        } else if(coverage >= minCoverage) {
            color = "green";
        } else if(coverage == 1) {
            color = "#ff0000";  // Red
        } else if(coverage == 2) {
            color = "#ff8000";  // Orange
        } else {
            color = "#ff80ff";  // Purple
        }
        s << " style=filled";
        s << " fillcolor=\"" << color << "\"";

        // Id, so we can use JavaScript code to manipulate the vertex.
        s << " id=vertex" << vertex.vertexId;

        // Tooltip.
        s << " tooltip=\"Coverage " << coverage << ", distance " << vertex.distance << "\"";

        // Write the label using Graphviz html-like functionality.
        s << " label=<<font><table border=\"0\">";

        // Kmer.
        s << "<tr><td colspan=\"3\"><b>";
        kmer.write(s, k);
        s << "</b></td></tr>";

        // Coverage.
        s << "<tr><td colspan=\"3\"><b>";
        s << "Coverage " << coverage;
        s << "</b></td></tr>";

        // Distance.
        s << "<tr><td colspan=\"3\" ";
        s << " href=\"\"";  // Necessary to activate tooltip.
        s << " id=\"vertexDistance" << vertex.vertexId << "\" tooltip=\"Click to recenter graph here\">";
        s << "<font color=\"blue\"><b><u>Distance " << vertex.distance;
        s << "</u></b></font></td></tr>";

        // Column headers.
        s << "<tr><td><b>Read</b></td><td><b>Ord</b></td><td><b>Pos</b></td></tr>";

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
            s << "#" << markerInfo.ordinal << "\"";
            s << "><font color=\"blue\"><b><u>" << markerInfo.ordinal << "</u></b></font></td>";

            // Position.
            s << "<td align=\"right\"><b>" << marker.position << "</b></td></tr>";
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
    const size_t consensus = edge.consensus();

    if(!detailed) {

        // Compact output.

        // Begin edge attributes.
        s << "[";

        s << "tooltip=\"Consensus " << consensus << ", coverage " << coverage << "\"";

        // Color is determined by consensus.
        string color;
        if(consensus >= minCoverage) {
            color = "black";
        } else if(consensus == 1) {
            color = "#ff000080";  // Red, half way transparent
        } else if(consensus == 2) {
            color = "#ff800080";  // Orange, half way transparent
        } else {
            color = "#ff80ff80";  // Purple, half way transparent
        }
        s << " fillcolor=\"" << color << "\"";
        s << " color=\"" << color << "\"";

        // Thickness is determined by coverage.
        s << " penwidth=";
        const auto oldPrecision = s.precision(4);
        s <<  0.5 * double(coverage);
        s.precision(oldPrecision);

        // End edge attributes.
        s << "]";

    } else {

        // Detailed output.

        // Begin edge attributes.
        s << "[";

        // const string tooltipText = "Consensus " + to_string(consensus) + ", coverage " +to_string(coverage);
        s << " tooltip=\" \"";
        s << " labeltooltip=\" \"";
        // s << " URL=\"#abcdef\"";   // Hack to convince graphviz to not ignore the labeltooltip.

        // Thickness is determined by coverage.
        s << " penwidth=";
        const auto oldPrecision = s.precision(4);
        s <<  0.5 * double(coverage);
        s.precision(oldPrecision);

        // Color is determined by consensus.
        string color;
        string fillColor;
        if(consensus >= minConsensus) {
            color = "black";
            fillColor = "green";
        } else if(consensus == 1) {
            color = "red";
            fillColor = color;
        } else if(consensus == 2) {
            color = "#ff8000";  // Orange
            fillColor = color;
        } else {
            color = "#ff80ff";  // Purple
            fillColor = color;
        }
        s << " fillcolor=\"" << fillColor << "\"";
        s << " color=\"" << color << "\"";

        // Label.
        s << " label=<<font color=\"black\">";
        s << "<table";
        s << " color=\"black\"";
        s << " bgcolor=\"" << fillColor << "\"";
        s << " border=\"0\"";
        s << " cellborder=\"1\"";
        s << " cellspacing=\"0\"";
        s << ">";

        // Consensus and coverage.
        s << "<tr><td colspan=\"4\"><b>Consensus " << consensus << "</b></td></tr>";
        s << "<tr><td colspan=\"4\"><b>Coverage " << coverage << "</b></td></tr>";

        // Header row.
        s <<
            "<tr>"
            "<td align=\"center\"><b>Read</b></td>"
            "<td align=\"center\"><b>Ord0</b></td>"
            "<td align=\"center\"><b>Ord1</b></td>"
            "<td align=\"center\"><b>Seq</b></td>"
            "</tr>";

        // Loop over the infos table for this edge.
        for(const auto& p: edge.infos) {
            const auto& sequence = p.first;
            const auto& infos = p.second;

            // Construct the string representing this sequence.
            string sequenceString;
            if(sequence.sequence.empty()) {
                sequenceString = to_string(sequence.overlappingBaseCount);
            } else {
                for(const shasta::Base base: sequence.sequence) {
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
                s << "#" << info.ordinals[1] << "\"";
                s << "><font color=\"blue\"><b><u>" << info.ordinals[0] << "</u></b></font></td>";

                s << "<td align=\"right\"";
                s << " href=\"exploreRead?readId&amp;" << info.orientedReadId.getReadId();
                s << "&amp;strand=" << info.orientedReadId.getStrand();
                s << "&amp;highlightMarker=" << info.ordinals[0];
                s << "&amp;highlightMarker=" << info.ordinals[1];
                s << "#" << info.ordinals[1] << "\"";
                s << "><font color=\"blue\"><b><u>" << info.ordinals[1] << "</u></b></font></td>";

                s << "<td align=\"center\"><b>";
                if(it == infos.begin()) {
                    s << sequenceString;
                } else {
                    s << "=";
                }
                s << "</b></td></tr>";
            }
        }

        s << "</table></font>> decorate=true";


        // End edge attributes.
        s << "]";
    }

}

