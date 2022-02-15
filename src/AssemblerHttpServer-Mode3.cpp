#include "Assembler.hpp"
#include "mode3.hpp"
#include "mode3-LocalAssemblyGraph.hpp"
using namespace shasta;
using namespace mode3;


void Assembler::exploreMode3AssemblyGraph(
    const vector<string>& request,
    ostream& html)
{
    SHASTA_ASSERT(assemblyGraph3Pointer);

    // Get the parameters for the request.
    uint64_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    uint64_t startSegmentId;
    const bool startSegmentIdIsPresent = getParameterValue(request, "startSegmentId", startSegmentId);

    uint32_t sizePixels = 600;
    getParameterValue(request, "sizePixels", sizePixels);



    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<td>Start segment"
        "<td><input type=text required name=startSegmentId size=8 style='text-align:center'"
        " value='" << (startSegmentIdIsPresent ? to_string(startSegmentId) : "") <<
        "'>"

        "<tr>"
        "<td>Maximum distance"
        "<td><input type=text name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance <<
        "'>"

        "<tr>"
        "<td>Graphics size in pixels"
        "<td><input type=text name=sizePixels size=8 style='text-align:center'"
        " value='" << sizePixels <<
        "'>"

        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";




    if(not startSegmentIdIsPresent) {
        return;
    }


    html << "<h1>Local assembly graph near segment " << startSegmentId << "</h1></p>";


    // Create the local assembly graph and write it to html in svg format.
    mode3::LocalAssemblyGraph localAssemblyGraph(
        *assemblyGraph3Pointer, startSegmentId, maxDistance);
    // localAssemblyGraph.writeGraphviz("LocalAssemblyGraph.dot");
    localAssemblyGraph.writeSvg1(html, sizePixels);

}



void Assembler::exploreMode3AssemblyGraphSegment(
    const vector<string>& request,
    ostream& html)
{
    SHASTA_ASSERT(assemblyGraph3Pointer);
    const mode3::AssemblyGraph& assemblyGraph3 = *assemblyGraph3Pointer;

    // Get the segment id from the request.
    uint64_t segmentId;
    const bool segmentIdIsPresent = getParameterValue(request, "segmentId", segmentId);



    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<td>Segment id"
        "<td><input type=text required name=segmentId size=8 style='text-align:center'"
        " value='" << (segmentIdIsPresent ? to_string(segmentId) : "") <<
        "'>"

        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";

    // If the segmentId was not specified, stop here.
    if(not segmentIdIsPresent) {
        return;
    }

    // Check tha we have a valid segmentId.
    if(segmentId >= assemblyGraph3.paths.size()) {
        html << "Invalid segment id. Maximum valid value is " <<
            assemblyGraph3.paths.size() - 1 << ".";
        return;
    }

    // Access the marker graph path for this segment.
    const auto path = assemblyGraph3.paths[segmentId];

    html << "<h1>Assembly graph segment " << segmentId << "</h1>";
    html << "<p>Path length is " << path.size() << ".";



    // Write the path in a table.
    html <<
        "<h2>Marker graph path for this segment</h2>"
        "<table>"
        "<tr>"
        "<th>Position"
        "<th>Edge"
        "<th>Source<br>vertex"
        "<th>Target<br>vertex";

    for(uint64_t position=0; position<path.size(); position++) {
        const MarkerGraphEdgeInfo& markerGraphEdgeInfo = path[position];

        // For now we are not generating virtual marker graph edges.
        SHASTA_ASSERT(not markerGraphEdgeInfo.isVirtual);

        const MarkerGraph::EdgeId edgeId = markerGraphEdgeInfo.edgeId;
        const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
        const MarkerGraph::VertexId vertexId0 = edge.source;
        const MarkerGraph::VertexId vertexId1 = edge.target;

        html << "<tr>"
            "<td class=centered>" << position <<
            "<td class=centered>" << edgeId <<
            "<td class=centered>" << vertexId0 <<
            "<td class=centered>" << vertexId1 <<
            "\n";



    }
    html << "</table>";

}



void Assembler::exploreMode3AssemblyGraphLink(
    const vector<string>& request,
    ostream& html)
{
    SHASTA_ASSERT(assemblyGraph3Pointer);

    // Get the link id from the request.
    uint64_t linkId;
    const bool linkIdIsPresent = getParameterValue(request, "linkId", linkId);



    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<td>Link id"
        "<td><input type=text required name=linkId size=8 style='text-align:center'"
        " value='" << (linkIdIsPresent ? to_string(linkId) : "") <<
        "'>"

        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";

    // If the segmentId was not specified, stop here.
    if(not linkIdIsPresent) {
        return;
    }


    html << "<h1>Assembly graph link " << linkId << "</h1>";

}
