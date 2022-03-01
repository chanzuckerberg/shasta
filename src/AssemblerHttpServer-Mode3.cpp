#ifdef SHASTA_HTTP_SERVER

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
    mode3::LocalAssemblyGraph::SvgOptions options(request);

    uint64_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    uint64_t startSegmentId;
    const bool startSegmentIdIsPresent = getParameterValue(request, "startSegmentId", startSegmentId);

    double timeout = 30.;
    getParameterValue(request, "timeout", timeout);



    // Write the form.
    html <<
        "<h2>Display the local assembly graph near a given segment</h2>"
        "<form>"
        "<table>"

        "<tr>"
        "<td>Start segment"
        "<td class=centered><input type=text required name=startSegmentId size=8 style='text-align:center'"
        " value='" << (startSegmentIdIsPresent ? to_string(startSegmentId) : "") <<
        "'>"

        "<tr>"
        "<td>Maximum distance in the assembly graph (edges)"
        "<td class=centered><input type=text name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance <<
        "'>"

        "<tr>"
        "<td>Timeout for graph layout (seconds)"
        "<td class=centered><input type=text name=timeout size=8 style='text-align:center'"
        " value='" << timeout <<
        "'>";

    options.addFormRows(html);

    html <<
        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";



    if(not startSegmentIdIsPresent) {
        return;
    }

    if(startSegmentId >= assemblyGraph3Pointer->paths.size()) {
        html << "<p>Invalid start segment id. Maximum valid value is " <<
            assemblyGraph3Pointer->paths.size() - 1;
        return;
    }
    if(options.referenceSegmentId >= assemblyGraph3Pointer->paths.size()) {
        html << "<p>Invalid reference segment id. Maximum valid value is " <<
            assemblyGraph3Pointer->paths.size() - 1;
        return;
    }


    html << "<h1>Local assembly graph near segment " << startSegmentId << "</h1></p>";

    // Create the local assembly graph.
    mode3::LocalAssemblyGraph localAssemblyGraph(
        markerGraph,
        *assemblyGraph3Pointer,
        startSegmentId, maxDistance);
    localAssemblyGraph.computeLayout(options, timeout);
    localAssemblyGraph.computeSegmentTangents();
    html << "<p>The local assembly graph has " <<
        num_vertices(localAssemblyGraph) << " segments and " <<
        num_edges(localAssemblyGraph) << " links.";

    // Allow manually highlighting selected vertices.
    html << R"stringDelimiter(
        <script>
        function highlightSegment()
        {
            // Get the segment id from the input field.
            inputField = document.getElementById("highlightInputField");
            segmentId = inputField.value;
            inputField.value = "";

            // Make it dashed and wider.
            var element = document.getElementById("Segment-" + segmentId);
            var thickness = element.getAttribute("stroke-width");
            element.style.strokeDasharray = 0.2 * thickness;
            element.setAttribute("stroke-width", 2. * thickness);
        }
        </script>
        <p>
        <input id=highlightInputField type=text onchange="highlightSegment()" size=10>
        Enter a segment id to highlight, then press Enter.
        To highlight multiple segments, enter them one at a time in the same way.
        <p>
        )stringDelimiter";

    // Display the local assembly graph.
    localAssemblyGraph.writeSvg(html, options);



    // To facilitate debugging and testing, also write a gfa file
    // that represents the LocalAssemblyGraph.
    localAssemblyGraph.writeGfa("LocalAssemblyGraph.gfa");

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
        "<h3>Display details of an assembly graph segment</h3>"
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

    // Check tha7 we have a valid segmentId.
    if(segmentId >= assemblyGraph3.paths.size()) {
        html << "Invalid segment id. Maximum valid value is " <<
            assemblyGraph3.paths.size() - 1 << ".";
        return;
    }

    // Access the marker graph path for this segment.
    const auto path = assemblyGraph3.paths[segmentId];

    // Get the oriented reads and average edge coverage.
    vector<OrientedReadId> orientedReadIds;
    const double averageEdgeCoverage = assemblyGraph3.findOrientedReadsOnSegment(segmentId, orientedReadIds);

    const auto oldPrecision = html.precision(1);
    const auto oldFlags = html.setf(std::ios_base::fixed, std::ios_base::floatfield);
    html <<
        "<h1>Assembly graph segment " << segmentId << "</h1>"
        "<p><table>"
        "<tr><th class=left>Length of marker graph path<td class=centered>" << path.size() <<
        "<tr><th class=left>Average marker graph edge coverage on path<td class=centered>" << averageEdgeCoverage <<
        "<tr><th class=left>Number of distinct oriented reads on path<td class=centered>" << orientedReadIds.size() <<
        "</table>";
    html.precision(oldPrecision);
    html.flags(oldFlags);



    // Write the path in a table.
    html <<
        "<h2>Marker graph path for this segment</h2>"
        "<table>"
        "<tr>"
        "<th>Position"
        "<th>Edge"
        "<th>Coverage"
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
            "<td class=centered>" <<
            "<a href='exploreMarkerGraphEdge?edgeId=" << edgeId <<
            "'>" << edgeId << "</a>"
            "<td class=centered>" << markerGraph.edgeMarkerIntervals.size(edgeId) <<
            "<td class=centered>" <<
            "<a href='exploreMarkerGraphVertex?vertexId=" << vertexId0 <<
            "'>" << vertexId0 << "</a>"
            "<td class=centered>" <<
            "<a href='exploreMarkerGraphVertex?vertexId=" << vertexId1 <<
            "'>" << vertexId1 << "</a>"
            "\n";



    }
    html << "</table>";

}



void Assembler::exploreMode3AssemblyGraphLink(
    const vector<string>& request,
    ostream& html)
{
    SHASTA_ASSERT(assemblyGraph3Pointer);
    const mode3::AssemblyGraph& assemblyGraph3 = *assemblyGraph3Pointer;

    // Get the link id from the request.
    uint64_t linkId;
    const bool linkIdIsPresent = getParameterValue(request, "linkId", linkId);



    // Write the form.
    html <<
        "<h3>Display details of an assembly graph link</h3>"
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

    const Link& link = assemblyGraph3.links[linkId];
    const auto transitions = assemblyGraph3.transitions[linkId];
    const uint64_t segmentId0 = link.segmentId0;
    const uint64_t segmentId1 = link.segmentId1;
    const auto path0 = assemblyGraph3.paths[segmentId0];
    const auto path1 = assemblyGraph3.paths[segmentId1];
    const uint64_t pathLength0 = path0.size();
    const uint64_t pathLength1 = path1.size();
    const MarkerGraph::VertexId vertexId0 = markerGraph.edges[path0.back().edgeId].target;
    const MarkerGraph::VertexId vertexId1 = markerGraph.edges[path1.front().edgeId].source;

    const double linkSeparation = mode3::linkSeparation(transitions, pathLength0);

    html <<
        "<h1>Assembly graph link " << linkId << "</h1>"
        "<p><table>"
        "<tr><th>Segment<th>Id<th>Path<br>length"
        "<tr><th class = left>Source segment<td class=centered>" << segmentId0 << "<td class=centered>" << pathLength0 <<
        "<tr><th class = left>Target segment<td class=centered>" << segmentId1 << "<td class=centered>" << pathLength1 <<
        "</table>";

    if(vertexId0 == vertexId1) {
        html << "<p>The paths of these segments are consecutive.";
    } else {
        html << "<p>The paths of these segments are not consecutive.";
    }


    const auto oldPrecision = html.precision(1);
    const auto oldFlags = html.setf(std::ios_base::fixed, std::ios_base::floatfield);
    html <<
        "<p><table>"
        "<tr><th class = left tooltip='Number of supporting transitions'>Coverage<td class=centered>" <<
        transitions.size() <<
        "<tr><th class = left>Average link separation<td class=centered>" <<
        linkSeparation <<
        "</table>";
    html.precision(oldPrecision);
    html.flags(oldFlags);


    html <<
        "<h2>Transitions</h2>"
        "<p><table><tr>"
        "<th class=centered>Oriented<br>read<br>id"
        "<th class=centered>Last<br>position<br>on segment<br>" << link.segmentId0 <<
        "<th class=centered>Last<br>ordinal<br>on segment<br>" << link.segmentId0 <<
        "<th class=centered>First<br>position<br>on segment<br>" << link.segmentId1 <<
        "<th class=centered>First<br>ordinal<br>on segment<br>" << link.segmentId1 <<
        "<th class=centered>Link<br>separation";


    for(const auto& p: transitions) {
        const OrientedReadId orientedReadId = p.first;
        const Transition& transition = p.second;
        const auto& pseudoPathEntry0 = transition[0];
        const auto& pseudoPathEntry1 = transition[1];

        SHASTA_ASSERT(pseudoPathEntry1.ordinals[0] >= pseudoPathEntry0.ordinals[1]);

        const int64_t linkSeparation =
            int64_t(pseudoPathEntry1.ordinals[0] - pseudoPathEntry0.ordinals[1]) -
            int64_t(pathLength0 - 1 - pseudoPathEntry0.position) -
            int64_t(pseudoPathEntry1.position);

        html <<
            "<tr><td class=centered>" << orientedReadId <<

            "<td class=centered>" << pseudoPathEntry0.position <<
            "<td class=centered>" << pseudoPathEntry0.ordinals[1] <<

            "<td class=centered>" << pseudoPathEntry1.position <<
            "<td class=centered>" << pseudoPathEntry1.ordinals[0] <<

            "<td class=centered>" << linkSeparation;
    }
    html << "</table>";




}

#endif
