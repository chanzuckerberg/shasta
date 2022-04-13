#ifdef SHASTA_HTTP_SERVER

#include "Assembler.hpp"
#include "mode3.hpp"
#include "mode3-LocalAssemblyGraph.hpp"
using namespace shasta;
using namespace mode3;

#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/right_open_interval.hpp>


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
        num_edges(localAssemblyGraph) << " links."
        "<p>";



    // Display the local assembly graph.
    html << "<div style='display: inline-block; vertical-align:top'>";
    localAssemblyGraph.writeSvg(html, options);
    html << "</div>";
    addSvgDragAndZoom(html);

    // Side panel.
    html << "<div style='display: inline-block'>";



    // Highlight a segment.
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
        Highlight segment
        <input id=highlightInputField type=text onchange="highlightSegment()" size=10>
        )stringDelimiter";



    // Zoom to a segment.
    html << R"stringDelimiter(
        <script>
        function zoomToSegment()
        {
            // Get the segment id from the input field.
            inputField = document.getElementById("zoomInputField");
            segmentId = inputField.value;
            inputField.value = "";

            // Find the bounding box and its center.
            var element = document.getElementById("Segment-" + segmentId);
            var box = element.getBBox();
            var xCenter = box.x + 0.5 * box.width;
            var yCenter = box.y + 0.5 * box.height;

            // Change the viewbox of the svg to be a bit larger than a square
            // containing the bounding box.
            var enlargeFactor = 5.;
            var size = enlargeFactor * Math.max(box.width, box.height);
            width = size;
            height = size;
            x = xCenter - 0.5 * size;
            y = yCenter - 0.5 * size;
            var svg = document.querySelector('svg');
            svg.setAttribute('viewBox', `${x} ${y} ${size} ${size}`);
            ratio = size / svg.getBoundingClientRect().width;

        }
        </script>
        <p>Zoom to segment
        <input id=zoomInputField type=text onchange="zoomToSegment()" size=10>
        )stringDelimiter";



    // Table that will be automatically updated when the mouse is on a segment.
    html << R"zzz(  
<table style='font-size:9'>
<tr><th class='left' style='width:16em'>Segment id<td id='segmentIdCell' class=centered style='width:8em'>
<tr><th class='left' style='width:16em'>Distance from start segment<td id='distanceCell' class=centered style='width:8em'>
<tr><th class='left' style='width:16m'>Path length<td id='pathLengthCell' class=centered style='width:8em'>
<tr><th class='left' style='width:16em'>Average edge coverage<td id='coverageCell' class=centered style='width:8em'>
<tr><th class='left' style='width:16em'>Number of oriented reads on this segment<td id='orientedReadsCell' class=centered style='width:8em'>
<tr><th class='left style='width:16em''>Number of oriented reads on this segment that are also in the reference segment
<td id='comonOrientedReadsCell' class=centered style='width:8em'>
<tr><th class='left' style='width:16em'>Number of oriented reads in this segment only or in the reference segment only
which are too short to appear on both.
<td id='tooShortCell' class=centered style='width:8em'>
<tr><th class='left' style='width:16em'>Number of oriented reads in this segment but not in the reference segment
<td id='missingFromReferenceSegmentCell' class=centered style='width:8em'>
<tr><th class='left' style='width:16em'>Number of oriented reads in the reference segment but not in this segment
<td id='missingFromDisplayedSegmentCell' class=centered style='width:8em'>
</table>
<script>
function onMouseEnterSegment(id, distance, pathLength, coverage, orientedReads,
    common, tooShort, missingFromReference, missingFromDisplayed)
{
    document.getElementById('segmentIdCell').innerHTML = id;
    document.getElementById('distanceCell').innerHTML = distance;
    document.getElementById('pathLengthCell').innerHTML = pathLength;
    document.getElementById('coverageCell').innerHTML = coverage;
    document.getElementById('orientedReadsCell').innerHTML = orientedReads;
    document.getElementById('comonOrientedReadsCell').innerHTML = common;
    document.getElementById('tooShortCell').innerHTML = tooShort;
    document.getElementById('missingFromReferenceSegmentCell').innerHTML = missingFromReference;
    document.getElementById('missingFromDisplayedSegmentCell').innerHTML = missingFromDisplayed;
}
function onMouseExitSegment()
{
    document.getElementById('segmentIdCell').innerHTML = '';
    document.getElementById('distanceCell').innerHTML = '';
    document.getElementById('pathLengthCell').innerHTML = '';
    document.getElementById('coverageCell').innerHTML = '';
    document.getElementById('orientedReadsCell').innerHTML = '';
    document.getElementById('comonOrientedReadsCell').innerHTML = '';
    document.getElementById('tooShortCell').innerHTML = '';
    document.getElementById('missingFromReferenceSegmentCell').innerHTML = '';
    document.getElementById('missingFromDisplayedSegmentCell').innerHTML = '';
}
</script>
    )zzz";



    // Change segment thickness
    html << R"stringDelimiter(
    <p><table>
    <tr><th class=left>Segment thickness<td>
    <button type='button' onClick='segmentThickness(0.5)' style='width:2em'>--</button>
    <button type='button' onClick='segmentThickness(0.8)' style='width:2em'>-</button>
    <button type='button' onClick='segmentThickness(1.25)' style='width:2em'>+</button>
    <button type='button' onClick='segmentThickness(2)' style='width:2em'>++</button>
        <script>
        function segmentThickness(factor)
        {
            const group = document.getElementById('LocalAssemblyGraph-segments');
            for (let i=0; i<group.children.length; i++) {
                path = group.children[i];
                if(path.tagName == 'path') {
                    path.setAttribute('stroke-width', factor * path.getAttribute('stroke-width'));
                }
            }
        }
        </script>
        )stringDelimiter";



    // Change link thickness
    html << R"stringDelimiter(
    <tr><th class=left>Link thickness<td>
    <button type='button' onClick='linkThickness(0.5)' style='width:2em'>--</button>
    <button type='button' onClick='linkThickness(0.8)' style='width:2em'>-</button>
    <button type='button' onClick='linkThickness(1.25)' style='width:2em'>+</button>
    <button type='button' onClick='linkThickness(2.)' style='width:2em'>++</button>
        <script>
        function linkThickness(factor)
        {
            const group1 = document.getElementById('LocalAssemblyGraph-links');
            for (let i=0; i<group1.children.length; i++) {
                group2 = group1.children[i];
                if(group2.tagName == 'g') {
                    for (let j=0; j<group2.children.length; j++) {
                        path = group2.children[j];
                        if(path.tagName == 'path') {
                            path.setAttribute('stroke-width', factor * path.getAttribute('stroke-width'));
                        }
                    }
                }
            }
        }
        </script>
     </table>
        )stringDelimiter";



    // End of side panel.
    html << "</div>";

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

    // Check that we have a valid segmentId.
    if(segmentId >= assemblyGraph3.paths.size()) {
        html << "Invalid segment id. Maximum valid value is " <<
            assemblyGraph3.paths.size() - 1 << ".";
        return;
    }

    // Access the marker graph path for this segment.
    const auto path = assemblyGraph3.paths[segmentId];

    // Get information about the oriented reads of this segment.
    mode3::AssemblyGraph::SegmentOrientedReadInformation orientedReads;
    assemblyGraph3.getOrientedReadsOnSegment(segmentId, orientedReads);

    const auto oldPrecision = html.precision(1);
    const auto oldFlags = html.setf(std::ios_base::fixed, std::ios_base::floatfield);
    html <<
        "<h1>Assembly graph segment " << segmentId << "</h1>"
        "<p><table>"
        "<tr><th class=left>Length of marker graph path<td class=centered>" << path.size() <<
        "<tr><th class=left>Average marker graph edge coverage on path<td class=centered>" <<
        orientedReads.averageCoverage <<
        "<tr><th class=left>Number of distinct oriented reads on path<td class=centered>" << orientedReads.infos.size() <<
        "</table>";
    html.precision(oldPrecision);
    html.flags(oldFlags);



    // Write the oriented reads in a table.
    html <<
        "<h2>Oriented reads on this segment</h2>"
        "<table>"
        "<tr>"
        "<th>Oriented<br>read"
        "<th>Average<br>offset";
    for(const auto& info: orientedReads.infos) {
        html<<
            "<tr>"
            "<td class=centered>" << info.orientedReadId <<
            "<td class=centered>" << info.averageOffset;
    }
    html << "</table>";



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



void Assembler::exploreMode3AssemblyGraphSegmentPair(
    const vector<string>& request,
    ostream& html)
{
    using boost::icl::discrete_interval;
    using boost::icl::intersects;
    using boost::icl::length;

    SHASTA_ASSERT(assemblyGraph3Pointer);
    const mode3::AssemblyGraph& assemblyGraph3 = *assemblyGraph3Pointer;

    // Get the segment ids from the request.
    uint64_t segmentId0;
    const bool segmentId0IsPresent = getParameterValue(request, "segmentId0", segmentId0);
    uint64_t segmentId1;
    const bool segmentId1IsPresent = getParameterValue(request, "segmentId1", segmentId1);



    // Write the form.
    html <<
        "<h3>Display details for a pair assembly graph segment</h3>"
        "<form>"
        "<table>"

        "<tr>"
        "<td>Segment id 0"
        "<td><input type=text required name=segmentId0 size=8 style='text-align:center'"
        " value='" << (segmentId0IsPresent ? to_string(segmentId0) : "") <<
        "'>"

        "<tr>"
        "<td>Segment id 1"
        "<td><input type=text required name=segmentId1 size=8 style='text-align:center'"
        " value='" << (segmentId1IsPresent ? to_string(segmentId1) : "") <<
        "'>"

        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";

    // If the segmentId's were not specified, stop here.
    if(not segmentId0IsPresent) {
        return;
    }
    if(not segmentId1IsPresent) {
        return;
    }

    // Check that we have valid segmentId's.
    if(segmentId0 >= assemblyGraph3.paths.size()) {
        html << "Invalid segment id. Maximum valid value is " <<
            assemblyGraph3.paths.size() - 1 << ".";
        return;
    }
    if(segmentId1 >= assemblyGraph3.paths.size()) {
        html << "Invalid segment id. Maximum valid value is " <<
            assemblyGraph3.paths.size() - 1 << ".";
        return;
    }


    // Get information about the oriented reads of these segments.
    mode3::AssemblyGraph::SegmentOrientedReadInformation orientedReads0;
    mode3::AssemblyGraph::SegmentOrientedReadInformation orientedReads1;
    assemblyGraph3.getOrientedReadsOnSegment(segmentId0, orientedReads0);
    assemblyGraph3.getOrientedReadsOnSegment(segmentId1, orientedReads1);
    const uint64_t length0 = assemblyGraph3.paths.size(segmentId0);
    const uint64_t length1 = assemblyGraph3.paths.size(segmentId1);

    // Estimate the offset between the segments and count missing
    // oriented reads.
    mode3::AssemblyGraph::SegmentPairInformation segmentPairInformation;
    assemblyGraph3.analyzeSegmentPair(
            segmentId0, segmentId1,
            orientedReads0, orientedReads1,
            markers, segmentPairInformation);
    const uint64_t commonCount = segmentPairInformation.commonCount;


    /// Write a table with general information about this pair of segments.
    html <<
        "<p>"
        "See caption after the tables."
        "<p>"
        "<table>"
        "<tr><th class=left>Length of segment " << segmentId0 <<
        "<td class=centered>" << length0 <<
        "<tr><th class=left>Length of segment " << segmentId1 <<
        "<td class=centered>" << length1 <<
        "<tr><th class=left>Number of common oriented reads"
        "<td class=centered>" << commonCount;
    if(commonCount) {
        html <<
            "<tr><th class=left>Estimated offset between segment " << segmentId0 <<
            " and segment " << segmentId1 <<
            "<td class=centered>" << segmentPairInformation.offset <<
            "<tr><th class=left>Number of oriented reads on segment " << segmentId0 <<
            " that are too short to appear on both segments"
            "<td class=centered>" << segmentPairInformation.tooShortCount[0] <<
            "<tr><th class=left>Number of oriented reads on segment " << segmentId1 <<
            " that are too short to appear on both segments"
            "<td class=centered>" << segmentPairInformation.tooShortCount[1] <<
            "<tr><th class=left>Number of oriented reads missing from segment " << segmentId0 <<
            "<td class=centered>" << segmentPairInformation.missingCount[0] <<
            "<tr><th class=left>Number of oriented reads missing from segment " << segmentId1 <<
            "<td class=centered>" << segmentPairInformation.missingCount[1];
    }
    html <<  "</table>";



    // Write a table with a row for each oriented read.
    html <<
        "<p>"
        "<table>"
        "<tr>"
        "<th>Oriented<br>read"
        "<th>Length"
        "<th>Average<br>offset of<br>oriented read<br>relative to<br>segment " << segmentId0 <<
        "<th>Average<br>offset of<br>oriented read<br>relative to<br>segment " << segmentId1 <<
        "<th>Estimated<br>offset of<br>segment " << segmentId1 <<
        "<br>relative to<br>segment " << segmentId0 <<
        "<th>Hypothetical<br>offset of<br>oriented read<br>relative to<br>segment " << segmentId0 <<
        "<th>Hypothetical<br>offset of<br>oriented read<br>relative to<br>segment " << segmentId1 <<
        "<th>Hypothetical<br>overlap of<br>oriented read<br>with<br>segment " << segmentId0 <<
        "<th>Hypothetical<br>overlap of<br>oriented read<br>with<br>segment " << segmentId1 <<
        "<th>On both<br>segments" <<
        "<th>Too<br>short" <<
        "<th>On segment<br>" << segmentId0 << "<br>only,<br>missing from<br>segment<br>" << segmentId1 <<
        "<th>On segment<br>" << segmentId1 << "<br>only,<br>missing from<br>segment<br>" << segmentId0;


    // Set up a joint loop over oriented reads in the two segments.
    const auto begin0 = orientedReads0.infos.begin();
    const auto begin1 = orientedReads1.infos.begin();
    const auto end0 = orientedReads0.infos.end();
    const auto end1 = orientedReads1.infos.end();
    auto it0 = begin0;
    auto it1 = begin1;



    while(true) {

        // At end of both segments.
        if(it0 == end0 and it1 == end1) {
            break;
        }



        // Only on segment 0.
        if((it1 == end1) or ((it0!=end0) and (it0->orientedReadId < it1->orientedReadId))) {
            const int64_t orientedReadLength = markers.size(it0->orientedReadId.getValue());
            html <<
                "<tr>"
                "<td class=centered>" <<
                "<a href='exploreRead?readId=" << it0->orientedReadId.getReadId() <<
                "&strand=" << it0->orientedReadId.getStrand() << "'>" << it0->orientedReadId << "</a>"
                "<td class=centered>" << orientedReadLength <<
                "<td class=centered>" << it0->averageOffset <<
                "<td>"
                "<td><td>";

            if(commonCount) {
                // Compute the hypothetical range of the oriented read relative
                // to the beginning of segment 1.
                const discrete_interval<int64_t> orientedReadRange1(
                    it0->averageOffset - segmentPairInformation.offset,
                    it0->averageOffset - segmentPairInformation.offset + orientedReadLength);
                const discrete_interval<int64_t> segment1Range(0, length1);
                const bool wouldOverlap = intersects(orientedReadRange1, segment1Range);
                html <<
                    "<td class=centered>" << orientedReadRange1.lower() <<
                    "<td><td class=centered>" << length(orientedReadRange1 & segment1Range);
                if(wouldOverlap) {
                    html << "<td><td><td class=centered>&#10003;<td>";
                } else {
                    html << "<td><td class=centered>&#10003;<td><td>";
                }
            } else {
                html << "<td><td><td><td><td><td><td>";
            }
            ++it0;
        }



        // Only on segment 1
        else if((it0 == end0) or ((it1!=end1) and (it1->orientedReadId < it0->orientedReadId))) {
            const int64_t orientedReadLength = markers.size(it1->orientedReadId.getValue());
            html <<
                "<tr>"
                "<td class=centered>" <<
                "<a href='exploreRead?readId=" << it1->orientedReadId.getReadId() <<
                "&strand=" << it1->orientedReadId.getStrand() << "'>" << it1->orientedReadId << "</a>"
                "<td class=centered>" << orientedReadLength <<
                "<td>"
                "<td class=centered>" << it1->averageOffset <<
                "<td>";

            if(commonCount) {
                // Compute the hypothetical range of the oriented read relative
                // to the beginning of segment 0.
                const discrete_interval<int64_t> orientedReadRange0(
                    it1->averageOffset + segmentPairInformation.offset,
                    it1->averageOffset + segmentPairInformation.offset + orientedReadLength);
                const discrete_interval<int64_t> segment0Range(0, length0);
                const bool wouldOverlap = intersects(orientedReadRange0, segment0Range);
                html <<
                    "<td class=centered>" << orientedReadRange0.lower() <<
                    "<td><td class=centered>" << length(orientedReadRange0 & segment0Range) << "<td>";
                if(wouldOverlap) {
                    html << "<td><td><td><td class=centered>&#10003;";
                } else {
                    html << "<td><td class=centered>&#10003;<td><td>";
                }

            } else {
                html << "<td><td><td><td><td><td><td><td>";
            }

            ++it1;
        }

        // On both segments.
        else {
            html <<
                "<tr>"
                "<td class=centered>" <<
                "<a href='exploreRead?readId=" << it0->orientedReadId.getReadId() <<
                "&strand=" << it0->orientedReadId.getStrand() << "'>" << it0->orientedReadId << "</a>"
                "<td class=centered>" << markers.size(it0->orientedReadId.getValue()) <<
                "<td class=centered>" << it0->averageOffset <<
                "<td class=centered>" << it1->averageOffset <<
                "<td class=centered>" << it0->averageOffset - it1->averageOffset <<
                "<td><td><td><td>"
                "<td class=centered>&#10003;<td><td><td>";

            ++it0;
            ++it1;
        }
    }
    html << "</table>";



    // Write a caption for this messy table.
    html <<
        "<p><ul>"
        "<li>All offsets, overlaps, and lengths are in markers."
        "<li>On both segments: this oriented read appears in both segments "
        "and was used to estimate their offset."
        "<li>Too short: this oriented read appears in only one of the two segments "
        "and, based on the estimated offset between the segments, "
        "is too short to appear in the other segment."
        "<li>On one segment only, missing from the other segment: "
        "this oriented read appears in one segment only and, "
        "based on the estimated offset between the segments, was expected to also "
        "appear in the other segment.";



}


#endif
