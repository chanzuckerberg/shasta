// Shasta.
#include "Assembler.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "mode3.hpp"
#include "mode3-AssemblyPath.hpp"
#include "mode3-LocalAssemblyGraph.hpp"
#include "PngImage.hpp"
using namespace shasta;
using namespace mode3;

// Boost library.
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/right_open_interval.hpp>

// Standard library.
#include "fstream.hpp"


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

    if(startSegmentId >= assemblyGraph3Pointer->markerGraphPaths.size()) {
        html << "<p>Invalid start segment id. Maximum valid value is " <<
            assemblyGraph3Pointer->markerGraphPaths.size() - 1;
        return;
    }
    if(options.referenceSegmentId >= assemblyGraph3Pointer->markerGraphPaths.size()) {
        html << "<p>Invalid reference segment id. Maximum valid value is " <<
            assemblyGraph3Pointer->markerGraphPaths.size() - 1;
        return;
    }


    html << "<h1>Local assembly graph near segment " << startSegmentId << "</h1></p>";



    // Create the local assembly graph, or reuse the last one, if possible.
    static shared_ptr<mode3::LocalAssemblyGraph> lastLocalAssemblyGraphPointer;
    static shared_ptr<mode3::LocalAssemblyGraph::SvgOptions> lastOptions;
    static uint64_t lastStartSegmentId = invalid<uint64_t>;
    static uint64_t lastMaxDistance = invalid<uint64_t>;
    const bool canReuse =
        lastLocalAssemblyGraphPointer and
        (startSegmentId == lastStartSegmentId) and
        (maxDistance == lastMaxDistance) and
        options.hasSameLayoutOptions(*lastOptions);
    if(canReuse) {
        cout << "Reusing the previous mode3::LocalAssemblyGraph." << endl;
    } else {
        lastLocalAssemblyGraphPointer = make_shared<mode3::LocalAssemblyGraph>(
            markerGraph,
            *assemblyGraph3Pointer,
            startSegmentId, maxDistance);
        lastOptions = make_shared<mode3::LocalAssemblyGraph::SvgOptions>(options);
        lastStartSegmentId = startSegmentId;
        lastMaxDistance = maxDistance;
        lastLocalAssemblyGraphPointer->computeLayout(options, timeout);
        lastLocalAssemblyGraphPointer->computeSegmentTangents();
    }
    mode3::LocalAssemblyGraph& localAssemblyGraph = *lastLocalAssemblyGraphPointer;

    html << "<p>The local assembly graph has " <<
        num_vertices(localAssemblyGraph) << " segments and " <<
        num_edges(localAssemblyGraph) << " links."
        "<p>";



    // Display the local assembly graph.
    localAssemblyGraph.writeHtml(html, options);

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
    if(segmentId >= assemblyGraph3.markerGraphPaths.size()) {
        html << "Invalid segment id. Maximum valid value is " <<
            assemblyGraph3.markerGraphPaths.size() - 1 << ".";
        return;
    }

    // Access the marker graph path for this segment.
    const auto path = assemblyGraph3.markerGraphPaths[segmentId];

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
        assemblyGraph3.segmentCoverage[segmentId] <<
        "<tr><th class=left>Number of distinct oriented reads on path<td class=centered>" << orientedReads.infos.size();

    // Write the incoming and outgoing links.
    html << "<tr><th class=left>Incoming links<td class=centered>";
    for(const uint64_t linkId: assemblyGraph3.linksByTarget[segmentId]) {
        html << "<a href='exploreMode3AssemblyGraphLink?linkId=" << linkId << "'>" << linkId << "</a> ";
    }
    html << "<tr><th class=left>Outgoing links<td class=centered>";
    for(const uint64_t linkId: assemblyGraph3.linksBySource[segmentId]) {
        html << "<a href='exploreMode3AssemblyGraphLink?linkId=" << linkId << "'>" << linkId << "</a> ";
    }


    html << "</table>";
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
        const MarkerGraphEdgeId& edgeId = path[position];
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

    const mode3::AssemblyGraph::Link& link = assemblyGraph3.links[linkId];
    const auto transitions = assemblyGraph3.transitions[linkId];
    const uint64_t segmentId0 = link.segmentId0;
    const uint64_t segmentId1 = link.segmentId1;
    const auto path0 = assemblyGraph3.markerGraphPaths[segmentId0];
    const auto path1 = assemblyGraph3.markerGraphPaths[segmentId1];
    const uint64_t pathLength0 = path0.size();
    const uint64_t pathLength1 = path1.size();

    html <<
        "<h1>Assembly graph link " << linkId << "</h1>"
        "<p><table>"
        "<tr><th>Segment<th>Id<th>Path<br>length"
        "<tr><th class = left>Source segment<td class=centered>" << segmentId0 << "<td class=centered>" << pathLength0 <<
        "<tr><th class = left>Target segment<td class=centered>" << segmentId1 << "<td class=centered>" << pathLength1 <<
         "</table>";

    if(link.segmentsAreAdjacent) {
        html << "<p>The paths of these segments are adjacent.";
    } else {
        html << "<p>The paths of these segments are not adjacent.";
    }


    const auto oldPrecision = html.precision(1);
    const auto oldFlags = html.setf(std::ios_base::fixed, std::ios_base::floatfield);
    html <<
        "<p><table>"
        "<tr><th class = left tooltip='Number of supporting transitions'>Coverage<td class=centered>" <<
        transitions.size() <<
        "<tr><th class = left>Average link separation<td class=centered>" <<
        link.separation <<
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
    if(segmentId0 >= assemblyGraph3.markerGraphPaths.size()) {
        html << "Invalid segment id. Maximum valid value is " <<
            assemblyGraph3.markerGraphPaths.size() - 1 << ".";
        return;
    }
    if(segmentId1 >= assemblyGraph3.markerGraphPaths.size()) {
        html << "Invalid segment id. Maximum valid value is " <<
            assemblyGraph3.markerGraphPaths.size() - 1 << ".";
        return;
    }


    // Get information about the oriented reads of these segments.
    mode3::AssemblyGraph::SegmentOrientedReadInformation orientedReads0;
    mode3::AssemblyGraph::SegmentOrientedReadInformation orientedReads1;
    assemblyGraph3.getOrientedReadsOnSegment(segmentId0, orientedReads0);
    assemblyGraph3.getOrientedReadsOnSegment(segmentId1, orientedReads1);
    const uint64_t length0 = assemblyGraph3.markerGraphPaths.size(segmentId0);
    const uint64_t length1 = assemblyGraph3.markerGraphPaths.size(segmentId1);

    // Estimate the offset between the segments and count missing
    // oriented reads.
    mode3::AssemblyGraph::SegmentPairInformation segmentPairInformation;
    assemblyGraph3.analyzeSegmentPair(
            segmentId0, segmentId1,
            orientedReads0, orientedReads1,
            markers, segmentPairInformation);
    const uint64_t commonCount = segmentPairInformation.commonCount;



    /// Write a table with information about this pair of segments.
    html <<
        "<p>"
        "<table>"

        "<tr>"
        "<th class=left>Segment id"
        "<td class=centered>" << segmentId0 <<
        "<td class=centered>" << segmentId1 <<

        "<tr title='Segment length in marker graph edges'>"
        "<th class=left>Length"
        "<td class=centered>" << length0 <<
        "<td class=centered>" << length1 <<

        "<tr title='Total number of oriented reads in this segment'>"
        "<th class=left>Total"
        "<td class=centered>" << segmentPairInformation.totalCount[0] <<
        "<td class=centered>" << segmentPairInformation.totalCount[1] <<

        "<tr title='Number of oriented reads present in both segments'>"
        "<th class=left>Common"
        "<td class=centered>" << segmentPairInformation.commonCount <<
        "<td class=centered>" << segmentPairInformation.commonCount;

    if(segmentPairInformation.commonCount > 0) {
        const auto oldPrecision = html.precision(2);
        const auto oldFlags = html.setf(std::ios_base::fixed, std::ios_base::floatfield);
        html <<
            "<tr title='Number of oriented reads in this segment that are too short to appear in the other segment'>"
            "<th class=left>Short"
            "<td class=centered>" << segmentPairInformation.shortCount[0] <<
            "<td class=centered>" << segmentPairInformation.shortCount[1] <<

            "<tr title='Number of oriented reads in this segment that are "
            "unexpectedly missing in the other segment'>"
            "<th class=left>Unexplained"
            "<td class=centered>" << segmentPairInformation.unexplainedCount[0] <<
            "<td class=centered>" << segmentPairInformation.unexplainedCount[1] <<

            "<tr title='Jaccard similarity without counting short reads'>"
            "<th class=left>Jaccard"
            "<td class=centered>" << segmentPairInformation.jaccard() <<
            "<td class=centered>" << segmentPairInformation.jaccard() <<

            "<tr title='Fraction of oriented reads in this segment that are "
            "unexpectedly missing in the other segment'>"
            "<th class=left>Unexplained fraction"
            "<td class=centered>" << segmentPairInformation.unexplainedFraction(0) <<
            "<td class=centered>" << segmentPairInformation.unexplainedFraction(1);
        html.precision(oldPrecision);
        html.flags(oldFlags);
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

}



void Assembler::exploreMode3MetaAlignment(
    const vector<string>& request,
    ostream& html)
{
    // Access the mode 3 assembly graph.
    SHASTA_ASSERT(assemblyGraph3Pointer);
    const mode3::AssemblyGraph& assemblyGraph3 = *assemblyGraph3Pointer;

    // Get the oriented read ids from the request.
    string orientedReadId0String;
    const bool orientedReadId0IsPresent = getParameterValue(request, "orientedReadId0", orientedReadId0String);
    string orientedReadId1String;
    const bool orientedReadId1IsPresent = getParameterValue(request, "orientedReadId1", orientedReadId1String);

    // Write the form.
    html <<
        "Enter the two oriented read ids:"
        "<form>"
        "<p><input type=text size=8 name=orientedReadId0 value='" <<
        (orientedReadId0IsPresent ? orientedReadId0String : "") << "'>"
        "<p><input type=text size=8 name=orientedReadId1 value='" <<
        (orientedReadId1IsPresent ? orientedReadId1String : "") << "'>"
        "<p><input type=submit value='Compute the meta-alignment'>"
        "</form>";

    // If the oriented reads are not present, do nothing.
    if(not(orientedReadId0IsPresent and orientedReadId1IsPresent)) {
        return;
    }
    const OrientedReadId orientedReadId0(orientedReadId0String);
    const OrientedReadId orientedReadId1(orientedReadId1String);

    html << "<h1>Meta-alignment of oriented reads " <<
        orientedReadId0 << " " << orientedReadId1 << "</h1>";

    // Access the pseudo-paths, that is the meta-sequences to be aligned.
    const auto pseudoPath0 = assemblyGraph3.assemblyGraphJourneys[orientedReadId0.getValue()];
    const auto pseudoPath1 = assemblyGraph3.assemblyGraphJourneys[orientedReadId1.getValue()];
    const int n0 = int(pseudoPath0.size());
    const int n1 = int(pseudoPath1.size());

    // Create a png file representing the alignment matrix.
    PngImage image = PngImage(int(n0), int(n1));
    for(int i0=0; i0<n0; i0++) {
        const uint64_t segmentId0 = pseudoPath0[i0].segmentId;
        for(int i1=0; i1<n1; i1++) {
            const uint64_t segmentId1 = pseudoPath1[i1].segmentId;
            if(segmentId0 == segmentId1) {
                image.setPixel(i0, i1, 255, 0, 0);
            }
        }
    }
    image.write("MetaAlignment.png");

    // Create a base64 version of the png file.
    const string command = "base64 MetaAlignment.png > MetaAlignment.png.base64";
    ::system(command.c_str());

    // Display the picture with the alignment.
    // image-rendering:crisp-edges; is currently supported on Firefox but not Chrome,
    // so Chrome will display blurry pictures.
    html <<
        "<h3>Alignment matrix</h3>"
        "<p><img "
        " style='width:" << 3*n0 << "px;height:auto;image-rendering:crisp-edges;'"
        "src=\"data:image/png;base64,";
    ifstream png("MetaAlignment.png.base64");
    html << png.rdbuf();
    html << "\"/>";

}



void Assembler::exploreMode3AssemblyPath(
    const vector<string>& request,
    ostream& html)
{
    // Get the parametersof the request.

    // The segment that the path will start from.
    string pathStartString;
    HttpServer::getParameterValue(request, "pathStart", pathStartString);

    // The path direction can be forward, backward, or bidirectional.
    string pathDirection = "bidirectional";
    HttpServer::getParameterValue(request, "pathDirection", pathDirection);



    // Write the form.
    html <<
        "<h2>Assembly path computation</h2>"
        "<form>"

        "Start the path at segment &nbsp;<input type=text name=pathStart required size=8 style='text-align:center'"
                " value='" << pathStartString << "'>"

        "<br><input type=radio name=pathDirection value=forward" <<
        (pathDirection=="forward" ? " checked=checked" : "") << "> Forward"
        "<br><input type=radio name=pathDirection value=backward" <<
        (pathDirection=="backward" ? " checked=checked" : "") << "> Backward"
        "<br><input type=radio name=pathDirection value=bidirectional" <<
        (pathDirection=="bidirectional" ? " checked=checked" : "") << "> Both directions" <<

        "<p><input type=submit value='Compute the path and assemble its sequence'>"
        "</form>";

    // If the path start was not specified, stop here.
    if(pathStartString.empty()) {
        return;
    }

    // Get the path start segment.
    uint64_t pathStart;
    try {
        pathStart = boost::lexical_cast<uint64_t>(pathStartString);
    } catch(std::exception&) {
        throw runtime_error("Invalid path start segment id.");
    }

    // Check that it is a valid segment id.
    const mode3::AssemblyGraph& assemblyGraph = *assemblyGraph3Pointer;
    if(pathStart >= assemblyGraph.markerGraphPaths.size()) {
        throw runtime_error("Invalid path start segment id. The assembly graph has " +
            to_string(assemblyGraph.markerGraphPaths.size()) + " segments.");
    }

    // Write a header.
    html << "<h1>Assembly path</h1>";



    // Compute the assembly path.
    AssemblyPath path;
    if(pathDirection == "forward" or pathDirection == "backward") {

        // Forward or backward.
        assemblyGraph.createAssemblyPath(pathStart,
            (pathDirection == "forward") ? 0 : 1, path);
        if(pathDirection == "backward") {
            reverse(path.segments.begin(), path.segments.end());
        }

    } else {

        // Bidirectional.
        AssemblyPath forwardPath;
        AssemblyPath backwardPath;
        assemblyGraph.createAssemblyPath(pathStart, 0, forwardPath);
        assemblyGraph.createAssemblyPath(pathStart, 1, backwardPath);

        // Stitch them together, making sure not to repeat the starting segment.
        path.segments.clear();
        copy(backwardPath.segments.rbegin(), backwardPath.segments.rend(), back_inserter(path.segments));
        copy(forwardPath.segments.begin() + 1, forwardPath.segments.end(), back_inserter(path.segments));

    }

    html << "<p>This assembly path was created starting at segment " << pathStart <<
        " and moving ";
    if(pathDirection == "forward") {
        html << "forward.";
    } else if(pathDirection == "backward") {
        html << "backward.";
    } else if(pathDirection == "bidirectional") {
        html << "in both directions.";
    }

    // Assemble sequence for this path.
    path.assemble(assemblyGraph);

    // Write path details to html.
    path.writeHtml(html, assemblyGraph);


}



void Assembler::exploreMode3LinkAssembly(
    const vector<string>& request,
    ostream& html)
{
    // Access the AssemblyGraph.
    using mode3::AssemblyGraph; // Hide shasta::AssemblyGraph;
    SHASTA_ASSERT(assemblyGraph3Pointer);
    const AssemblyGraph& assemblyGraph = *assemblyGraph3Pointer;

    // Get the parameters of the request.
    uint64_t linkId = invalid<uint64_t>;
    getParameterValue(request, "linkId", linkId);
    SHASTA_ASSERT(linkId < assemblyGraph.links.size());
    uint64_t previousPrimarySegmentId = invalid<uint64_t>;
    getParameterValue(request, "previousPrimarySegmentId", previousPrimarySegmentId);
    SHASTA_ASSERT(previousPrimarySegmentId < assemblyGraph.markerGraphPaths.size());
    uint64_t nextPrimarySegmentId = invalid<uint64_t>;
    getParameterValue(request, "nextPrimarySegmentId", nextPrimarySegmentId);
    SHASTA_ASSERT(nextPrimarySegmentId < assemblyGraph.markerGraphPaths.size());

    // Access the link.
    if(linkId >= assemblyGraph.links.size()) {
        html << "Invalid link id. There are " << assemblyGraph.links.size() <<
            " links in the assembly graph.";
        return;
    }
    const AssemblyGraph::Link& link = assemblyGraph.links[linkId];

    // If this is a trivial link, there is nothing to show.
    if(link.segmentsAreAdjacent) {
        html << "This is a trivial link. No assembly is required.";
        return;
    }



    html << "<h1>Details of link assembly</h1>";

    // Create the segments and assemble them.
    AssemblyPathSegment segment0(link.segmentId0, false);
    AssemblyPathSegment segment1(link.segmentId1, false);
    assembleMarkerGraphPath(
        assemblyGraph.readRepresentation,
        assemblyGraph.k,
        assemblyGraph.markers,
        assemblyGraph.markerGraph,
        assemblyGraph.markerGraphPaths[segment0.id],
        false,
        segment0.assembledSegment);
    assembleMarkerGraphPath(
        assemblyGraph.readRepresentation,
        assemblyGraph.k,
        assemblyGraph.markers,
        assemblyGraph.markerGraph,
        assemblyGraph.markerGraphPaths[segment1.id],
        false,
        segment1.assembledSegment);

    // Create the AssemblyPathLink.
    AssemblyPathLink assemblyPathLink;
    assemblyPathLink.id = linkId;
    assemblyPathLink.isTrivial = false;
    assemblyPathLink.previousPrimarySegmentId = previousPrimarySegmentId;
    assemblyPathLink.nextPrimarySegmentId = nextPrimarySegmentId;

    // Do the assembly.
    AssemblyPath::assembleNonTrivialLink(
        assemblyGraph,
        segment0,
        segment1,
        assemblyPathLink,
        html);
}
