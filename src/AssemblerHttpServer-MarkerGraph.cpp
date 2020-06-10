#ifdef SHASTA_HTTP_SERVER

// Shasta.
#include "Assembler.hpp"
#include "ConsensusCaller.hpp"
#include "InducedAlignment.hpp"
#include "LocalMarkerGraph.hpp"
#include "platformDependent.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Spoa.
#include "spoa/spoa.hpp"

// Standard library.
#include "chrono.hpp"
#include "iterator.hpp"



void Assembler::exploreMarkerGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the request parameters.
    LocalMarkerGraphRequestParameters requestParameters;
    getLocalMarkerGraphRequestParameters(request, requestParameters);

    // Write the form and the color legend.
    html << "<h3>Display a local subgraph of the global marker graph</h3>";
    html << "<div style='clear:both; display:table;'>";
    html << "<div style='float:left;margin:10px;'>";
    requestParameters.writeForm(html, markerGraph.vertexCount());
    html << "</div>";
    html << "<div style='float:left;margin:10px;'>";
    LocalMarkerGraph::writeColorLegend(html);
    html << "</div>";
    html << "</div>";

    // If any required values are missing, stop here.
    if(requestParameters.hasMissingRequiredParameters()) {
        return;
    }



    // Validity checks.
    if(requestParameters.vertexId > markerGraph.vertexCount()) {
        html << "<p>Invalid vertex id " << requestParameters.vertexId;
        html << ". Must be between 0 and " << markerGraph.vertexCount()-1 << " inclusive.";
        return;
    }



    // Create the local marker graph.
    LocalMarkerGraph graph(
        uint32_t(assemblerInfo->k),
        reads,
        readRepeatCounts,
        markers,
        markerGraph.vertexTable,
        *consensusCaller);
    const auto createStartTime = steady_clock::now();
    if(!extractLocalMarkerGraphUsingStoredConnectivity(
        requestParameters.vertexId,
        requestParameters.maxDistance,
        requestParameters.timeout,
        requestParameters.useWeakEdges,
        requestParameters.usePrunedEdges,
        requestParameters.useSuperBubbleEdges,
        graph)) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    if(num_vertices(graph) == 0) {
        html << "<p>The local marker graph is empty.";
        return;
    }
    vector< pair<shasta::Base, int> > sequence;
    const auto createFinishTime = steady_clock::now();
    if(requestParameters.timeout>0 && seconds(createFinishTime - createStartTime) > requestParameters.timeout) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }

    html <<
        "<script>\n"
        "function positionAtVertex(vertexId) {\n"
        "var element = document.getElementById('vertex' + vertexId);\n"
        "var r = element.getBoundingClientRect();\n"
        "window.scrollBy((r.left + r.right - window.innerWidth) / 2, (r.top + r.bottom - window.innerHeight) / 2);\n"
        "}\n"
        "</script>\n";



    // Write it out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    graph.write(
        dotFileName,
        requestParameters.maxDistance,
        requestParameters.addLabels,
        requestParameters.useDotLayout,
        requestParameters.vertexScalingFactor,
        requestParameters.arrowScalingFactor);

    // Compute layout in svg format.
    const string command =
        timeoutCommand() + " " + to_string(requestParameters.timeout - int(seconds(createFinishTime - createStartTime))) +
        " dot -O -T svg " + dotFileName +
        " -Gsize=" + to_string(requestParameters.sizePixels/72.);
    const int commandStatus = ::system(command.c_str());
    if(WIFEXITED(commandStatus)) {
        const int exitStatus = WEXITSTATUS(commandStatus);
        if(exitStatus == 124) {
            html << "<p>Timeout for graph layout exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
            filesystem::remove(dotFileName);
            return;
        }
        else if(exitStatus!=0 && exitStatus!=1) {    // sfdp returns 1 all the time just because of the message about missing triangulation.
            // filesystem::remove(dotFileName);
            throw runtime_error("Error " + to_string(exitStatus) + " running graph layout command: " + command);
        }
    } else if(WIFSIGNALED(commandStatus)) {
        const int signalNumber = WTERMSIG(commandStatus);
        throw runtime_error("Signal " + to_string(signalNumber) + " while running graph layout command: " + command);
    } else {
        throw runtime_error("Abnormal status " + to_string(commandStatus) + " while running graph layout command: " + command);

    }
    // Remove the .dot file.
    // filesystem::remove(dotFileName);



    // Write the graph.
    html <<
        "<h2>Marker graph near marker graph vertex " << requestParameters.vertexId <<
        "</h2>";

    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << svgFile.rdbuf();
    svgFile.close();

    // Remove the .svg file.
    filesystem::remove(svgFileName);



    // Make the vertices clickable: left click recenters
    // the graph at that vertex, right click shows vertex details.
    html << "<script>\n";
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        const LocalMarkerGraphVertex& vertex = graph[v];
        SHASTA_ASSERT(!vertex.markerInfos.empty());
        const string url =
            "exploreMarkerGraph?vertexId=" + to_string(vertex.vertexId) +
            "&maxDistance=" + to_string(requestParameters.maxDistance) +
            "&sizePixels=" + to_string(requestParameters.sizePixels) +
            "&timeout=" + to_string(requestParameters.timeout) +
            (requestParameters.addLabels ? "&addLabels=on" : "") +
            "&layout=" + (requestParameters.useDotLayout ? "dot" : "sfdp") +
            (requestParameters.useWeakEdges ? "&useWeakEdges=on" : "") +
            (requestParameters.usePrunedEdges ? "&usePrunedEdges=on" : "") +
            (requestParameters.useSuperBubbleEdges ? "&useSuperBubbleEdges=on" : "");
        html <<
            "element = document.getElementById('vertex" << vertex.vertexId << "');\n"
            "element.onclick = function() {location.href='" << url << "';};\n"
            "element.style.cursor = \"pointer\";\n";

        // Add a right click to show details.
        const string detailUrl =
            "exploreMarkerGraphVertex?vertexId=" + to_string(vertex.vertexId);
        html <<
            "element.oncontextmenu = function() {window.open('" << detailUrl << "');"
            "return false;};\n";
    }
    html << "</script>\n";



    // Make the edges clickable: left click recenters
    // the graph at the source vertex of that edge, right click shows edge details.
    html << "<script>\n";
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph) {
        const LocalMarkerGraphEdge& edge = graph[e];
        const LocalMarkerGraph::vertex_descriptor v0 = source(e, graph);
        const LocalMarkerGraphVertex& vertex0 = graph[v0];
        const string url =
            "exploreMarkerGraph?vertexId=" + to_string(vertex0.vertexId) +
            "&maxDistance=" + to_string(requestParameters.maxDistance) +
            "&sizePixels=" + to_string(requestParameters.sizePixels) +
            "&timeout=" + to_string(requestParameters.timeout) +
            (requestParameters.addLabels ? "&addLabels=on" : "") +
            "&layout=" + (requestParameters.useDotLayout ? "dot" : "sfdp") +
            (requestParameters.useWeakEdges ? "&useWeakEdges=on" : "") +
            (requestParameters.usePrunedEdges ? "&usePrunedEdges=on" : "") +
            (requestParameters.useSuperBubbleEdges ? "&useSuperBubbleEdges=on" : "");
        html <<
            "element = document.getElementById('edge" << edge.edgeId << "');\n"
            "element.onclick = function() {location.href='" << url << "';};\n"
            "element.style.cursor = \"pointer\";\n";

        // Add a right click to show details.
        const string detailUrl =
            "exploreMarkerGraphEdge?edgeId=" + to_string(edge.edgeId);
        html <<
            "element.oncontextmenu = function() {window.open('" << detailUrl << "');"
            "return false;};\n";
    }
    html << "</script>\n";



    // Position the start vertex at the center of the window.
    html <<
        "<script>\n"
        "positionAtVertex(" << requestParameters.vertexId << ");\n"
        "</script>\n";
}



// Extract  from the request the parameters for the display
// of the local marker graph.
void Assembler::getLocalMarkerGraphRequestParameters(
    const vector<string>& request,
    LocalMarkerGraphRequestParameters& parameters) const
{
    parameters.vertexId = 0;
    parameters.vertexIdIsPresent = getParameterValue(
        request, "vertexId", parameters.vertexId);

    parameters.maxDistance = 0;
    parameters.maxDistanceIsPresent = getParameterValue(
        request, "maxDistance", parameters.maxDistance);

    string addLabelsString;
    parameters.addLabels = getParameterValue(
        request, "addLabels", addLabelsString);

    string layoutString;
    getParameterValue(
        request, "layout", layoutString);
    parameters.useDotLayout = true;
    if(layoutString == "sfdp") {
        parameters.useDotLayout = false;
    }

    string useWeakEdgesString;
    parameters.useWeakEdges = getParameterValue(
        request, "useWeakEdges", useWeakEdgesString);

    string usePrunedEdgesString;
    parameters.usePrunedEdges = getParameterValue(
        request, "usePrunedEdges", usePrunedEdgesString);

    string useSuperBubbleEdgesString;
    parameters.useSuperBubbleEdges = getParameterValue(
        request, "useSuperBubbleEdges", useSuperBubbleEdgesString);

    parameters.sizePixels = 800;
    parameters.sizePixelsIsPresent = getParameterValue(
        request, "sizePixels", parameters.sizePixels);

    parameters.vertexScalingFactor = 1;
    parameters.vertexScalingFactorIsPresent = getParameterValue(
        request, "vertexScalingFactor", parameters.vertexScalingFactor);

    parameters.arrowScalingFactor = 1;
    parameters.arrowScalingFactorIsPresent = getParameterValue(
        request, "arrowScalingFactor", parameters.arrowScalingFactor);

    parameters.timeout = 30;
    parameters.timeoutIsPresent = getParameterValue(
        request, "timeout", parameters.timeout);

}



void Assembler::LocalMarkerGraphRequestParameters::writeForm(
    ostream& html,
    MarkerGraph::VertexId vertexCount) const
{
    html <<
        "<form>"

        "<table>"

        "<tr title='Start vertex id between 0 and " << vertexCount << "'>"
        "<td>Start vertex id"
        "<td><input type=text required name=vertexId size=8 style='text-align:center'"
        << (vertexIdIsPresent ? ("value='"+to_string(vertexId)+"'") : "") <<
        ">"

        "<tr title='Maximum distance from start vertex (number of edges)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        << (maxDistanceIsPresent ? ("value='" + to_string(maxDistance)+"'") : " value='6'") <<
        ">"

        "<tr title='Check for to add labels to vertices and edges'>"
        "<td>Labels"
        "<td class=centered><input type=checkbox name=addLabels"
        << (addLabels ? " checked=checked" : "") <<
        ">"

        "<tr title='Check for to add labels to vertices and edges'>"
        "<td>Graph layout"
        "<td>"
        "<span title='Best for small subgraphs'><input type=radio name=layout value=dot"
        << (useDotLayout ? " checked=checked" : "") <<
        ">Dot</span><br>"
        "<span title='Best for large subgraphs, without labels'><input type=radio name=layout value=sfdp"
        << (!useDotLayout ? " checked=checked" : "") <<
        ">Sfdp</span>"

        "<tr title='Check to include in the local marker graph "
        "edges that were removed during transitive reduction'>"
        "<td>Edges removed during transitive reduction"
        "<td class=centered>"
        "<input type=checkbox name=useWeakEdges" <<
        (useWeakEdges ? " checked=checked" : "") << ">"

        "<tr title='Check to include in the local marker graph "
        "edges that were removed during pruning'>"
        "<td>Edges removed during pruning"
        "<td class=centered>"
        "<input type=checkbox name=usePrunedEdges" <<
        (usePrunedEdges ? " checked=checked" : "") << ">"

        "<tr title='Check to include in the local marker graph "
        "edges that were removed during bubble/superbubble detection'>"
        "<td>Edges removed during bubble/superbubble detection"
        "<td class=centered>"
        "<input type=checkbox name=useSuperBubbleEdges" <<
        (useSuperBubbleEdges ? " checked=checked" : "") << ">"

        "<tr title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "<td>Graphics size in pixels"
        "<td><input type=text required name=sizePixels size=8 style='text-align:center'"
        << (sizePixelsIsPresent ? (" value='" + to_string(sizePixels)+"'") : " value='800'") <<
        ">"

        "<tr>"
        "<td>Vertex scaling factor (sfdp only)"
        "<td><input type=text required name=vertexScalingFactor size=8 style='text-align:center'" <<
        " value='" + vertexScalingFactorString() + "'>" <<

        "<tr>"
        "<td>Edge arrow scaling factor"
        "<td><input type=text required name=arrowScalingFactor size=8 style='text-align:center'" <<
        " value='" + arrowScalingFactorString() + "'>" <<


        "<tr title='Maximum time allowed (seconds) for graph creation and layout, or 0 if unlimited'>"
        "<td>Timeout (seconds) for graph creation and layout"
        "<td><input type=text required name=timeout size=8 style='text-align:center'"
        << (timeoutIsPresent ? (" value='" + to_string(timeout)+"'") : " value='30'") <<
        ">"
        "</table>"



        "<br><input type=submit value='Display'>"
        "</form>";
}



bool Assembler::LocalMarkerGraphRequestParameters::hasMissingRequiredParameters() const
{
    return
        !vertexIdIsPresent ||
        !maxDistanceIsPresent ||
        !timeoutIsPresent;
}



string Assembler::LocalMarkerGraphRequestParameters::vertexScalingFactorString() const
{
    if(vertexScalingFactorIsPresent) {
        std::ostringstream s;
        s << vertexScalingFactor;
        return s.str();
    } else {
        return "1";
    }
}



string Assembler::LocalMarkerGraphRequestParameters::arrowScalingFactorString() const
{
    if(arrowScalingFactorIsPresent) {
        std::ostringstream s;
        s << arrowScalingFactor;
        return s.str();
    } else {
        return "1";
    }
}



void Assembler::exploreMarkerGraphVertex(const vector<string>& request, ostream& html)
{
    // Get the vertex id.
    MarkerGraph::VertexId vertexId = 0;
    const bool vertexIdIsPresent = getParameterValue(request, "vertexId", vertexId);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Show details for marker graph vertex'> "
        "<input type=text name=vertexId required" <<
        (vertexIdIsPresent ? (" value=" + to_string(vertexId)) : "") <<
        " size=8 title='Enter a vertex id between 0 and " << markerGraph.vertexCount()-1 << "'>";
    html << "</form>";

    // If the vertex id missing or invalid, stop here.
    if(!vertexIdIsPresent || !vertexIdIsPresent) {
        return;
    }
    if(vertexId >= markerGraph.vertexCount()) {
        html << "<p>Invalid vertex id. Must be less than " << markerGraph.vertexCount() << ".";
        return;
    }

    // Access the markers of this vertex.
    span<MarkerId> markerIds = markerGraph.vertices[vertexId];
    const size_t markerCount = markerIds.size();
    SHASTA_ASSERT(markerCount > 0);

    // Get the marker sequence.
    const KmerId kmerId = markers.begin()[markerIds[0]].kmerId;
    const size_t k = assemblerInfo->k;
    const Kmer kmer(kmerId, k);



    // Extract the information we need.
    vector<OrientedReadId> orientedReadIds(markerCount);
    vector<uint32_t> ordinals(markerCount);
    vector< vector<uint8_t> > repeatCounts(markerCount, vector<uint8_t>(k));
    for(size_t j=0; j<markerCount; j++) {
        const MarkerId markerId = markerIds[j];
        const CompressedMarker& marker = markers.begin()[markerId];
        tie(orientedReadIds[j], ordinals[j]) = findMarkerId(markerId);

        // Get the repeat count for this marker at each of the k positions.
        for(size_t i=0; i<k; i++) {
            Base base;
            tie(base, repeatCounts[j][i]) =
                getOrientedReadBaseAndRepeatCount(orientedReadIds[j], uint32_t(marker.position+i));
            SHASTA_ASSERT(base == kmer[i]);
        }
    }



    // Find all the repeat counts represented.
    std::set<size_t> repeatCountsSet;
    for(const auto& v: repeatCounts) {
        for(const auto r: v) {
            repeatCountsSet.insert(r);
        }
    }



    // Compute consensus repeat counts at each of the k positions.
    vector<size_t> consensusRepeatCounts(k);
    const auto storedConsensusRepeatCounts =
        markerGraph.vertexRepeatCounts.begin() + k * vertexId;
    for(size_t i=0; i<k; i++) {

        Coverage coverage;
        for(size_t j=0; j<markerCount; j++) {
            coverage.addRead(
                AlignedBase(kmer[i]),
                orientedReadIds[j].getStrand(),
                repeatCounts[j][i]);
        }

        const Consensus consensus = (*consensusCaller)(coverage);
        SHASTA_ASSERT(Base(consensus.base) == kmer[i]);
        consensusRepeatCounts[i] = consensus.repeatCount;

        // Check that this repeat count agrees with what was
        // computed during the assembly.
        if(consensusRepeatCounts[i] != storedConsensusRepeatCounts[i]) {
            html << "<p><b>Stored consensus repeat counts do not agree with "
                "the values computed on the fly.</b>" << endl;
        }
    }




    // Compute concordant and discordant coverage at each position.
    vector<size_t> concordantCoverage(k, 0);
    vector<size_t> discordantCoverage(k, 0);
    for(size_t i=0; i<k; i++) {
        for(size_t j=0; j<markerCount; j++) {
            if(repeatCounts[j][i] == consensusRepeatCounts[i]) {
                ++concordantCoverage[i];
            } else {
                ++discordantCoverage[i];
            }
        }
    }


    // Page title.
    const string titleUrl =
        "exploreMarkerGraph?vertexId=" + to_string(vertexId) +
        "&maxDistance=3"
        "&useWeakEdges=on"
        "&usePrunedEdges=on"
        "&useSuperBubbleEdges=on"
        "&sizePixels=800"
        "&timeout=30";
    html << "<h1>Marker graph vertex <a href='" << titleUrl << "'> "<< vertexId << "</a></h1>";



    // Table with summary information for this vertex.
    html <<
        "<table>"
        "<tr><th class=left>Coverage<td class=centered>" << markerCount <<
        "<tr><th class=left>Marker sequence (run-length)" <<
        "<td class=centered style='font-family:monospace'>";
    kmer.write(html, assemblerInfo->k);

    // Write a row with consensus repeat counts.
    html <<
        "<tr><th class=left>Consensus repeat counts"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<k; i++) {
        const size_t repeatCount = consensusRepeatCounts[i];
        if(repeatCount < 10) {
            html << repeatCount;
        } else {
            html << "*";
        }
    }

    // Write a row with the consensus raw sequence.
    html <<
        "<tr><th class=left>Consensus raw sequence"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<k; i++) {
        const Base base = kmer[i];
        const size_t repeatCount = consensusRepeatCounts[i];
        for(size_t r=0; r<repeatCount; r++) {
            html << base;
        }
    }

    // Write a row with outgoing edges.
    html <<
        "<tr><th class=left>Edges starting at this vertex<td class=centered>";
    for(const MarkerGraph::EdgeId edgeId: markerGraph.edgesBySource[vertexId]) {
        html << "<a href='exploreMarkerGraphEdge?edgeId=" << edgeId << "'";
        if(markerGraph.edges[edgeId].wasRemoved()) {
            html << " style='color:#a8b9ea'";
        }
        html << ">" << edgeId << "</a> ";
    }
    html <<
        "<tr><th class=left>Edges ending at this vertex<td class=centered>";
    for(const MarkerGraph::EdgeId edgeId: markerGraph.edgesByTarget[vertexId]) {
        html << "<a href='exploreMarkerGraphEdge?edgeId=" << edgeId << "'";
        if(markerGraph.edges[edgeId].wasRemoved()) {
            html << " style='color:#a8b9ea'";
        }
        html << ">" << edgeId << "</a> ";
    }


    html << "</table>";



    // Table with one row for each marker.
    html <<
        "<h3>Markers of this vertex</h3>"
        "<table>"
        "<tr>"
        "<th>Oriented<br>read"
        "<th title='Ordinal of the marker in the oriented read'>Ordinal"
        "<th>Repeat<br>counts";
    for(size_t j=0; j<markerIds.size(); j++) {
        const OrientedReadId orientedReadId = orientedReadIds[j];
        const ReadId readId = orientedReadId.getReadId();
        const Strand strand = orientedReadId.getStrand();
        const uint32_t ordinal = ordinals[j];

        // Oriented read id.
        html <<
            "<tr>"
            "<td class=centered>"
            "<a href='exploreRead"
            "?readId=" << readId <<
            "&strand=" << strand << "'>" <<
            orientedReadId << "</a>";

        // Marker ordinal.
        html <<
            "<td class=centered>" <<
            "<a href='exploreRead"
            "?readId=" << readId <<
            "&strand=" << strand <<
            "&highlightMarker=" << ordinal <<
            "'>" <<
            ordinal << "</a>";

        // Repeat counts.
        html << "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<k; i++) {
            const uint8_t repeatCount = repeatCounts[j][i];
            if(repeatCount < 10) {
                html << int(repeatCount);
            } else {
                html << "*";
            }
        }
    }



    // Write a row with consensus repeat counts.
    html <<
        "<tr><th colspan=2 class=left>Consensus repeat counts"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<k; i++) {
        const size_t repeatCount = consensusRepeatCounts[i];
        if(repeatCount < 10) {
            html << repeatCount;
        } else {
            html << "*";
        }
    }



    // Write a row with the marker sequence (run-length).
    html <<
        "<tr><th colspan=2 class=left>Run-length sequence"
        "<td class=centered style='font-family:monospace'>";
    kmer.write(html, assemblerInfo->k);



    // Write rows with coverage information for each represented repeat value.
    for(size_t repeatCount: repeatCountsSet) {
        html <<
            "<tr><th colspan=2 class=left>Coverage for repeat count " << repeatCount <<
            "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<k; i++) {

            // Compute coverage for this repeat count, at this position.
            size_t coverage = 0;
            for(size_t j=0; j<markerCount; j++) {
                if(repeatCounts[j][i] == repeatCount) {
                    coverage++;
                }
            }

            // Write it out.
            if(coverage == 0) {
                html << ".";
            } else if(coverage < 10) {
                html << coverage;
            } else {
                html << "*";
            }
        }
    }



    // Write a row with concordant coverage.
    html <<
        "<tr><th colspan=2 class=left>Concordant repeat count coverage"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<assemblerInfo->k; i++) {
        const size_t coverage = concordantCoverage[i];
        if(coverage == 0) {
            html << ".";
        } else if(coverage < 10) {
            html << coverage;
        } else {
            html << "*";
        }
    }



    // Write a row with discordant coverage.
    html <<
        "<tr><th colspan=2 class=left>Discordant repeat count coverage"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<assemblerInfo->k; i++) {
        const size_t coverage = discordantCoverage[i];
        if(coverage == 0) {
            html << ".";
        } else if(coverage < 10) {
            html << coverage;
        } else {
            html << "*";
        }
    }



    // Write a row with the consensus raw sequence.
    html <<
        "<tr><th colspan=2 class=left>Consensus raw sequence"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<k; i++) {
        const Base base = kmer[i];
        const size_t repeatCount = consensusRepeatCounts[i];
        for(size_t r=0; r<repeatCount; r++) {
            html << base;
        }
    }

    html << "</table>";
}



void Assembler::exploreMarkerGraphEdge(const vector<string>& request, ostream& html)
{
    // Get the edge id.
    MarkerGraph::EdgeId edgeId = 0;
    const bool edgeIdIsPresent = getParameterValue(request, "edgeId", edgeId);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Show details for marker graph edge'> "
        "<input type=text name=edgeId required" <<
        (edgeIdIsPresent ? (" value=" + to_string(edgeId)) : "") <<
        " size=8 title='Enter an edge id between 0 and " << markerGraph.edges.size()-1 << "'>";
    html << "</form>";

    // If the edge id missing or invalid, stop here.
    if(!edgeIdIsPresent || !edgeIdIsPresent) {
        return;
    }
    if(edgeId >= markerGraph.edges.size()) {
        html << "<p>Invalid edge id. Must be less than " << markerGraph.edges.size() << ".";
        return;
    }

    // Access the edge.
    const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
    array<MarkerGraph::VertexId, 2> vertexIds = {edge.source, edge.target};
    const size_t markerCount = edge.coverage;

    // The marker intervals of this edge.
    const span<MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
    SHASTA_ASSERT(markerIntervals.size() == markerCount);

    // The length of each marker sequence.
    const size_t k = assemblerInfo->k;



    // Extract the sequences and repeat counts for each marker interval.
    // This includes sequences and repeat counts for the flanking markers.
    vector< vector<Base> > sequences(markerCount);
    vector< vector<uint8_t> > repeatCounts(markerCount);
    for(size_t j=0; j!=markerCount; j++) {
        const MarkerInterval& markerInterval = markerIntervals[j];
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Get the two markers.
        const CompressedMarker& marker0 = orientedReadMarkers[markerInterval.ordinals[0]];
        const CompressedMarker& marker1 = orientedReadMarkers[markerInterval.ordinals[1]];

        // Get the position range, including the flanking markers.
        const uint32_t positionBegin = marker0.position;
        const uint32_t positionEnd = marker1.position + uint32_t(k);

        // Store the bases and repeat counts in this interval.
        for(uint32_t position=positionBegin; position!=positionEnd; ++position) {
            Base base;
            uint8_t repeatCount;
            tie(base, repeatCount) = getOrientedReadBaseAndRepeatCount(orientedReadId, position);
            sequences[j].push_back(base);
            repeatCounts[j].push_back(repeatCount);
        }
    }


    // Access stored consensus for this edge.
    const int storedConsensusOverlappingBaseCount = int(markerGraph.edgeConsensusOverlappingBaseCount[edgeId]);
    const auto storedConsensus = markerGraph.edgeConsensus[edgeId];



    // To compute consensus, use the same code used during assembly.
    const uint32_t markerGraphEdgeLengthThresholdForConsensus = 1000;
    vector<Base> spoaSequence;
    vector<uint32_t> spoaRepeatCounts;
    uint8_t spoaOverlappingBaseCount;
    ComputeMarkerGraphEdgeConsensusSequenceUsingSpoaDetail spoaDetail;
    
    const spoa::AlignmentType alignmentType = spoa::AlignmentType::kNW;
    const int8_t match = 1;
    const int8_t mismatch = -1;
    const int8_t gap = -1;
    auto spoaAlignmentEngine = spoa::createAlignmentEngine(alignmentType, match, mismatch, gap);
    auto spoaAlignmentGraph = spoa::createGraph();
    
    computeMarkerGraphEdgeConsensusSequenceUsingSpoa(
        edgeId,
        markerGraphEdgeLengthThresholdForConsensus,
        spoaAlignmentEngine,
        spoaAlignmentGraph,
        spoaSequence,
        spoaRepeatCounts,
        spoaOverlappingBaseCount,
        spoaDetail,
        0);



    // Page title.
    const string titleUrl =
        "exploreMarkerGraph?vertexId=" + to_string(vertexIds[0]) +
        "&maxDistance=3"
        "&useWeakEdges=on"
        "&usePrunedEdges=on"
        "&useSuperBubbleEdges=on"
        "&sizePixels=800"
        "&timeout=30";
    html << "<h1>Marker graph edge <a href='" << titleUrl << "'> "<< edgeId << "</a></h1>";



    // Table to summarize this edge.
    html <<
        "<table style='display:block;white-space:nowrap;'>"
        "<tr><th class=left>Source vertex<td class=centered>" <<
        "<a href='exploreMarkerGraphVertex?vertexId=" << vertexIds[0] << "'>" << vertexIds[0] << "</a>"
        "<tr><th class=left>Target vertex<td class=centered>" <<
        "<a href='exploreMarkerGraphVertex?vertexId=" << vertexIds[1] << "'>" << vertexIds[1] << "</a>"
        "<tr><th class=left>Coverage<td class=centered>" << markerCount <<
        "<tr><th class=left>Removed during transitive reduction?<td class=centered>" <<
        (edge.wasRemovedByTransitiveReduction ? "Yes" : "No") <<
        "<tr><th class=left>Removed during pruning?<td class=centered>" <<
        (edge.wasPruned ? "Yes" : "No") <<
        "<tr><th class=left>Removed during bubble/superbubble removal?<td class=centered>" <<
        (edge.isSuperBubbleEdge ? "Yes" : "No");

    if(edge.wasAssembled) {
        html << "<tr><th class=left>Assembled?<td class=centered>Yes" ;
        if(storedConsensusOverlappingBaseCount>0 || storedConsensus.size()==0) {
            html <<
                "<tr><th class=left>Stored consensus: number of overlapping bases<td class=centered>" <<
                storedConsensusOverlappingBaseCount;
        } else {
            html <<
                "<tr><th class=left>Stored consensus: sequence (run-length)<td class=centered style='font-family:monospace'>";
            for(size_t i=0; i<storedConsensus.size(); i++) {
                html << storedConsensus[i].first;
            }
            html <<
                "<tr><th class=left>Stored consensus: repeat counts<td class=centered style='font-family:monospace'>";
            for(size_t i=0; i<storedConsensus.size(); i++) {
                if(storedConsensus[i].second >= 10) {
                    html << "*";
                } else {
                    html << int(storedConsensus[i].second);
                }
            }
            html <<
                "<tr><th class=left>Stored consensus: sequence (raw)<td class=centered style='font-family:monospace'>";
            for(size_t i=0; i<storedConsensus.size(); i++) {
                const Base base = storedConsensus[i].first;
                const int n = int(storedConsensus[i].second);
                for(int j=0; j<n; j++) {
                    html << base;
                }
            }
        }
    } else {
        html << "<tr><th class=left>Assembled?<td class=centered>No";
    }

    html << "</table>";



    // Assembly details section.
    html <<
        "<h3>Assembly details for this marker graph edge</h3>";
    if(spoaDetail.hasLongMarkerInterval) {
        html <<
            "<p>This edge has a long marker interval. "
            "Alignment computation was not performed. "
            "Instead, the consensus was taken equal to the edge sequence "
            "for the oriented read with the shortest marker interval, highlighted below.";
    } else if(spoaDetail.assemblyMode == 1) {
        html << "<p>This edge is dominated by oriented reads with "
            "overlapping bases between the flanking markers. "
            "Reads with one or more intervening bases between "
            "the flanking markers were not used.";
    } else {
        SHASTA_ASSERT(spoaDetail.assemblyMode == 2);
        html << "<p>This edge is dominated by oriented reads with "
            "intervening bases between the flanking markers. "
            "The run-length sequences of oriented reads are shown aligned. "
            "Reads with one or more overlapping bases between "
            "the flanking markers were not used.";
    }



    // Table of assembly details.
    html <<
        "<p><table>"
        "<tr>"
        "<th rowspan=2>Oriented<br>read"
        "<th colspan=2>Flanking<br>markers"
        "<th rowspan=2>Markers<br>skipped"
        "<th rowspan=2>Used<br>for<br>assembly?"
        "<th colspan=2>Sequence<br>between<br>flanking<br>markers<br>(run-length)"
        "<th rowspan=2>Sequence<br>intervening<br>between<br>flanking<br>markers<br>(raw)"
        "<tr>"
        "<th>Left"
        "<th>Right"
        "<th>Overlapping"
        "<th>Intervening";

    // Write one row for each oriented read.
    for(size_t j=0; j<markerCount; j++) {
        const MarkerInterval& markerInterval = markerIntervals[j];
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const ReadId readId = orientedReadId.getReadId();
        const Strand strand = orientedReadId.getStrand();
        const bool hasOverlappingBases = (sequences[j].size() < 2*k);
        const bool hasInterveningBases = (sequences[j].size() > 2*k);

        // Figure out if this read was used for assembly.
        bool wasUsed;
        if(spoaDetail.hasLongMarkerInterval) {
            wasUsed = (j == spoaDetail.iShortest);
        } else if(spoaDetail.assemblyMode == 1) {
            wasUsed = !hasInterveningBases;
        } else {
            SHASTA_ASSERT(spoaDetail.assemblyMode == 2);
            wasUsed = !hasOverlappingBases;
        }

        // Oriented read id.
        html <<
            "<tr>"
            "<td class=centered>"
            "<a href='exploreRead"
            "?readId=" << readId <<
            "&strand=" << strand << "'>" <<
            orientedReadId << "</a>";

        // Ordinals.
        for(size_t m=0; m<2; m++) {
            html <<
                "<td class=centered>" <<
                "<a href='exploreRead"
                "?readId=" << readId <<
                "&strand=" << strand <<
                "&highlightMarker=" << markerInterval.ordinals[m] <<
                "'>" <<
                markerInterval.ordinals[m] << "</a>";
        }

        // Ordinal skip.
        html << "<td class=centered>";
        if(spoaDetail.hasLongMarkerInterval && j == spoaDetail.iShortest) {
            html << "<span style='background-color:yellow'>";
        }
        html << markerInterval.ordinals[1] - 1 - markerInterval.ordinals[0];
        if(spoaDetail.hasLongMarkerInterval && j == spoaDetail.iShortest) {
            html << "</span>";
        }

        // Check if this read was used for assembly.
        html << "<td class=centered>";
        if(wasUsed) {
            html << "&check;";
        }


        // Number of overlapping bases.
        if(sequences[j].size() <= 2*k) {
            html << "<td class=centered>";
            const size_t overlappingBaseCount = 2*k - sequences[j].size();
            if(overlappingBaseCount > 0) {
                html << overlappingBaseCount;
            }
            html << "<td><td>";
        } else {
            html << "<td>";



            // Sequence (run-length). For assembly mode 2,
            // it is written aligned.
            html << "<td class=centered style='font-family:monospace;white-space:nowrap'>";
            if(spoaDetail.hasLongMarkerInterval || spoaDetail.assemblyMode == 1) {

                // Write run-length sequence, without worrying about aligning it.
                for(size_t index=k; index<sequences[j].size()-k; index++) {
                    html << sequences[j][index];
                }
                html << "<br>";
                for(size_t index=k; index<sequences[j].size()-k; index++) {
                    html << int(repeatCounts[j][index]);
                }
            } else {

                // Write alignment sequence, aligned.
                const int row = spoaDetail.alignmentRow[j];
                SHASTA_ASSERT(row>=0 && row<int(spoaDetail.msa.size()));
                const string& msaRow = spoaDetail.msa[row];
                size_t position = k;
                for(const char c: msaRow) {
                    const AlignedBase base = AlignedBase::fromCharacter(c);
                    if(base.isGap()) {
                        html << "-";
                    } else {
                        html << sequences[j][position++];
                    }
                }
                html << "<br>";
                position = k;
                for(const char c: msaRow) {
                    const AlignedBase base = AlignedBase::fromCharacter(c);
                    if(base.isGap()) {
                        html << "-";
                    } else {
                        const uint32_t repeatCount = repeatCounts[j][position++];
                        if(repeatCount < 10) {
                            html << int(repeatCount);
                        } else {
                            html << "*";
                        }
                    }
                }
            }



            // Sequence (raw).
            html << "<td class=centered style='font-family:monospace;white-space:nowrap'>";
            for(size_t index=k; index<sequences[j].size()-k; index++) {
                const size_t repeatCount = repeatCounts[j][index];
                for(size_t l=0; l<repeatCount; l++) {
                    html << sequences[j][index];
                }
            }
        }
    }



    // Write one row with consensus sequence.
    html << "<tr><th colspan=5>Consensus";

    // Consensus sequence (run-length).
    if(spoaDetail.hasLongMarkerInterval || spoaDetail.assemblyMode==1) {
        html << "<td class=centered>";
        if(spoaOverlappingBaseCount > 0) {
            html << int(spoaOverlappingBaseCount);
        }
        html << "<td class=centered style='font-family:monospace'>";
        for(const Base base: spoaSequence) {
            html << base;
        }
        html << "<br>";
        for(const uint32_t repeatCount: spoaRepeatCounts) {
            html << repeatCount;
        }
    } else {
        SHASTA_ASSERT(spoaDetail.assemblyMode==2);
        html <<
            "<td class=centered>" << spoaOverlappingBaseCount <<
            "<td class=centered style='font-family:monospace'>";
        for(const AlignedBase base: spoaDetail.alignedConsensus) {
            html << base;
        }
        html << "<br>";
        for(const uint8_t repeatCount: spoaDetail.alignedRepeatCounts) {
            if(repeatCount == 0) {
                html << "-";
            } else {
                html << int(repeatCount);
            }
        }
    }

    // Write consensus sequence (raw).
    html << "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<spoaSequence.size(); i++) {
        const Base base = spoaSequence[i];
        const int repeatCount = spoaRepeatCounts[i];
        for(int l=0; l<repeatCount; l++) {
            html << base;
        }
    }

    html << "</table>";

}



void Assembler::exploreMarkerGraphInducedAlignment(
    const vector<string>& request,
    ostream& html)
{
    html <<
        "<h1>Display the induced alignment matrix of two oriented reads</h1>"
        "<p>The marker graph induces an effective alignment between each pair "
        "of oriented reads which can be obtained by following each of the oriented reads "
        "in the marker graph. Aligned markers are those that are on the same vertex. "
        "The induced alignment matrix of two oriented reads <i>x</i> and <i>y</i> "
        "with <i>n<sub>x</sub></i> and <i>n<sub>y</sub></i> markers is an "
        "<i>n<sub>x</sub></i>&times;<i>n<sub>y</sub></i> matrix. "
        "Element <i>ij</i> of the matrix is 1 if marker <i>i</i> of <i>x</i> "
        "and marker <i>j</i> of <i>y</i> "
        "are on the same marker graph vertex and 0 otherwise.";

    // Get the read ids and strands from the request.
    ReadId readId0 = 0;
    const bool readId0IsPresent = getParameterValue(request, "readId0", readId0);
    Strand strand0 = 0;
    const bool strand0IsPresent = getParameterValue(request, "strand0", strand0);
    ReadId readId1 = 0;
    const bool readId1IsPresent = getParameterValue(request, "readId1", readId1);
    Strand strand1 = 0;
    const bool strand1IsPresent = getParameterValue(request, "strand1", strand1);
    string ordinalType = "ordinals";
    getParameterValue(request, "ordinalType", ordinalType);

    // Write the form.
    html <<
        "<p>Display the induced alignment matrix of these two reads:"
        "<form>"
        "<input type=text name=readId0 required size=8 " <<
        (readId0IsPresent ? "value="+to_string(readId0) : "") <<
        " title='Enter a read id between 0 and " << reads.size()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand0", strand0IsPresent && strand0==0, strand0IsPresent && strand0==1);
    html <<
        "<br><input type=text name=readId1 required size=8 " <<
        (readId1IsPresent ? "value="+to_string(readId1) : "") <<
        " title='Enter a read id between 0 and " << reads.size()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand1", strand1IsPresent && strand1==0, strand1IsPresent && strand1==1);

    html <<
        "<p>Plot alignment matrix using "
        "<input type=radio required name=ordinalType value='ordinals'" <<
        (ordinalType == "ordinals" ? " checked=on" : "") <<
        ">ordinals "
        "<input type=radio required name=ordinalType value='compressedOrdinals'" <<
        (ordinalType == "compressedOrdinals" ? " checked=on" : "") <<
        ">compressed ordinals"
        "<p><input type=submit value='Display induced alignment'></form>";


     // If the readId's or strand's are missing, stop here.
     if(!readId0IsPresent || !strand0IsPresent || !readId1IsPresent || !strand1IsPresent) {
         return;
     }

     // Compute the induced alignment.
     const OrientedReadId orientedReadId0(readId0, strand0);
     const OrientedReadId orientedReadId1(readId1, strand1);
     InducedAlignment inducedAlignment;
     computeInducedAlignment(orientedReadId0, orientedReadId1, inducedAlignment);
     fillCompressedOrdinals(orientedReadId0, orientedReadId1, inducedAlignment);
     html << "The induced alignment has " << inducedAlignment.data.size() <<
         " marker pairs.";

     // Write the alignment matrix to a png file.
     inducedAlignment.writePngImage(
         uint32_t(markers.size(orientedReadId0.getValue())),
         uint32_t(markers.size(orientedReadId1.getValue())),
         ordinalType == "compressedOrdinals",
         "Alignment.png");

     // Create a base64 version of the png file.
     const string command = "base64 Alignment.png > Alignment.png.base64";
     ::system(command.c_str());

     // Write the image to html.
     html <<
         "<p><img id=\"alignmentMatrix\" onmousemove=\"updateTitle(event)\" "
         "src=\"data:image/png;base64,";
         ifstream png("Alignment.png.base64");
         html << png.rdbuf();
         html << "\"/>"
             "<script>"
             "function updateTitle(e)"
             "{"
             "    var element = document.getElementById(\"alignmentMatrix\");"
             "    var rectangle = element.getBoundingClientRect();"
             "    var x = e.clientX - Math.round(rectangle.left);"
             "    var y = e.clientY - Math.round(rectangle.top);"
             "    element.title = " <<
             "\"" << orientedReadId0 << " marker \" + x + \", \" + "
             "\"" << orientedReadId1 << " marker \" + y;"
             "}"
             "</script>";



     // Write the induced alignment in a table.
     html <<
         "<p><table><tr><th>Vertex"
         "<th>Ordinal<br> in " << orientedReadId0 <<
         "<th>Ordinal<br> in " << orientedReadId1 <<
         "<th>Compressed<br>ordinal<br> in " << orientedReadId0 <<
         "<th>Compressed<br>ordinal<br> in " << orientedReadId1;
    for(const InducedAlignmentData& d: inducedAlignment.data) {
         html <<
             "<tr><td class=centered>" << d.vertexId <<
             "<td class=centered>" << d.ordinal0 <<
             "<td class=centered>" << d.ordinal1 <<
             "<td class=centered>" << d.compressedOrdinal0 <<
             "<td class=centered>" << d.compressedOrdinal1;
     }
     html << "</table>";



}



void Assembler::exploreMarkerCoverage(
    const vector<string>& request,
    ostream&html)
{
    html <<
        "<h1>Marker coverage of an oriented read</h1>"
        "<p>For each marker of a read, marker coverage is defined "
        "as the coverage (number of oriented reads) for the "
        "marker graph vertex associated with that marker, "
        "or zero if there is no marker graph vertex "
        "associated with the marker.";

    // Get the parameters from the request.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);
    int width = 600;
    getParameterValue(request, "width", width);
    int height = 400;
    getParameterValue(request, "height", height);

    // Write the form.
    html <<
        "<form><table>"
        "<tr><td>Read id<td class=centered>"
        "<input type=text name=readId required style='text-align:center'" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " size=8 title='Enter a read id between 0 and " << reads.size()-1 << "'>"
        "<tr><td>Strand<td class=centered>";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);
    html <<
        "<tr><td>Plot width<td class=centered>"
        "<input type=text name=width style='text-align:center' size=8 value='" << width << "'>"
        "<tr><td>Plot height<td class=centered>"
        "<input type=text name=height style='text-align:center' size=8 value='" << height << "'>"
        "</table><input type=submit value='Plot'></form>";

    // If the readId or strand are missing, stop here.
    if(!readIdIsPresent || !strandIsPresent) {
        return;
    }

    const OrientedReadId orientedReadId(readId, strand);
    html << "<h2>Marker coverage of oriented read " << orientedReadId << "</h2>";

    std::ostringstream gnuplotCommands;
    gnuplotCommands <<
        "set border linewidth 1\n"
        "set xtics out nomirror\n"
        "set mxtics 10\n"
        "set ytics out nomirror\n"
        "set grid xtics mxtics ytics linestyle 1 linewidth 1 linecolor rgb '#e0e0e0'\n"
        "plot '-' with points pointtype 7 pointsize 0.5 linecolor rgb '#0000ff' notitle\n";

    const uint32_t markerCount = uint32_t(markers.size(orientedReadId.getValue()));
    for(uint32_t ordinal=0; ordinal<markerCount; ordinal++) {
        const MarkerGraph::VertexId vertexId =
            getGlobalMarkerGraphVertex(orientedReadId, ordinal);
        if(vertexId == MarkerGraph::invalidCompressedVertexId) {
            gnuplotCommands << ordinal << " " << "0\n";
        } else {
            const uint64_t coverage = markerGraph.vertexCoverage(vertexId);
            gnuplotCommands << ordinal << " " << coverage << "\n";
        }
    }

    gnuplotCommands << "e\n";
    writeGnuPlotPngToHtml(html, width, height, gnuplotCommands.str());
}


#endif

