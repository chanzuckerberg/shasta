// Shasta.
#include "Assembler.hpp"
#include "AssemblyGraph.hpp"
#include "compressAlignment.hpp"
#include "ConsensusCaller.hpp"
#include "Coverage.hpp"
#include "hsv.hpp"
#include "InducedAlignment.hpp"
#include "LocalMarkerGraph.hpp"
#include "MarkerConnectivityGraph.hpp"
#include "MurmurHash2.hpp"
#include "platformDependent.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Spoa.
#include "spoa/spoa.hpp"

// Standard library.
#include "chrono.hpp"
#include "filesystem.hpp"
#include "fstream.hpp"
#include "iterator.hpp"
#include <queue>



void Assembler::exploreMarkerGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the request parameters.
    LocalMarkerGraphRequestParameters requestParameters;
    getLocalMarkerGraphRequestParameters(request, requestParameters);

    // Write the form.
    html << "<h1>Display a local subgraph of the global marker graph</h3>";
    requestParameters.writeForm(html, markerGraph.vertexCount());

    // If any required values are missing, stop here.
    if(requestParameters.hasMissingRequiredParameters()) {
        return;
    }

    // Some sanity checks.
    if(requestParameters.vertexRedCoverage >= requestParameters.vertexGreenCoverage) {
        html << "</div><span style='color:purple'>"
            "Red vertex coverage must be less than green vertex coverage.</span></div>";
        return;
    }
    if(requestParameters.edgeRedCoverage >= requestParameters.edgeGreenCoverage) {
        html << "</div><span style='color:purple'>"
            "Red edge coverage must be less than green edge coverage.</span></div>";
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
        assemblerInfo->readRepresentation,
        uint32_t(assemblerInfo->k),
        assemblerInfo->assemblyMode,
        getReads(),
        markers,
        markerGraph.vertexTable,
        *consensusCaller);
    const auto createStartTime = steady_clock::now();
    if(!extractLocalMarkerGraph(
        requestParameters.vertexId,
        requestParameters.maxDistance,
        requestParameters.timeout,
        requestParameters.minVertexCoverage,
        requestParameters.minEdgeCoverage,
        requestParameters.useWeakEdges,
        requestParameters.usePrunedEdges,
        requestParameters.useSuperBubbleEdges,
        requestParameters.useLowCoverageCrossEdges,
        requestParameters.useRemovedSecondaryEdges,
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
    graph.write(dotFileName, requestParameters);


    // Compute layout in svg format.
    const string command =
        timeoutCommand() + " " + to_string(requestParameters.timeout - int(seconds(createFinishTime - createStartTime))) +
        " dot -O -T svg " + dotFileName;
    const int commandStatus = ::system(command.c_str());
    if(WIFEXITED(commandStatus)) {
        const int exitStatus = WEXITSTATUS(commandStatus);
        if(exitStatus == 124) {
            html << "<p>Timeout for graph layout exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
            std::filesystem::remove(dotFileName);
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
    std::filesystem::remove(dotFileName);



    // Write the title.
    html <<
        "<h2>Marker graph near marker graph vertex " << requestParameters.vertexId <<
        "</h2>";



    // Color legend for vertices when colored by distance.
    if(requestParameters.vertexColoring == "byDistance") {
        html << "<h3>Color legend for vertices</h3>";
        LocalMarkerGraph::writeColorLegendVerticesByDistance(html);
    }



    // Color legend for vertices when colored by coverage.
    if(requestParameters.vertexColoring == "byCoverage") {
        html <<
            "<h3>Color legend for vertices</h3><table>"
            "<tr><td class=left>Coverage";
        for(uint64_t coverage=requestParameters.vertexRedCoverage;
            coverage<= requestParameters.vertexGreenCoverage; coverage++) {
            html << "<td class=centered>";
            if(coverage == requestParameters.vertexRedCoverage) {
                html << "&leq;";
            }
            if(coverage == requestParameters.vertexGreenCoverage) {
                html << "&geq;";
            }
            html << coverage;
        }
        html << "<tr><td class=left>Color";
        for(uint64_t coverage=requestParameters.vertexRedCoverage;
            coverage<= requestParameters.vertexGreenCoverage; coverage++) {
            double h =
                double(coverage - requestParameters.vertexRedCoverage) /
                double(requestParameters.vertexGreenCoverage - requestParameters.vertexRedCoverage);
            h = max(h, 0.);
            h = min(h, 1.);
            const double hue = 120. * h;
            double S, L;
            tie(S, L) = hsvToHsl(1., 0.9);
            html << "<td style='height:15px;background-color:hsl(" << hue << "," <<
                uint64_t(100. * S) << "%," <<
                uint64_t(100. * L) << "%)'>";
        }

        html << "</table>";
    }



    // Color legends for edges when colored by flags.
    if(requestParameters.edgeColoring == "byFlags") {
        if(requestParameters.highlightedOrientedReads.empty()) {
            html << "<h3>Color legend for edges lines and arrows</h3>";
            graph.writeColorLegendEdgeArrowsByFlags(html);
        }
        if(requestParameters.edgeLabels > 0) {
            html << "<h3>Color legend for edge labels</h3>";
            graph.writeColorLegendEdgeLabelsByFlags(html);
        }
    }



    // Color legend for edges when colored by coverage.
    if(requestParameters.edgeColoring == "byCoverage" and requestParameters.highlightedOrientedReads.empty()) {
        html <<
            "<h3>Color legend for edges</h3><table>"
            "<tr><td class=left>Coverage";
        for(uint64_t coverage=requestParameters.edgeRedCoverage;
            coverage<= requestParameters.edgeGreenCoverage; coverage++) {
            html << "<td class=centered>";
            if(coverage == requestParameters.edgeRedCoverage) {
                html << "&leq;";
            }
            if(coverage == requestParameters.edgeGreenCoverage) {
                html << "&geq;";
            }
            html << coverage;
        }
        html << "<tr><td class=left>Color";
        for(uint64_t coverage=requestParameters.edgeRedCoverage;
            coverage<= requestParameters.edgeGreenCoverage; coverage++) {
            double h =
                double(coverage - requestParameters.edgeRedCoverage) /
                double(requestParameters.edgeGreenCoverage - requestParameters.edgeRedCoverage);
            h = max(h, 0.);
            h = min(h, 1.);
            const double hue = 120. * h;
            double S, L;
            tie(S, L) = hsvToHsl(1., 0.9);
            html << "<td style='height:15px;background-color:hsl(" << hue << "," <<
                uint64_t(100. * S) << "%," <<
                uint64_t(100. * L) << "%)'>";
        }

        html << "</table>";
    }



    // If there are highlighted oriented reads, write a legend with their colors.
    if(not requestParameters.highlightedOrientedReads.empty()) {
        double S, L;
        tie(S, L) = hsvToHsl(requestParameters.S, requestParameters.V);
        html << "<h3>Color legend for highlighted oriented reads</h3><table><tr>";
        for(const auto& p:requestParameters.highlightedOrientedReads) {
            html << "<td class=centered>" << p.first;
        }
        html << "<tr>";
        for(const auto& p:requestParameters.highlightedOrientedReads) {
            const uint64_t hue = uint64_t(p.second * 360.);
            html << "<td style='height:15px;background-color:hsl(" << hue << "," <<
                uint64_t(100. * S) << "%," <<
                uint64_t(100. * L) << "%)'>";
        }
        html << "</table>";

    }



    // Buttons to resize the svg locally.
    html << "<br>";
    addScaleSvgButtons(html, requestParameters.sizePixels);

    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << "<div id=svgDiv style='display:none'>"; // Make it invisible until after we scale it.
    html << svgFile.rdbuf();
    svgFile.close();

    // Remove the .svg file.
    std::filesystem::remove(svgFileName);

    // Scale to desired size, then make it visible.
    html <<
        "</div>"
        "<script>"
        "var element = document.getElementsByTagName('svg')[0];"
        "w0 = element.getAttribute('width');"
        "h0 = element.getAttribute('height');"
        "element.setAttribute('width', " << requestParameters.sizePixels << ");"
        "w1 = element.getAttribute('width');"
        "h1 = element.getAttribute('height');"
        "element.setAttribute('height', h0 * (w1 / w0));"
        "document.getElementById('svgDiv').setAttribute('style', 'display:block');"
        "</script>";



    // Make the vertices clickable: Ctrl-click recenters
    // the graph at that vertex, right click shows vertex details.
    html << "<script>\n";
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        const LocalMarkerGraphVertex& vertex = graph[v];
        SHASTA_ASSERT(!vertex.markerInfos.empty());
        const string url = requestParameters.urlForVertex(vertex.vertexId);
        html <<
            "element = document.getElementById('vertex" << vertex.vertexId << "');\n"
            "element.onclick = function() {if(!event.ctrlKey) {return;} location.href='" << url <<
            "' + '&sizePixels=' + sizePixels;};\n"
            "element.style.cursor = \"default\";\n";

        // Add a right click to show details.
        const string detailUrl =
            "exploreMarkerGraphVertex?vertexId=" + to_string(vertex.vertexId);
        html <<
            "element.oncontextmenu = function() {window.open('" << detailUrl << "');"
            "return false;};\n";
    }
    html << "</script>\n";



    // Make the edges clickable: Ctrl-click recenters
    // the graph at the source vertex of that edge, right click shows edge details.
    html << "<script>\n";
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph) {
        const LocalMarkerGraphEdge& edge = graph[e];
        const LocalMarkerGraph::vertex_descriptor v0 = source(e, graph);
        const LocalMarkerGraphVertex& vertex0 = graph[v0];
        const string url = requestParameters.urlForVertex(vertex0.vertexId);
        html <<
            "element = document.getElementById('edge" << edge.edgeId << "');\n"
            "element.onclick = function() {if(!event.ctrlKey) {return;} location.href='" << url <<
            "' + '&sizePixels=' + sizePixels;};\n"
            "element.style.cursor = \"default\";\n";

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

    parameters.layoutMethod = "dotLr";
    getParameterValue(request, "layoutMethod", parameters.layoutMethod);

    parameters.minVertexCoverage = 0;
    parameters.minVertexCoverageIsPresent = getParameterValue(
        request, "minVertexCoverage", parameters.minVertexCoverage);

    parameters.minEdgeCoverage = 0;
    parameters.minEdgeCoverageIsPresent = getParameterValue(
        request, "minEdgeCoverage", parameters.minEdgeCoverage);

    string useWeakEdgesString;
    parameters.useWeakEdges = getParameterValue(
        request, "useWeakEdges", useWeakEdgesString);

    string usePrunedEdgesString;
    parameters.usePrunedEdges = getParameterValue(
        request, "usePrunedEdges", usePrunedEdgesString);

    string useSuperBubbleEdgesString;
    parameters.useSuperBubbleEdges = getParameterValue(
        request, "useSuperBubbleEdges", useSuperBubbleEdgesString);

    string useLowCoverageCrossEdgesString;
    parameters.useLowCoverageCrossEdges = getParameterValue(
        request, "useLowCoverageCrossEdges", useLowCoverageCrossEdgesString);

    string useRemovedSecondaryEdgesString;
    parameters.useRemovedSecondaryEdges = getParameterValue(
        request, "useRemovedSecondaryEdges", useRemovedSecondaryEdgesString);

    parameters.sizePixels = 800;
    parameters.sizePixelsIsPresent = getParameterValue(
        request, "sizePixels", parameters.sizePixels);

    parameters.vertexScalingFactor = 1;
    parameters.vertexScalingFactorIsPresent = getParameterValue(
        request, "vertexScalingFactor", parameters.vertexScalingFactor);

    parameters.edgeThicknessScalingFactor = 1;
    parameters.edgeThicknessScalingFactorIsPresent = getParameterValue(
        request, "edgeThicknessScalingFactor", parameters.edgeThicknessScalingFactor);

    parameters.arrowScalingFactor = 1;
    parameters.arrowScalingFactorIsPresent = getParameterValue(
        request, "arrowScalingFactor", parameters.arrowScalingFactor);

    parameters.edgeThickness = "byCoverage";
    getParameterValue(request, "edgeThickness", parameters.edgeThickness);

    parameters.timeout = 30;
    parameters.timeoutIsPresent = getParameterValue(
        request, "timeout", parameters.timeout);

    parameters.vertexLabels = 0;
    getParameterValue(request, "vertexLabels", parameters.vertexLabels);

    parameters.edgeLabels = 0;
    getParameterValue(request, "edgeLabels", parameters.edgeLabels);

    parameters.vertexColoring = "byDistance";
    getParameterValue(request, "vertexColoring", parameters.vertexColoring);

    parameters.edgeColoring = "byFlags";
    getParameterValue(request, "edgeColoring", parameters.edgeColoring);

    parameters.vertexRedCoverage = 1;
    getParameterValue(request, "vertexRedCoverage", parameters.vertexRedCoverage);

    parameters.vertexGreenCoverage = 10;
    getParameterValue(request, "vertexGreenCoverage", parameters.vertexGreenCoverage);

    parameters.vertexRedCoveragePerStrand = 1;
    getParameterValue(request, "vertexRedCoveragePerStrand", parameters.vertexRedCoveragePerStrand);

    parameters.vertexGreenCoveragePerStrand = 2;
    getParameterValue(request, "vertexGreenCoveragePerStrand", parameters.vertexGreenCoveragePerStrand);

    parameters.edgeRedCoverage = 1;
    getParameterValue(request, "edgeRedCoverage", parameters.edgeRedCoverage);

    parameters.edgeGreenCoverage = 10;
    getParameterValue(request, "edgeGreenCoverage", parameters.edgeGreenCoverage);

    parameters.edgeRedCoveragePerStrand = 1;
    getParameterValue(request, "edgeRedCoveragePerStrand", parameters.edgeRedCoveragePerStrand);

    parameters.edgeGreenCoveragePerStrand = 2;
    getParameterValue(request, "edgeGreenCoveragePerStrand", parameters.edgeGreenCoveragePerStrand);

    getParameterValue(request, "highlightedOrientedReads", parameters.highlightedOrientedReadsString);
    parameters.parseHighlightedOrientedReads();

}


// This parses highlightedOrientedReadsString and creates
// highlightedOrientedReads. Each oriented read is assigned a hue
// via hashing of the OrientedReadId. This way, an oriented read
// is always highlighted in the same color.
void LocalMarkerGraphRequestParameters::parseHighlightedOrientedReads()
{
    highlightedOrientedReads.clear();
    if(highlightedOrientedReadsString.empty()) {
        return;
    }

    vector<string> tokens;
    boost::algorithm::split(tokens, highlightedOrientedReadsString, boost::algorithm::is_any_of(" "));

    for(const string& token: tokens) {
        const OrientedReadId orientedReadId = OrientedReadId(token);

        // Hash the OrientedReadId to create a Hue value in (0,1).
        // This way the same OrientedReadId always gets the same color.
        const uint32_t hashValue = MurmurHash2(&orientedReadId, sizeof(OrientedReadId), 751);
        const double hue = double(hashValue) / double(std::numeric_limits<uint32_t>::max());

        highlightedOrientedReads.insert(make_pair(orientedReadId, hue));
    }
}



void LocalMarkerGraphRequestParameters::writeForm(
    ostream& html,
    MarkerGraph::VertexId vertexCount) const
{
    html <<
        "<form>"


        "<h3>Local marker graph creation</h3>"
        "<table>"

        "<tr title='Start vertex id between 0 and " << vertexCount << "'>"
        "<td>Start vertex id"
        "<td class=centered><input type=text required name=vertexId size=8 style='text-align:center'"
        << (vertexIdIsPresent ? ("value='"+to_string(vertexId)+"'") : "") <<
        ">"

        "<tr title='Maximum distance from start vertex (number of edges)'>"
        "<td>Maximum distance"
        "<td class=centered><input type=text required name=maxDistance size=8 style='text-align:center'"
        << (maxDistanceIsPresent ? ("value='" + to_string(maxDistance)+"'") : " value='6'") <<
        ">"

        "<tr>"
        "<td>Minimum vertex coverage"
        "<td class=centered><input type=text required name=minVertexCoverage size=8 style='text-align:center'"
        << (minVertexCoverageIsPresent ? ("value='" + to_string(minVertexCoverage)+"'") : " value='0'") <<
        ">"

        "<tr>"
        "<td>Minimum edge coverage"
        "<td class=centered><input type=text required name=minEdgeCoverage size=8 style='text-align:center'"
        << (minEdgeCoverageIsPresent ? ("value='" + to_string(minEdgeCoverage)+"'") : " value='0'") <<
        ">"

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

        "<tr title='Check to include in the local marker graph "
        "edges that were removed as low coverage cross edges'>"
        "<td>Edges removed as low coverage cross edges"
        "<td class=centered>"
        "<input type=checkbox name=useLowCoverageCrossEdges" <<
        (useLowCoverageCrossEdges ? " checked=checked" : "") << ">"

        "<tr title='Check to include in the local marker graph "
        "secondary edges that were removed during splitting'>"
        "<td>Secondary edges removed during splitting"
        "<td class=centered>"
        "<input type=checkbox name=useRemovedSecondaryEdges" <<
        (useRemovedSecondaryEdges ? " checked=checked" : "") << ">"

        "</table>"
        "<h3>Graphics</h3>"

        "<table>"
        "<tr>"
        "<td colspan=2>Width in pixels"
        "<td class=centered><input type=text required name=sizePixels size=8 style='text-align:center'"
        << (sizePixelsIsPresent ? (" value='" + to_string(sizePixels)+"'") : " value='800'") <<
        ">"

        "<tr>"
        "<td colspan=2>Graph layout method"
        "<td class=left>"
        "<span title='Best for small subgraphs'><input type=radio name=layoutMethod value=dotLr"
        << (layoutMethod=="dotLr" ? " checked=checked" : "") <<
        ">Dot, left to right</span><br>"
        "<span title='Best for small subgraphs with labels'><input type=radio name=layoutMethod value=dotTb"
        << (layoutMethod=="dotTb" ? " checked=checked" : "") <<
        ">Dot, top to bottom</span><br>"
        "<span title='Best for large subgraphs, without labels'><input type=radio name=layoutMethod value=sfdp"
        << (layoutMethod=="sfdp" ? " checked=checked" : "") <<
        ">Sfdp</span>"

        "<tr>"
        "<td colspan=2>Highlight oriented reads"
        "<td class=centered><input type=text name=highlightedOrientedReads size=12"
        << (highlightedOrientedReadsString.empty() ? "" : (" value='" + highlightedOrientedReadsString + "'")) <<
        " title='Enter one or more oriented reads separated by spaces, for example \"432-0 1256-1\"'"
        "</textarea>"

        "<tr title='Maximum time allowed (seconds) for graph creation and layout, or 0 if unlimited'>"
        "<td colspan=2>Timeout (seconds) for graph creation and layout"
        "<td class=centered><input type=text required name=timeout size=8 style='text-align:center'"
        << (timeoutIsPresent ? (" value='" + to_string(timeout)+"'") : " value='30'") <<
        ">"

        "<tr>"
        "<td rowspan=4 class=centered>Vertices"
        "<td>Labels"
        "<td><input type=radio name=vertexLabels value=0" <<
        ((vertexLabels==0) ? " checked=checked" : "") << ">None"
        "<br><input type=radio name=vertexLabels value=1" <<
        ((vertexLabels==1) ? " checked=checked" : "") << ">Terse"
        "<br><input type=radio name=vertexLabels value=2" <<
        ((vertexLabels==2) ? " checked=checked" : "") << ">Verbose"

        "<tr><td>Coloring"
        "<td>"
        "<input type=radio name=vertexColoring value=none"
        << (vertexColoring=="none" ? " checked=checked" : "") <<
        ">None<br>"
        "<input type=radio name=vertexColoring value=byCoverage"
        << (vertexColoring=="byCoverage" ? " checked=checked" : "") <<
        ">By coverage<br>"
        "<input type=radio name=vertexColoring value=byDistance"
        << (vertexColoring=="byDistance" ? " checked=checked" : "") <<
        ">By distance"

        "<tr><td>Color by coverage"
        "<td><table style='margin-left:auto;margin-right:auto'>"
        "<tr><td class=centered>Total<br>coverage<td class=centered>Strand<br>coverage<td class=left>Color"
        "<tr>"
        "<td class=centered><input type=text name=vertexRedCoverage size=4 style='text-align:center'"
        " value='" << vertexRedCoverage << "'>"
        "<td class=centered><input type=text name=vertexRedCoveragePerStrand size=4 style='text-align:center'"
        " value='" << vertexRedCoveragePerStrand << "'>"
        "<td class=centered style='background-color:hsl(0,100%,45%)'>"
        "<tr>"
        "<td class=centered><input type=text name=vertexGreenCoverage size=4 style='text-align:center'" <<
        " value='" << vertexGreenCoverage << "'>"
        "<td class=centered><input type=text name=vertexGreenCoveragePerStrand size=4 style='text-align:center'"
        " value='" << vertexGreenCoveragePerStrand << "'>"
        "<td class=centered style='background-color:hsl(120,100%,45%)'></table>"

        "<tr>"
        "<td>Scaling factor"
        "<td class=centered><input type=text required name=vertexScalingFactor size=8 style='text-align:center'" <<
        " value='" + vertexScalingFactorString() + "'>" <<

        "<tr>"
        "<td rowspan=5 class=centered>Edges"
        "<td>Labels"
        "<td><input type=radio name=edgeLabels value=0" <<
        ((edgeLabels==0) ? " checked=checked" : "") << ">None"
        "<br><input type=radio name=edgeLabels value=1" <<
        ((edgeLabels==1) ? " checked=checked" : "") << ">Terse"
        "<br><input type=radio name=edgeLabels value=2" <<
        ((edgeLabels==2) ? " checked=checked" : "") << ">Verbose"

        "<tr><td>Coloring"
        "<td>"
        "<input type=radio name=edgeColoring value=none"
        << (edgeColoring=="none" ? " checked=checked" : "") <<
        ">None<br>"
        "<input type=radio name=edgeColoring value=byCoverage"
        << (edgeColoring=="byCoverage" ? " checked=checked" : "") <<
        ">By coverage<br>"
        "<input type=radio name=edgeColoring value=byFlags"
        << (edgeColoring=="byFlags" ? " checked=checked" : "") <<
        ">By flags"

        "<tr><td>Color by coverage"
        "<td><table style='margin-left:auto;margin-right:auto'>"
        "<tr><td class=centered>Total<br>coverage<td class=centered>Strand<br>coverage<td class=left>Color"
        "<tr>"
        "<td class=centered><input type=text name=edgeRedCoverage size=4 style='text-align:center'"
        " value='" << edgeRedCoverage << "'>"
        "<td class=centered><input type=text name=edgeRedCoveragePerStrand size=4 style='text-align:center'"
        " value='" << edgeRedCoveragePerStrand << "'>"
        "<td class=centered style='background-color:hsl(0,100%,45%)'>"
        "<tr>"
        "<td class=centered><input type=text name=edgeGreenCoverage size=4 style='text-align:center'" <<
        " value='" << edgeGreenCoverage << "'>"
        "<td class=centered><input type=text name=edgeGreenCoveragePerStrand size=4 style='text-align:center'"
        " value='" << edgeGreenCoveragePerStrand << "'>"
        "<td class=centered style='background-color:hsl(120,100%,45%)'></table>"

        "<tr>"
        "<td>Thickness"
        "<td class=left>"
        "<input type=radio name=edgeThickness value=constant"
        << (edgeThickness=="constant" ? " checked=checked" : "") <<
        ">Constant"
        "<br>"
        "<input type=radio name=edgeThickness value=byCoverage"
        << (edgeThickness=="byCoverage" ? " checked=checked" : "") <<
        ">Proportional to coverage"

        "<tr>"
        "<td >Thickness scaling factor"
        "<td class=centered><input type=text required name=edgeThicknessScalingFactor size=8 style='text-align:center'" <<
        " value='" + edgeThicknessScalingFactorString() + "'>" <<

        "<tr>"
        "<td>Arrow scaling factor"
        "<td class=centered><input type=text required name=arrowScalingFactor size=8 style='text-align:center'" <<
        " value='" + arrowScalingFactorString() + "'>" <<

        "</table>"



        "<br><input type=submit value='Display'>"
        "</form>";
}



bool LocalMarkerGraphRequestParameters::hasMissingRequiredParameters() const
{
    return
        !vertexIdIsPresent ||
        !maxDistanceIsPresent ||
        !timeoutIsPresent;
}



string LocalMarkerGraphRequestParameters::vertexScalingFactorString() const
{
    if(vertexScalingFactorIsPresent) {
        std::ostringstream s;
        s << vertexScalingFactor;
        return s.str();
    } else {
        return "1";
    }
}



string LocalMarkerGraphRequestParameters::arrowScalingFactorString() const
{
    if(arrowScalingFactorIsPresent) {
        std::ostringstream s;
        s << arrowScalingFactor;
        return s.str();
    } else {
        return "1";
    }
}



string LocalMarkerGraphRequestParameters::edgeThicknessScalingFactorString() const
{
    if(edgeThicknessScalingFactorIsPresent) {
        std::ostringstream s;
        s << edgeThicknessScalingFactor;
        return s.str();
    } else {
        return "1";
    }
}



string LocalMarkerGraphRequestParameters::url() const
{
    return
        string("exploreMarkerGraph") +
        "?vertexId=" + to_string(vertexId) +
        "&maxDistance=" + to_string(maxDistance) +
        "&minVertexCoverage=" + to_string(minVertexCoverage) +
        "&minEdgeCoverage=" + to_string(minEdgeCoverage) +
        "&layoutMethod=" + layoutMethod +
        (useWeakEdges ? "&useWeakEdges=on" : "") +
        (usePrunedEdges ? "&usePrunedEdges=on" : "") +
        (useSuperBubbleEdges ? "&useSuperBubbleEdges=on" : "") +
        (useLowCoverageCrossEdges ? "&useLowCoverageCrossEdges=on" : "") +
        "&vertexScalingFactor=" + to_string(vertexScalingFactor) +
        "&edgeThicknessScalingFactor=" + to_string(edgeThicknessScalingFactor) +
        "&arrowScalingFactor=" + to_string(arrowScalingFactor) +
        "&timeout=" + to_string(timeout) +
        "&vertexLabels=" + vertexLabelsString() +
        "&edgeLabels=" + edgeLabelsString() +
        "&vertexColoring=" + vertexColoring +
        "&vertexRedCoverage=" + to_string(vertexRedCoverage) +
        "&vertexGreenCoverage=" + to_string(vertexGreenCoverage) +
        "&edgeColoring=" + edgeColoring +
        "&edgeRedCoverage=" + to_string(edgeRedCoverage) +
        "&edgeGreenCoverage=" + to_string(edgeGreenCoverage) +
        "&highlightedOrientedReads=" + highlightedOrientedReadsString;

}



string LocalMarkerGraphRequestParameters::urlForVertex(uint64_t newVertexId) const
{
    LocalMarkerGraphRequestParameters newParameters = *this;
    newParameters.vertexId = newVertexId;
    return newParameters.url();
}



string LocalMarkerGraphRequestParameters::vertexLabelsString() const
{
    switch(vertexLabels) {
        case 0: return "none";
        case 1: return "terse";
        case 2: return "verbose";
        default: SHASTA_ASSERT(0);
    }
}



string LocalMarkerGraphRequestParameters::edgeLabelsString() const
{
    switch(edgeLabels) {
        case 0: return "none";
        case 1: return "terse";
        case 2: return "verbose";
        default: SHASTA_ASSERT(0);
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
    span<MarkerId> markerIds = markerGraph.getVertexMarkerIds(vertexId);
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
            if(assemblerInfo->readRepresentation == 1) {
                tie(base, repeatCounts[j][i]) =
                    reads->getOrientedReadBaseAndRepeatCount(orientedReadIds[j], uint32_t(marker.position+i));
            } else {
                base = reads->getOrientedReadBase(orientedReadIds[j], uint32_t(marker.position+i));
                repeatCounts[j][i] = 1;
            }
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
    const bool consensusIsAvailable = markerGraph.vertexRepeatCounts.isOpen;
    vector<size_t> consensusRepeatCounts(k);
    if(consensusIsAvailable) {
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
    }




    // Compute concordant and discordant coverage at each position.
    vector<size_t> concordantCoverage(k, 0);
    vector<size_t> discordantCoverage(k, 0);
    if(consensusIsAvailable) {
        for(size_t i=0; i<k; i++) {
            for(size_t j=0; j<markerCount; j++) {
                if(repeatCounts[j][i] == consensusRepeatCounts[i]) {
                    ++concordantCoverage[i];
                } else {
                    ++discordantCoverage[i];
                }
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
        "&useLowCoverageCrossEdges=on"
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


    if(consensusIsAvailable) {
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
    const MarkerGraph::VertexId vertexIdRc = markerGraph.reverseComplementVertex[vertexId];
    html << "<tr><th class=left>Reverse complement vertex<td class=centered>"
        "<a href='exploreMarkerGraphVertex?vertexId=" << vertexIdRc << "'>" <<
        vertexIdRc << "</a> ";


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
            "&amp;showMarkers=on"
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


    if(consensusIsAvailable) {
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
    }



    // Write a row with the marker sequence (run-length).
    html <<
        "<tr><th colspan=2 class=left>Run-length sequence"
        "<td class=centered style='font-family:monospace'>";
    kmer.write(html, assemblerInfo->k);



    // Write rows with coverage information for each represented repeat value.
    if(consensusIsAvailable) {
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

            if(assemblerInfo->readRepresentation == 1) {
                Base base;
                uint8_t repeatCount;
                tie(base, repeatCount) = reads->getOrientedReadBaseAndRepeatCount(orientedReadId, position);
                sequences[j].push_back(base);
                repeatCounts[j].push_back(repeatCount);
            } else {
                const Base base = reads->getOrientedReadBase(orientedReadId, position);
                sequences[j].push_back(base);
                repeatCounts[j].push_back(1);

            }
        }
    }


    // Access stored consensus for this edge.
    const bool consensusIsAvailable =
        markerGraph.edgeConsensus.isOpen() and
        markerGraph.edgeConsensusOverlappingBaseCount.isOpen;
    int storedConsensusOverlappingBaseCount = 0;
    span< pair<Base, uint8_t> > storedConsensus;
    if(consensusIsAvailable)   {
        storedConsensusOverlappingBaseCount = int(markerGraph.edgeConsensusOverlappingBaseCount[edgeId]);
        storedConsensus = markerGraph.edgeConsensus[edgeId];
    }



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
    auto spoaAlignmentEngine = spoa::AlignmentEngine::Create(alignmentType, match, mismatch, gap);
    spoa::Graph spoaAlignmentGraph = {};
    spoaAlignmentGraph.Clear();

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
        "&useLowCoverageCrossEdges=on"
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
        "<tr><th class=left>Coverage<td class=centered>" << markerCount;

    if(assemblerInfo->assemblyMode == 2) {
        html <<
        "<tr><th class=left>Secondary edge?<td class=centered>" <<
        (edge.isSecondary ? "Yes" : "No") <<
        "<tr><th class=left>Removed while splitting of secondary edges?<td class=centered>" <<
        (edge.wasRemovedWhileSplittingSecondaryEdges ? "Yes" : "No");
    } else {
        html <<
        "<tr><th class=left>Removed during transitive reduction?<td class=centered>" <<
        (edge.wasRemovedByTransitiveReduction ? "Yes" : "No") <<
        "<tr><th class=left>Removed during pruning?<td class=centered>" <<
        (edge.wasPruned ? "Yes" : "No") <<
        "<tr><th class=left>Removed during bubble/superbubble removal?<td class=centered>" <<
        (edge.isSuperBubbleEdge ? "Yes" : "No") <<
        "<tr><th class=left>Removed as a low coverage cross edge?<td class=centered>" <<
        (edge.isLowCoverageCrossEdge ? "Yes" : "No");
    }

    // Usage of this edge in the assembly graph.
    if(assemblyGraphPointer and assemblyGraphPointer->markerToAssemblyTable.isOpen()) {
        const auto& locations = assemblyGraphPointer->markerToAssemblyTable[edgeId];
        for(const pair<AssemblyGraph::EdgeId, uint32_t>& location: locations) {
            const AssemblyGraph::EdgeId assemblyGraphEdgeId = location.first;
            const uint32_t positionInAssemblyGraphEdge = location.second;
            html << "<tr><th class=left>In assembly graph edge<td class=centered>" <<
                assemblyGraphEdgeId << " at position " <<
                positionInAssemblyGraphEdge;
        }
    }


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

    const MarkerGraph::EdgeId edgeIdRc = markerGraph.reverseComplementEdge[edgeId];
    html << "<tr><th class=left>Reverse complement edge<td class=centered>"
        "<a href='exploreMarkerGraphEdge?edgeId=" << edgeIdRc << "'>" <<
        edgeIdRc << "</a> ";

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
                "&amp;showMarkers=on"
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
        " title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand0", strand0IsPresent && strand0==0, strand0IsPresent && strand0==1);
    html <<
        "<br><input type=text name=readId1 required size=8 " <<
        (readId1IsPresent ? "value="+to_string(readId1) : "") <<
        " title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"
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
        " size=8 title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"
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



// This shows a table that follows a reads and it alignments in the marker graph.
void Assembler::followReadInMarkerGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the request parameters.
    ReadId readId0 = 0;
    const bool readId0IsPresent = getParameterValue(request, "readId", readId0);
    Strand strand0 = 0;
    const bool strand0IsPresent = getParameterValue(request, "strand", strand0);
    uint32_t firstOrdinal = 0;
    const bool firstOrdinalIsPresent = getParameterValue(request, "firstOrdinal", firstOrdinal);
    uint32_t lastOrdinal = 0;
    const bool lastOrdinalIsPresent = getParameterValue(request, "lastOrdinal", lastOrdinal);
    string whichAlignments = "ReadGraphAlignments";
    getParameterValue(request, "whichAlignments", whichAlignments);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Follow this oriented read and its alignments in the marker graph'> "
        "<br>Read <input type=text name=readId required" <<
        (readId0IsPresent ? (" value=" + to_string(readId0)) : "") <<
        " size=8 title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand", strand0IsPresent && strand0==0, strand0IsPresent && strand0==1);
    html <<
        " First marker ordinal: <input type=text name=firstOrdinal required" <<
        (firstOrdinalIsPresent ? (" value=" + to_string(firstOrdinal)) : "") << ">"
        " Last marker ordinal: <input type=text name=lastOrdinal required" <<
        (lastOrdinalIsPresent ? (" value=" + to_string(lastOrdinal)) : "") << ">";
    html << "<br><input type=radio name=whichAlignments value=AllAlignments" <<
        (whichAlignments=="AllAlignments" ? " checked=checked" : "") << "> All alignments";
    html << "<br><input type=radio name=whichAlignments value=ReadGraphAlignments" <<
        (whichAlignments=="ReadGraphAlignments" ? " checked=checked" : "") <<
        "> Only alignments used in the read graph.";
    html << "</form>";

    // If a required parameter is missing, stop here.
    if(!readId0IsPresent || !strand0IsPresent || !firstOrdinalIsPresent || !lastOrdinalIsPresent) {
        return;
    }
    const OrientedReadId orientedReadId0(readId0, strand0);
    const uint32_t markerCount0 = uint32_t(markers.size(orientedReadId0.getValue()));
    const uint32_t ordinal0Begin = firstOrdinal;
    const uint32_t ordinal0End = min(markerCount0, lastOrdinal + 1);




    // Loop over alignment involving this oriented read, as stored in the
    // alignment table.
    Alignment alignment;
    vector<OrientedReadId> orientedReadIds1;
    vector<bool> isInReadGraph;
    vector< vector<uint32_t> > alignedOrdinals1Matrix; // alignedOrdinals1Matrix[i][ordinal0] = ordinal1;
    const auto alignmentTable0 = alignmentTable[orientedReadId0.getValue()];
    for(const auto alignmentId: alignmentTable0) {
        const AlignmentData& ad = alignmentData[alignmentId];

        // If this alignment is not in the read graph and only read graph alignments
        // were requested, skip it.
        if((whichAlignments=="ReadGraphAlignments") and (not ad.info.isInReadGraph)) {
            continue;
        }

        // The alignment is stored with its first read on strand 0.
        OrientedReadId alignmentOrientedReadId0(ad.readIds[0], 0);
        OrientedReadId alignmentOrientedReadId1(ad.readIds[1],
            ad.isSameStrand ? 0 : 1);

        // Access the alignment and decompress it.
        const span<char> compressedAlignment = compressedAlignments[alignmentId];
        const span<const char> constCompressedAlignment(compressedAlignment.begin(), compressedAlignment.end());
        decompress(constCompressedAlignment, alignment);
        SHASTA_ASSERT(alignment.ordinals.size() == ad.info.markerCount);

        // Swap the reads, if necessary.
        bool swapReads = false;
        if(alignmentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
            swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
            swapReads = true;
        }

        // Reverse complement, if necessary.
        bool reverseComplement = false;
        if(alignmentOrientedReadId0 != orientedReadId0) {
            alignmentOrientedReadId0.flipStrand();
            alignmentOrientedReadId1.flipStrand();
            reverseComplement = true;
        }
        SHASTA_ASSERT(alignmentOrientedReadId0 == orientedReadId0);
        const OrientedReadId orientedReadId1 = alignmentOrientedReadId1;
        const uint32_t markerCount1 = uint32_t(markers.size(orientedReadId1.getValue()));
        orientedReadIds1.push_back(orientedReadId1);
        isInReadGraph.push_back(ad.info.isInReadGraph);

        // Store aligned ordinals of orientedReadId1.
        alignedOrdinals1Matrix.resize(orientedReadIds1.size());
        vector<uint32_t>& alignedOrdinals1 = alignedOrdinals1Matrix.back();
        alignedOrdinals1.resize(markerCount0, std::numeric_limits<uint32_t>::max());
        for(const auto& ordinals: alignment.ordinals) {
            uint32_t ordinal0 = ordinals[0];
            uint32_t ordinal1 = ordinals[1];
            if(swapReads) {
                swap(ordinal0, ordinal1);
            }
            if(reverseComplement) {
                ordinal0 = markerCount0 - 1 - ordinal0;
                ordinal1 = markerCount1 - 1 - ordinal1;
            }
            SHASTA_ASSERT(alignedOrdinals1[ordinal0] == std::numeric_limits<uint32_t>::max());
            alignedOrdinals1[ordinal0] = ordinal1;
        }


    }



    // Write the page header.
    html << "<h1>Follow oriented read " << orientedReadId0 <<
        " and its alignments in the marker graph</h1>"
        "<p>This follows oriented read " << orientedReadId0 <<
        " and its alignments in the marker graph.";
    if(whichAlignments=="ReadGraphAlignments") {
        html << " You selected to only display alignments that "
            "correspond to a read graph edge.";
    } else {
        html << " Alignments that correspond to a read graph edge are displayed "
            "with a light blue background.";
    }


    // Write the table header.
    html << "<table><tr><th colspan=2 style='background-color:Beige'>" << orientedReadId0;
    for(uint64_t i=0; i<orientedReadIds1.size(); i++) {
        html << "<th";
        if(isInReadGraph[i]) {
            html << " style='background-color:LightCyan'";
        }
        html << " colspan=2>" << orientedReadIds1[i];
    }
    html << "<tr>"
        "<th style='background-color:Beige'>Marker<br>ordinal"
        "<th style='background-color:Beige'>Marker<br>graph<br>vertex";
    for(uint64_t i=0; i<orientedReadIds1.size(); i++) {
        html << "<th";
        if(isInReadGraph[i]) {
            html << " style='background-color:LightCyan'";
        }
        html << ">";
        html << "Marker<br>ordinal<th";
        if(isInReadGraph[i]) {
            html << " style='background-color:LightCyan'";
        }
        html << ">Marker<br>graph<br>vertex";
    }



    // Write a table row for each marker in orientedReadId0.
    for(uint32_t ordinal0=ordinal0Begin; ordinal0<ordinal0End; ordinal0++) {
        const MarkerId markerId0 = getMarkerId(orientedReadId0, ordinal0);
        const MarkerGraph::CompressedVertexId vertexId0 = markerGraph.vertexTable[markerId0];

        // Write ordinal and marker graph vertex for orientedReadId0.
        html << "<tr><td class=centered title='" << orientedReadId0 << " ordinal'"
            " style='background-color:Beige'>" << ordinal0 <<
            "<td class=centered title='" << orientedReadId0 << " marker graph vertex'";
        if(vertexId0 == MarkerGraph::invalidCompressedVertexId) {
            html << " style='background-color:Beige'>";
        }
        else {
            html << " style='background-color:Aquamarine'>" <<
                "<a href='exploreMarkerGraphVertex?vertexId=" << vertexId0 << "'>" << vertexId0 << "</a>";
        }


        // Loop over aligned reads.
        for(uint64_t i=0; i<orientedReadIds1.size(); i++) {
            const OrientedReadId orientedReadId1 = orientedReadIds1[i];
            const vector<uint32_t>& alignedOrdinals1 = alignedOrdinals1Matrix[i];

            // Marker ordinal.
            html << "<td class=centered title='" << orientedReadId1 << " ordinal'";
            if(isInReadGraph[i]) {
                html << " style='background-color:LightCyan'";
            }
            html << ">";
            const uint32_t ordinal1 = alignedOrdinals1[ordinal0];
            if(ordinal1 != std::numeric_limits<uint32_t>::max()) {
                html << ordinal1;
            }



            // Marker graph vertex.
            if(ordinal1 == std::numeric_limits<uint32_t>::max()) {

                // There is no aligned ordinal.
                html << "<td class=centered title='" << orientedReadId1 << " marker graph vertex'";
                if(isInReadGraph[i]) {
                    html << " style='background-color:LightCyan'";
                }
                html << ">";

            } else {

                // There is an aligned ordinal. Look for a marker graph vertex.
                const MarkerId markerId1 = getMarkerId(orientedReadId1, ordinal1);
                const MarkerGraph::CompressedVertexId vertexId1 = markerGraph.vertexTable[markerId1];
                html << "<td class=centered title='" << orientedReadId1 << " marker graph vertex'";
                if(vertexId1 == MarkerGraph::invalidCompressedVertexId) {

                    // There is no marker graph vertex.
                    if(isInReadGraph[i]) {
                        html << " style='background-color:LightCyan'";
                    }
                    html << ">";

                } else {

                    // There is a marker graph vertex.
                    if(vertexId0 != MarkerGraph::invalidCompressedVertexId) {
                        // Color based on whether or not it is the same as vertexId0.
                        if(vertexId1 == vertexId0) {
                            html << "style='background-color:Aquamarine'";
                        } else {
                            html << "style='background-color:LightPink'";
                        }
                    } else {
                        if(isInReadGraph[i]) {
                            html << " style='background-color:LightCyan'";
                        }
                    }
                    html << ">" <<
                        "<a href='exploreMarkerGraphVertex?vertexId=" << vertexId1 << "'>" << vertexId1 << "</a>";
                }
            }
        }
    }

    // Finish the table.
    html << "</table>";

}



void Assembler::exploreMarkerConnectivity(
    const vector<string>& request,
    ostream& html)
{
    // Get the request parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);
    uint32_t ordinal = 0;
    const bool ordinalIsPresent = getParameterValue(request, "ordinal", ordinal);
    string whichAlignments = "ReadGraphAlignments";
    getParameterValue(request, "whichAlignments", whichAlignments);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Explore connectivity of this marker'> "
        "<br>Read <input type=text name=readId required" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " size=8 title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);
    html <<
        "<br>Marker ordinal <input type=text name=ordinal required";
    if(ordinalIsPresent) {
        html << " value=" << ordinal;
    }
    html <<
        ">"
        "<br><input type=radio name=whichAlignments value=AllAlignments" <<
        (whichAlignments=="AllAlignments" ? " checked=checked" : "") << "> Use all alignments";
    html << "<br><input type=radio name=whichAlignments value=ReadGraphAlignments" <<
        (whichAlignments=="ReadGraphAlignments" ? " checked=checked" : "") <<
        "> Only use alignments in the read graph.";
    html << "</form>";
    const bool useReadGraphAlignmentsOnly = (whichAlignments == "ReadGraphAlignments");

    // If the required parameters are missing, stop here.
    if(not(readIdIsPresent and strandIsPresent and ordinalIsPresent)) {
        return;
    }
    const OrientedReadId orientedReadId(readId, strand);

    // Check the ordinal.
    const uint64_t markerCount = markers.size(orientedReadId.getValue());
    if(ordinal >= markerCount) {
        html << "<p>" << orientedReadId << " has " << markerCount << " markers.";
        return;
    }



    // Create an undirected graph in which each vertex represents a marker
    // and aligned markers are joined by an edge.
    MarkerConnectivityGraph graph;
    createMarkerConnectivityGraph(orientedReadId, ordinal, useReadGraphAlignmentsOnly, graph);



    // Count how many times each oriented read appears.
    std::map<OrientedReadId, uint64_t> frequencyMap;
    BGL_FORALL_VERTICES(v, graph, MarkerConnectivityGraph) {
        const OrientedReadId orientedReadId = graph[v].first;
        ++frequencyMap[orientedReadId];
    }


    // Write the graph out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    ofstream dotFile(dotFileName);
    dotFile << "graph MarkerConnectivity {\n";
    BGL_FORALL_VERTICES(v, graph, MarkerConnectivityGraph) {
        const MarkerDescriptor markerDescriptor = graph[v];
        const OrientedReadId orientedReadId1 = markerDescriptor.first;
        const uint32_t ordinal1 = markerDescriptor.second;
        dotFile << "\"" << orientedReadId1 << "-" << ordinal1 << "\""
            " [label=\"" << orientedReadId1 << "\\n" << ordinal1 <<
            "\"";
        if(frequencyMap[orientedReadId1] != 1) {
            dotFile << " style=filled fillcolor=pink";
        } else {
            dotFile << " style=filled fillcolor=cornsilk";
        }
        dotFile << "];\n";
    }
    BGL_FORALL_EDGES(e, graph, MarkerConnectivityGraph) {
        const auto v0 = source(e, graph);
        const auto v1 = target(e, graph);
        const MarkerDescriptor markerDescriptor0 = graph[v0];
        const MarkerDescriptor markerDescriptor1 = graph[v1];
        dotFile
            << "\"" << markerDescriptor0.first << "-" << markerDescriptor0.second << "\"--"
            << "\"" << markerDescriptor1.first << "-" << markerDescriptor1.second << "\";\n";
    }
    dotFile << "}\n";
    dotFile.close();



    // Use graphviz to render it to svg.
    const string command = timeoutCommand() + " 30 sfdp -O -T svg " + dotFileName +
        " -Goverlap=false -Gsplines=true -Gsmoothing=triangle";
    const int commandStatus = ::system(command.c_str());
    if(WIFEXITED(commandStatus)) {
        const int exitStatus = WEXITSTATUS(commandStatus);
        if(exitStatus == 124) {
            html << "<p>Timeout for graph layout exceeded.";
            std::filesystem::remove(dotFileName);
            return;
        }
        else if(exitStatus!=0) {
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
    std::filesystem::remove(dotFileName);

    // Buttons to resize the svg locally.
    const int sizePixels = 800;
    addScaleSvgButtons(html, sizePixels);
    html << "<br>Found " << num_vertices(graph) << " markers.";

    // Display the svg file.
    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << "<div id=svgDiv style='display:none'>"; // Make it invisible until after we scale it.
    html << svgFile.rdbuf();
    svgFile.close();

    // Scale to desired size, then make it visible.
    html <<
        "</div>"
        "<script>"
        "var svgElement = document.getElementsByTagName('svg')[0];"
        "svgElement.setAttribute('width', " << sizePixels << ");"
        "document.getElementById('svgDiv').setAttribute('style', 'display:block');"
        "</script>";

    // Remove the .svg file.
    std::filesystem::remove(svgFileName);
}

