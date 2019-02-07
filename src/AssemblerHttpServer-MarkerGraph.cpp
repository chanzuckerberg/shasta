// Shasta.
#include "Assembler.hpp"
#include "ConsensusCaller.hpp"
#include "LocalMarkerGraph.hpp"
using namespace ChanZuckerberg;
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



void Assembler::exploreMarkerGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the request parameters.
    LocalMarkerGraphRequestParameters requestParameters;
    getLocalMarkerGraphRequestParameters(request, requestParameters);

    // Write the form.
    requestParameters.writeForm(html, globalMarkerGraphVertices.size());

    // If any required values are missing, stop here.
    if(requestParameters.hasMissingRequiredParameters()) {
        return;
    }



    // Validity checks.
    if(requestParameters.vertexId > globalMarkerGraphVertices.size()) {
        html << "<p>Invalid vertex id " << requestParameters.vertexId;
        html << ". Must be between 0 and " << globalMarkerGraphVertices.size()-1 << " inclusive.";
        return;
    }



    // Create the local marker graph.
    LocalMarkerGraph graph(
        uint32_t(assemblerInfo->k),
        reads,
        assemblerInfo->useRunLengthReads,
        readRepeatCounts,
        markers,
        globalMarkerGraphVertex,
        *consensusCaller);
    const auto createStartTime = steady_clock::now();
    if(requestParameters.useStoredConnectivity) {
        if(!extractLocalMarkerGraphUsingStoredConnectivity(
            requestParameters.vertexId,
            requestParameters.maxDistance,
            requestParameters.timeout,
            requestParameters.useWeakEdges,
            requestParameters.usePrunedEdges,
            requestParameters.useShortCycleEdges,
            requestParameters.useBubbleEdges,
            requestParameters.useBubbleReplacementEdges,
            requestParameters.useSuperBubbleEdges,
            requestParameters.useSuperBubbleReplacementEdges,
            graph)) {
            html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
            return;
        }
    } else {
        if(!extractLocalMarkerGraph(
            requestParameters.vertexId,
            requestParameters.maxDistance,
            requestParameters.timeout,
            graph)) {
            html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
            return;
        }
    }
    if(num_vertices(graph) == 0) {
        html << "<p>The local marker graph is empty.";
        return;
    }
    graph.approximateTopologicalSort();
    vector< pair<shasta::Base, int> > sequence;
    if( requestParameters.showOptimalSpanningTree ||
        requestParameters.portionToDisplay=="spanningTree" ||
        requestParameters.portionToDisplay=="optimalPath" ||
        requestParameters.portionToDisplay=="clippedOptimalPath" ||
        requestParameters.showAssembledSequence) {
        graph.computeOptimalSpanningTree();
        graph.computeOptimalSpanningTreeBestPath();
        if(requestParameters.showAssembledSequence) {
            graph.assembleDominantSequence(requestParameters.maxDistance, sequence);
        }
    }
    const auto createFinishTime = steady_clock::now();
    if(seconds(createFinishTime - createStartTime) > requestParameters.timeout) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }

    html <<
        "<script>\n"
        "function positionAtVertex(vertexId) {\n";
    if(requestParameters.detailed) {
        html <<
            "var element = document.getElementById('a_vertexDistance' + vertexId);\n";
    } else {
        html <<
            "var element = document.getElementById('vertex' + vertexId);\n";
    }
    html <<
        "var r = element.getBoundingClientRect();\n"
        "window.scrollBy((r.left + r.right - window.innerWidth) / 2, (r.top + r.bottom - window.innerHeight) / 2);\n"
        "}\n"
        "</script>\n";



    // Assembled sequence from the graph, if requested.
    // This must be done while the graph is still intact.
    if(requestParameters.showAssembledSequence) {

        const string fastaSequenceName =
            "AssembledSequence-" + to_string(requestParameters.vertexId) +
            "-" + to_string(requestParameters.maxDistance);
        const string fastaFileName = fastaSequenceName + ".fa";
        string fastaString = ">" + fastaSequenceName + " length " + to_string(sequence.size()) + "\\n";
        for(const auto& p: sequence) {
            char c = p.first.character();
            const auto coverage = p.second;
            if(coverage < 10) {
                c = char(std::tolower(c));
            }
            fastaString += c;
        }
        fastaString += "\\n";

        html << "<h2>Sequence assembled from this local marker graph</h2>";


        if(!assemblerInfo->useRunLengthReads) {
            html <<
                "<h4 style='margin:0'>"
                "Assembled sequence";;
            if(assemblerInfo->useRunLengthReads) {
                html << " in run-length representation";
            }
            html << " (" << sequence.size() << " bases).</h4>";
            html << "<br>Coverage/consensus ("
                "blank if &ge; 10) is shown in blue. Sequence base is lower case if coverage &lt; 10. "
                "<a id=fastaDownload>Download in FASTA format</a>"
                "<script>"
                "var element = document.getElementById('fastaDownload');"
                "element.setAttribute('href', 'data:text/plain;charset=utf-8,' +"
                "encodeURIComponent('" << fastaString << "'));"
                "element.setAttribute('download', '" << fastaFileName << "');"
                "</script>";

            // Labels.
            if(assemblerInfo->useRunLengthReads) {
                html << "<pre  style='margin:0' title='Position on run-length representation of assembled sequence'>";
            } else {
                html << "<pre  style='margin:0' title='Position on assembled sequence'>";
            }
            for(size_t position=0; position<sequence.size(); position+=10) {
                const string label = to_string(position);
                html << label;
                for(size_t i=0; i<10-label.size(); i++) {
                    html << " ";
                }
            }
            html << "</pre>";

            // Scale.
            if(assemblerInfo->useRunLengthReads) {
                html << "<pre  style='margin:0' title='Position on run-length representation of assembled sequence'>";
            } else {
                html << "<pre  style='margin:0' title='Position on assembled sequence'>";
            }
            for(size_t position=0; position<sequence.size(); position++) {
                if((position%10)==0) {
                    html << "|";
                } else if((position%5)==0) {
                    html << "+";
                } else {
                    html << ".";
                }
            }
            html << "</pre>";

            // Sequence.
            if(assemblerInfo->useRunLengthReads) {
                html << "<pre  style='margin:0' title='Run-length representation of assembled sequence'>";
            } else {
                html << "<pre  style='margin:0' title='Assembled sequence'>";
            }
            for(const auto& p: sequence) {
                const auto coverage = p.second;
                if(coverage < 10) {
                    html << char(std::tolower(p.first.character()));
                } else {
                    html << p.first;
                }
            }
            html << "</pre>";

            // Coverage.
            if(assemblerInfo->useRunLengthReads) {
                html << "<pre  style='margin:0;color:blue' title='Coverage/consensus for run-length representation of assembled sequence'>";
            } else {
                html << "<pre  style='margin:0;color:blue' title='Assembled sequence coverage/consensus'>";
            }
            for(const auto& p: sequence) {
                const auto coverage = p.second;
                if(coverage<10) {
                    html << coverage;
                } else {
                    html << " ";
                }
            }
            html << "</pre>";
        }



        if(assemblerInfo->useRunLengthReads) {
            graph.assembleDominantSequenceUsingSeqan(html);
        }

#if 0
        // Also show alignments of oriented reads to assembled sequence.
        html << "<br><h4 style='margin:0'>Marker alignments of assembled sequence to oriented reads</h4>";
        if(assemblerInfo->useRunLengthReads) {
            html << "<br>All read and assembled sequence in this table is "
                "in run-length representation.";
            showLocalMarkerGraphAlignments(html, graph, requestParameters);
        }
#endif
    }



    // Remove the portion of the graph that we don't want to display.
    if(requestParameters.portionToDisplay=="spanningTree") {
        graph.removeNonSpanningTreeEdges();
    } else if(requestParameters.portionToDisplay=="optimalPath") {
        graph.removeAllExceptOptimalPath();
    } else if(requestParameters.portionToDisplay=="clippedOptimalPath") {
        graph.removeAllExceptClippedOptimalPath();
    } else if(requestParameters.portionToDisplay=="none") {
        graph.clear();
    }

    // Write it out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = "/dev/shm/" + uuid + ".dot";
    graph.write(
        dotFileName,
        requestParameters.minCoverage,
        requestParameters.maxDistance,
        requestParameters.detailed,
        requestParameters.showVertexId);

    // Compute layout in svg format.
    const string command =
        "timeout " + to_string(requestParameters.timeout - seconds(createFinishTime - createStartTime)) +
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
            filesystem::remove(dotFileName);
            throw runtime_error("Error " + to_string(exitStatus) + " running graph layout command: " + command);
        }
    } else if(WIFSIGNALED(commandStatus)) {
        const int signalNumber = WTERMSIG(commandStatus);
        throw runtime_error("Signal " + to_string(signalNumber) + " while running graph layout command: " + command);
    } else {
        throw runtime_error("Abnormal status " + to_string(commandStatus) + " while running graph layout command: " + command);

    }
    // Remove the .dot file.
    filesystem::remove(dotFileName);


    // Write the graph.
    if(requestParameters.portionToDisplay == "none") {
        html << "<p>As requested, the graph was not displayed.";
    } else {
        // Write the page title.
        const string legendName =
            requestParameters.detailed ?
            "MarkerGraphLegend-Detailed.html" :
            "MarkerGraphLegend-Compact.html";
        html <<
            "<h2>Marker graph near marker graph vertex " << requestParameters.vertexId <<
            " <a href='docs/" << legendName << "'>(see legend)</a></h2>";

        const string svgFileName = dotFileName + ".svg";
        ifstream svgFile(svgFileName);
        html << svgFile.rdbuf();
        svgFile.close();

        // Remove the .svg file.
        filesystem::remove(svgFileName);

        // Make the vertices clickable to recompute the graph with the
        // same parameters, but starting at the clicked vertex.
        // For a detailed graph, only the "Distance" label of each vertex
        // is made clickable.
        html << "<script>\n";
        BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
            const LocalMarkerGraphVertex& vertex = graph[v];
            CZI_ASSERT(!vertex.markerInfos.empty());
            const auto& markerInfo = vertex.markerInfos.front();
            const string url =
                "exploreMarkerGraph?readId=" + to_string(markerInfo.orientedReadId.getReadId()) +
                "&strand=" + to_string(markerInfo.orientedReadId.getStrand()) +
                "&ordinal="  + to_string(markerInfo.ordinal) +
                "&maxDistance=" + to_string(requestParameters.maxDistance) +
                "&minCoverage=" + to_string(requestParameters.minCoverage) +
                "&sizePixels=" + to_string(requestParameters.sizePixels) +
                "&timeout=" + to_string(requestParameters.timeout) +
                "&portionToDisplay=" + requestParameters.portionToDisplay +
                (requestParameters.detailed ? "&detailed=on" : "") +
                (requestParameters.showVertexId ? "&showVertexId=on" : "") +
                (requestParameters.showOptimalSpanningTree ? "&showOptimalSpanningTree=on" : "");
            if(requestParameters.detailed) {
                html <<
                    "document.getElementById('a_vertexDistance" << vertex.vertexId <<
                    "').onclick = function() {location.href='" << url << "';};\n";
            } else {
                html <<
                    "document.getElementById('vertex" << vertex.vertexId <<
                    "').onclick = function() {location.href='" << url << "';};\n";

                // We are displaying the graph in compact mode.
                // Add a right click to recenter and show detailed.
                const string detailUrl =
                    "exploreMarkerGraph?readId=" + to_string(markerInfo.orientedReadId.getReadId()) +
                    "&strand=" + to_string(markerInfo.orientedReadId.getStrand()) +
                    "&ordinal="  + to_string(markerInfo.ordinal) +
                    "&maxDistance=1" +
                    "&minCoverage=" + to_string(requestParameters.minCoverage) +
                    "&sizePixels=" + to_string(requestParameters.sizePixels) +
                    "&timeout=" + to_string(requestParameters.timeout) +
                    "&portionToDisplay=" + requestParameters.portionToDisplay +
                    "&detailed=on" +
                    (requestParameters.showVertexId ? "&showVertexId=on" : "") +
                    (requestParameters.showOptimalSpanningTree ? "&showOptimalSpanningTree=on" : "");
                html <<
                    "document.getElementById('vertex" << vertex.vertexId <<
                    "').oncontextmenu = function() {location.href='" << detailUrl << "';"
                    "return false;};\n";
            }
        }
        html << "</script>\n";



        // Position the start vertex at the center of the window.
        html <<
            "<script>\n"
            "positionAtVertex(" << requestParameters.vertexId << ");\n"
            "</script>\n";
    }
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

    string detailedString;
    parameters.detailed = getParameterValue(
        request, "detailed", detailedString);

    string useStoredConnectivityString;
    parameters.useStoredConnectivity = getParameterValue(
        request, "useStoredConnectivity", useStoredConnectivityString);

    string useWeakEdgesString;
    parameters.useWeakEdges = getParameterValue(
        request, "useWeakEdges", useWeakEdgesString);

    string usePrunedEdgesString;
    parameters.usePrunedEdges = getParameterValue(
        request, "usePrunedEdges", usePrunedEdgesString);

    string useShortCycleEdgesString;
    parameters.useShortCycleEdges = getParameterValue(
        request, "useShortCycleEdges", useShortCycleEdgesString);

    string useBubbleEdgesString;
    parameters.useBubbleEdges = getParameterValue(
        request, "useBubbleEdges", useBubbleEdgesString);

    string useBubbleReplacementEdgesString;
    parameters.useBubbleReplacementEdges = getParameterValue(
        request, "useBubbleReplacementEdges", useBubbleReplacementEdgesString);

    string useSuperBubbleEdgesString;
    parameters.useSuperBubbleEdges = getParameterValue(
        request, "useSuperBubbleEdges", useSuperBubbleEdgesString);

    string useSuperBubbleReplacementEdgesString;
    parameters.useSuperBubbleReplacementEdges = getParameterValue(
        request, "useSuperBubbleReplacementEdges", useSuperBubbleReplacementEdgesString);

    string showVertexIdString;
    parameters.showVertexId = getParameterValue(
        request, "showVertexId", showVertexIdString);

    string showOptimalSpanningTreeString;
    parameters.showOptimalSpanningTree = getParameterValue(
        request, "showOptimalSpanningTree", showOptimalSpanningTreeString);

    string showAssembledSequenceString;
    parameters.showAssembledSequence = getParameterValue(
        request, "showAssembledSequence", showAssembledSequenceString);

    parameters.minCoverage = 0;
    parameters.minCoverageIsPresent = getParameterValue(
        request, "minCoverage", parameters.minCoverage);

    parameters.sizePixels = 800;
    parameters.sizePixelsIsPresent = getParameterValue(
        request, "sizePixels", parameters.sizePixels);

    parameters.timeout = 30;
    parameters.timeoutIsPresent = getParameterValue(
        request, "timeout", parameters.timeout);

    parameters.portionToDisplay = "all";
    getParameterValue(
        request, "portionToDisplay", parameters.portionToDisplay);

}



void Assembler::LocalMarkerGraphRequestParameters::writeForm(
    ostream& html,
    GlobalMarkerGraphVertexId vertexCount) const
{
    html <<
        "<h3>Display a local subgraph of the global marker graph</h3>"
        "<form>"

        "<table>"

        "<tr title='Vertex id between 0 and " << vertexCount << "'>"
        "<td>Vertex id"
        "<td><input type=text required name=vertexId size=8 style='text-align:center'"
        << (vertexIdIsPresent ? ("value='"+to_string(vertexId)+"'") : "") <<
        ">"

        "<tr title='Maximum distance from start vertex (number of edges)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        << (maxDistanceIsPresent ? ("value='" + to_string(maxDistance)+"'") : " value='6'") <<
        ">"

        "<tr title='Check for detailed graph with labels'>"
        "<td>Detailed"
        "<td class=centered><input type=checkbox name=detailed"
        << (detailed ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to use stored connectivity of the global marker graph "
        "instead of constructing it on the fly'>"
        "<td>Use stored connectivity"
        "<td class=centered><input type=checkbox name=useStoredConnectivity"
        << (useStoredConnectivity ? " checked=checked" : "") <<
        "><td class=left>"
        "Options only used when \"Use stored connectivity\" is checked:<br>"
        "<input type=checkbox name=useWeakEdges" << (useWeakEdges ? " checked=checked" : "") << ">Weak edges"
        "<br><input type=checkbox name=usePrunedEdges" << (usePrunedEdges ? " checked=checked" : "") << ">Pruned edges"
        "<br><input type=checkbox name=useShortCycleEdges" << (useShortCycleEdges ? " checked=checked" : "") << ">Short cycle edges"
        "<br><input type=checkbox name=useBubbleEdges" << (useBubbleEdges ? " checked=checked" : "") << ">Bubble edges"
        "<br><input type=checkbox name=useBubbleReplacementEdges" << (useBubbleReplacementEdges ? " checked=checked" : "") << ">Bubble replacement edges"
        "<br><input type=checkbox name=useSuperBubbleEdges" << (useSuperBubbleEdges ? " checked=checked" : "") << ">Superbubble edges"
        "<br><input type=checkbox name=useSuperBubbleReplacementEdges" << (useSuperBubbleReplacementEdges ? " checked=checked" : "") << ">Superbubble replacement edges"

        "<tr title='Check to show vertex ids (only useful for debugging)'>"
        "<td>Show vertex ids"
        "<td class=centered><input type=checkbox name=showVertexId"
        << (showVertexId ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to show in purple the optimal spanning tree of the local marker graph"
        " (always shown if assembled sequence is also selected)'>"
        "<td>Show optimal spanning tree"
        "<td class=centered><input type=checkbox name=showOptimalSpanningTree"
        << (showOptimalSpanningTree ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to show sequence assembled from this local marker graph'>"
        "<td>Show assembled sequence"
        "<td class=centered><input type=checkbox name=showAssembledSequence"
        << (showAssembledSequence ? " checked=checked" : "") <<
        ">"


        "<tr title='Minimum coverage (number of markers) for a vertex or edge to be considered strong. "
        "Affects the coloring of vertices and edges.'>"
        "<td>Coverage threshold"
        "<td><input type=text required name=minCoverage size=8 style='text-align:center'"
        << (minCoverageIsPresent ? (" value='" + to_string(minCoverage)+"'") : " value='3'") <<
        ">"

        "<tr title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "<td>Graphics size in pixels"
        "<td><input type=text required name=sizePixels size=8 style='text-align:center'"
        << (sizePixelsIsPresent ? (" value='" + to_string(sizePixels)+"'") : " value='800'") <<
        ">"

        "<tr>"
        "<td>Timeout (seconds) for graph creation and layout"
        "<td><input type=text required name=timeout size=8 style='text-align:center'"
        << (timeoutIsPresent ? (" value='" + to_string(timeout)+"'") : " value='30'") <<
        ">"
        "</table>"



        // Radio buttons to choose what portion of the graph to display.
        "Portion of the graph to display:"

        "<br><input type=radio name=portionToDisplay value=all" <<
        (portionToDisplay=="all" ? " checked=on" : "") <<
        ">Entire graph"

        "<br><input type=radio name=portionToDisplay value=spanningTree" <<
        (portionToDisplay=="spanningTree" ? " checked=on" : "") <<
        ">Optimal spanning tree only"

        "<br><input type=radio name=portionToDisplay value=optimalPath" <<
        (portionToDisplay=="optimalPath" ? " checked=on" : "") <<
        ">Optimal spanning tree best path only"

        "<br><input type=radio name=portionToDisplay value=clippedOptimalPath" <<
        (portionToDisplay=="clippedOptimalPath" ? " checked=on" : "") <<
        "><span "
        " title='Longest subpath of the best path in the optimal spanning tree "
        "that does not contain any vertices at maximum distance. "
        "This subpath is used to assemble sequence.'" <<
        ">Path used to assemble sequence</span>"

        "<br><input type=radio name=portionToDisplay value=none" <<
        (portionToDisplay=="none" ? " checked=on" : "") <<
        ">Nothing"



        "<br><input type=submit value='Display'>"

        " <span style='background-color:#e0e0e0' title='"
        "Fill this form to display a local subgraph of the global marker graph starting at the "
        "vertex containing the specified marker and extending out to a given distance "
        "(number of edges) from the start vertex."
        "The marker is specified by its oriented read id (read id, strand) and marker ordinal. "
        "The marker ordinal is the sequential index of the marker in the specified oriented read "
        "(the first marker is marker 0, the second marker is marker 1, and so on).'>"
        "Mouse here to explain form</span>"
        "</form>";
}



bool Assembler::LocalMarkerGraphRequestParameters::hasMissingRequiredParameters() const
{
    return
        !vertexIdIsPresent ||
        !maxDistanceIsPresent ||
        !minCoverageIsPresent ||
        !timeoutIsPresent;
}



void Assembler::showLocalMarkerGraphAlignments(
    ostream& html,
    const LocalMarkerGraph& graph,
    const LocalMarkerGraphRequestParameters& requestParameters
    )
{
    using vertex_descriptor = LocalMarkerGraph::vertex_descriptor;
    using edge_descriptor = LocalMarkerGraph::edge_descriptor;
    const size_t k = assemblerInfo->k;



    // Create a table with one column for each vertex
    // (in topologically sorted order), plus an additional column
    // between each pair of successive vertices.
    html << "<table style='table-layout:auto;white-space:nowrap;font-family:monospace;font-size:10px'>";


    // Row with vertex ids.
    if(requestParameters.showVertexId) {
        html << "<tr title='Vertex'><th style='text-align:left'>Vertex";
        for(size_t i=0; i<graph.topologicallySortedVertices.size(); i++) {
            const vertex_descriptor v = graph.topologicallySortedVertices[i];
            const LocalMarkerGraphVertex& vertex = graph[v];
            const auto vertexId = graph[v].vertexId;

            // Color.
            string color;
            if(vertex.distance == int(requestParameters.maxDistance)) {
                color = "cyan";
            } else if(vertex.distance == 0) {
                color = "lightGreen";
            } else if(vertex.markerInfos.size() >= requestParameters.minCoverage) {
                color = "green";
            } else {
                color = "red";
            }

            html << "<td style='background-color:" << color << ";text-align:center;font-size:8px'";
            if(requestParameters.portionToDisplay != "none") {
                html <<
                    " title='Click to position graph display at this vertex.'"
                    " onClick='positionAtVertex(" << vertexId << ")'";
            }
            html << ">";
            html << vertexId;

            // Add a column before the next vertex.
            if(i != graph.topologicallySortedVertices.size()-1) {
                html << "<td>";
            }
        }
    }



    // Row with vertex ranks.
    html << "<tr title='Rank'><th style='text-align:left' title='"
        "Vertex rank according to computed approximate topological sort."
        "'>Rank";
    for(size_t i=0; i<graph.topologicallySortedVertices.size(); i++) {
        const vertex_descriptor v = graph.topologicallySortedVertices[i];
        const LocalMarkerGraphVertex& vertex = graph[v];

        // Color.
        string color;
        if(vertex.distance == int(requestParameters.maxDistance)) {
            color = "cyan";
        } else if(vertex.distance == 0) {
            color = "lightGreen";
        } else if(vertex.markerInfos.size() >= requestParameters.minCoverage) {
            color = "green";
        } else {
            color = "red";
        }

        html << "<td style='background-color:" << color << ";text-align:center;cursor:pointer'";
        if(requestParameters.portionToDisplay != "none") {
            html <<
                " title='Vertex rank. Click to position graph display at this vertex.'"
                " onClick='positionAtVertex(" << vertex.vertexId << ")'";
        }
        html << ">";
        html << vertex.rank;

        // Add a column before the next vertex.
        if(i != graph.topologicallySortedVertices.size()-1) {
            html << "<td>";
        }
    }



    // Row with vertex distances.
    html << "<tr title='Distance'><th style='text-align:left' title='"
        "Vertex distance (number of edges) from start vertex."
        "'>Distance";
    for(size_t i=0; i<graph.topologicallySortedVertices.size(); i++) {
        const vertex_descriptor v = graph.topologicallySortedVertices[i];
        const LocalMarkerGraphVertex& vertex = graph[v];

        // Color.
        string color;
        if(vertex.distance == int(requestParameters.maxDistance)) {
            color = "cyan";
        } else if(vertex.distance == 0) {
            color = "lightGreen";
        } else if(vertex.markerInfos.size() >= requestParameters.minCoverage) {
            color = "green";
        } else {
            color = "red";
        }

        html << "<td style='background-color:" << color << ";text-align:center;cursor:pointer'";
        if(requestParameters.portionToDisplay != "none") {
            html <<
                " title='Vertex distance (number of edges) from start vertex. "
                "Click to position graph display at this vertex.'"
                " onClick='positionAtVertex(" << vertex.vertexId << ")'";
        }
        html << ">";
        html << vertex.distance;

        // Add a column before the next vertex.
        if(i != graph.topologicallySortedVertices.size()-1) {
            html << "<td>";
        }
    }



    // Row with marker k-mers.
    html << "<tr title='Marker k-mer'><th style='text-align:left'>Marker k-mer";
    for(size_t i=0; i<graph.topologicallySortedVertices.size(); i++) {
        const vertex_descriptor v = graph.topologicallySortedVertices[i];
        const LocalMarkerGraphVertex& vertex = graph[v];

        // Color.
        string color;
        if(vertex.distance == int(requestParameters.maxDistance)) {
            color = "cyan";
        } else if(vertex.distance == 0) {
            color = "lightGreen";
        } else if(vertex.markerInfos.size() >= requestParameters.minCoverage) {
            color = "green";
        } else {
            color = "red";
        }

        const KmerId kmerId = graph.getKmerId(v);
        const Kmer kmer(kmerId, k);

        html << "<td style='background-color:" << color << ";text-align:center;cursor:pointer'";
        if(requestParameters.portionToDisplay != "none") {
            html <<
                " title='Click to position graph display at this vertex.'"
                " onClick='positionAtVertex(" << vertex.vertexId << ")'";
        }
        html << ">";
        kmer.write(html, k);

        // Add a column before the next vertex.
        if(i != graph.topologicallySortedVertices.size()-1) {
            html << "<td>";
        }
    }



    // Row with assembled sequence.
    html <<
        "<tr title='Assembled sequence' style='background-color:pink'>"
        "<th style='text-align:left'>Assembled sequence";
    const vector<edge_descriptor>& assemblyPath = graph.localAssemblyPath;

    // Verify that rank increases along the assembly path.
    bool assemblyPathViolatesRankOrder = false;
    for(size_t i=0; i<assemblyPath.size(); i++) {
        const edge_descriptor e = assemblyPath[i];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        if(graph[v1].rank <= graph[v0].rank) {
            assemblyPathViolatesRankOrder = true;
            break;
        }
    }
    if(assemblyPathViolatesRankOrder) {
        html <<
            "<td colspan=" << 2*graph.topologicallySortedVertices.size()-1 << ">"
            "Assembly path violates rank order.";
    }

    // If necessary, add an empty cell at the beginning, covering up to and excluding
    // the vertex with lowest rank.
    const edge_descriptor firstEdge = assemblyPath.front();
    const vertex_descriptor firstVertex = source(firstEdge, graph);
    const size_t lowestRank = graph[firstVertex].rank;
    if(lowestRank > 0) {
        html << "<td colspan=" << 2*lowestRank << ">";
    }

    // Write the marker of the first vertex.
    const KmerId firstKmerId = graph.getKmerId(firstVertex);
    const Kmer firstKmer(firstKmerId, k);
    html << "<td>";
    for(size_t i=0; i<k; i++) {
        const Base base = firstKmer[i];
        char c = base.character();
        if(graph[firstVertex].markerInfos.size() < 10) {
            c = char(std::tolower(c));
        }
        html << c;
    }

    // To add assembled sequence to the alignment table,
    // loop over edges of the assembly path.
    // The code here is similar to LocalMarkerGraph::assembleDominantSequence.
    for(const edge_descriptor e: assemblyPath) {
        const LocalMarkerGraphEdge& edge = graph[e];
        const auto consensus = edge.consensus();
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const LocalMarkerGraphVertex& vertex0 = graph[v0];
        const LocalMarkerGraphVertex& vertex1 = graph[v1];
        const size_t rank0 = vertex0.rank;
        const size_t rank1 = vertex1.rank;
        CZI_ASSERT(rank0 < rank1);  // We checked for this above.

        // Write the edge sequence, if any.
        CZI_ASSERT(!edge.infos.empty());
        const auto& p = edge.infos.front();
        // const auto coverage = p.second.size();
        const auto& edgeSequence = p.first.sequence;
        html << "<td colspan=" << 2*(rank1-rank0) - 1 << ">";
        for(const shasta::Base b: edgeSequence) {
            char c = b.character();
            if(consensus < 10) {
                c = char(std::tolower(c));
            }
            html << c;
        }

        // Write the sequence of the target vertex.
        const KmerId kmerId1 = graph.getKmerId(v1);
        const Kmer kmer1(kmerId1, k);
        const auto coverage1 = graph[v1].markerInfos.size();
        html << "<td style='text-align:right'>";
        for(size_t i=size_t(p.first.overlappingBaseCount); i<k; i++) {
            const Base b = kmer1[i];
            char c = b.character();
            if(coverage1 < 10) {
                c = char(std::tolower(c));
            }
            html << c;
        }
    }

    // If necessary, also add an empty cell at the end.
    const edge_descriptor lastEdge = assemblyPath.back();
    const vertex_descriptor lastVertex = target(lastEdge, graph);
    const size_t highestRank = graph[lastVertex].rank;
    if(highestRank < graph.topologicallySortedVertices.size()-1) {
        html << "<td colspan=" << 2*(graph.topologicallySortedVertices.size()-1 -highestRank) << ">";
    }



    // Repeat the  loop over edges of the assembly path
    // to add a row with assembled position.
    html << "<tr style='background-color:pink' title='Assembled position'>"
        "<th style='text-align:left'>Assembled position";
    if(lowestRank > 0) {
        html << "<td colspan=" << 2*lowestRank << ">";
    }
    size_t position = 0;
    html << "<td class=centered>" << position;
    position += k;
    for(const edge_descriptor e: assemblyPath) {
        const LocalMarkerGraphEdge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const LocalMarkerGraphVertex& vertex0 = graph[v0];
        const LocalMarkerGraphVertex& vertex1 = graph[v1];
        const size_t rank0 = vertex0.rank;
        const size_t rank1 = vertex1.rank;
        CZI_ASSERT(rank0 < rank1);  // We checked for this above.
        const auto& p = edge.infos.front();
        const auto& edgeSequence = p.first.sequence;

        html <<
            "<td class=centered colspan=" << 2*(rank1-rank0) - 1 <<
            ">" << position;
        position += edgeSequence.size();
        html << "<td class=centered>" << position + p.first.overlappingBaseCount;
        position += k - p.first.overlappingBaseCount;
    }
    if(highestRank < graph.topologicallySortedVertices.size()-1) {
        html << "<td colspan=" << 2*(graph.topologicallySortedVertices.size()-1 -highestRank) << ">";
    }



    // Repeat the  loop over edges of the assembly path
    // to add a row for vertex and edge coverage.
    html << "<tr style='background-color:pink' title='Coverage'><th style='text-align:left'>Coverage";
    if(lowestRank > 0) {
        html << "<td colspan=" << 2*lowestRank << ">";
    }
    html << "<td class=centered>" << graph[firstVertex].markerInfos.size();
    for(const edge_descriptor e: assemblyPath) {
        const LocalMarkerGraphEdge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const LocalMarkerGraphVertex& vertex0 = graph[v0];
        const LocalMarkerGraphVertex& vertex1 = graph[v1];
        const size_t rank0 = vertex0.rank;
        const size_t rank1 = vertex1.rank;
        CZI_ASSERT(rank0 < rank1);  // We checked for this above.

        html <<
            "<td class=centered colspan=" << 2*(rank1-rank0) - 1 <<
            ">" << edge.coverage() << "<td class=centered>" << vertex1.markerInfos.size();
    }
    if(highestRank < graph.topologicallySortedVertices.size()-1) {
        html << "<td colspan=" << 2*(graph.topologicallySortedVertices.size()-1 -highestRank) << ">";
    }



    // Repeat the  loop over edges of the assembly path
    // to add a row for edge consensus.
    html << "<tr style='background-color:pink' title='Consensus'><th style='text-align:left'>Consensus";
    if(lowestRank > 0) {
        html << "<td colspan=" << 2*lowestRank << ">";
    }
    html << "<td>"; // First vertex
    for(const edge_descriptor e: assemblyPath) {
        const LocalMarkerGraphEdge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const LocalMarkerGraphVertex& vertex0 = graph[v0];
        const LocalMarkerGraphVertex& vertex1 = graph[v1];
        const size_t rank0 = vertex0.rank;
        const size_t rank1 = vertex1.rank;
        CZI_ASSERT(rank0 < rank1);  // We checked for this above.

        html <<
            "<td class=centered colspan=" << 2*(rank1-rank0) - 1 <<
            ">" << edge.consensus() << "<td>";
    }
    if(highestRank < graph.topologicallySortedVertices.size()-1) {
        html << "<td colspan=" << 2*(graph.topologicallySortedVertices.size()-1 -highestRank) << ">";
    }



    // Add a row for each oriented read.
    for(const OrientedReadId orientedReadId: graph.orientedReadIds) {
        html <<
            "<tr title='Oriented read " << orientedReadId <<
            "'><td><a href='exploreRead?readId&amp;" << orientedReadId.getReadId() <<
            "&amp;strand=" << orientedReadId.getStrand() << "'>" <<
            orientedReadId << "</a>";

        // Find the graph vertices corresponding to this oriented read,
        // and the corresponding ordinals.
        vector< pair<uint32_t, vertex_descriptor> > orientedReadVertices;
        graph.getOrientedReadVertices(orientedReadId, orientedReadVertices);

#if 0
        // Write out details of the vertices.
        cout << orientedReadId << " vertices:" << endl;
        for(size_t i=0; i<orientedReadVertices.size(); i++) {
            const auto& p = orientedReadVertices[i];
            const uint32_t ordinal = p.first;
            const vertex_descriptor v = p.second;
            const LocalMarkerGraphVertex& vertex = graph[v];
            cout << "Vertex " << vertex.vertexId;
            cout << " distance " << vertex.distance;
            cout << " rank " << vertex.rank;
            cout << " ordinal " << ordinal;
            if(i>0 && graph[orientedReadVertices[i-1].second].rank >= vertex.rank) {
                cout << " rank order violation";
            }
            cout << endl;
        }
#endif

        // Get rid of the vertices that are inconsistent with rank ordering.
        vector< pair<uint32_t, vertex_descriptor> > rankEnforcedOrientedReadVertices;
        graph.enforceRankOrder(orientedReadVertices, rankEnforcedOrientedReadVertices);

#if 0
        // Write out details of the vertices.
        for(size_t i=0; i<rankEnforcedOrientedReadVertices.size(); i++) {
            const auto& p = rankEnforcedOrientedReadVertices[i];
            const uint32_t ordinal = p.first;
            const vertex_descriptor v = p.second;
            const LocalMarkerGraphVertex& vertex = graph[v];
            cout << "Vertex " << vertex.vertexId;
            cout << " distance " << vertex.distance;
            cout << " rank " << vertex.rank;
            cout << " ordinal " << ordinal;
            cout << endl;
            if(i>0) {
                CZI_ASSERT(graph[rankEnforcedOrientedReadVertices[i-1].second].rank < vertex.rank);
            }
        }
#endif

        // If we have no vertices, just write an empty cell covering the entire width.
        if(rankEnforcedOrientedReadVertices.empty()) {
            html << "<td colspan=" << 2*graph.topologicallySortedVertices.size()-1 << ">";
            continue;
        }

        // If necessary, add an empty cell at the beginning, covering up to and excluding
        // the vertex with lowest rank.
        const size_t lowestRank = graph[rankEnforcedOrientedReadVertices.front().second].rank;
        if(lowestRank > 0) {
            html << "<td colspan=" << 2*lowestRank << ">";
        }



        // Loop over the vertices.
        size_t skipBaseCount = 0;
        for(size_t i=0; i<rankEnforcedOrientedReadVertices.size(); i++) {

            // Get the information we need about this vertex.
            const auto& p0 = rankEnforcedOrientedReadVertices[i];
            const uint32_t ordinal0 = p0.first;
            const vertex_descriptor v0 = p0.second;
            const LocalMarkerGraphVertex& vertex0 = graph[v0];
            const size_t rank0 = vertex0.rank;
            const KmerId kmerId0 = graph.getKmerId(v0);
            const Kmer kmer0(kmerId0, k);
            const MarkerId markerId0 = getMarkerId(orientedReadId, ordinal0);
            const CompressedMarker& marker0 = markers.begin()[markerId0];

            // Write the k-mer.
            html << "<td style='text-align:right' title='Oriented read " << orientedReadId <<
                " marker ordinal " << ordinal0 <<
                ", positions " << marker0.position << "-" << marker0.position+k-1;
            if(skipBaseCount>0) {
                html << " (first " << skipBaseCount << " bases overlap with previous marker and are not displayed)";
            }
            html <<
                "'>"
                "<a style='text-decoration:none;color:black' href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&amp;strand=" << orientedReadId.getStrand() <<
                "&amp;highlightMarker=" << ordinal0 <<
                "'>";
            for(size_t j=skipBaseCount; j<k; j++) {
                html << kmer0[j];
            }
            html << "</a>";

            // If this is the last vertex, we are done.
            if(i == rankEnforcedOrientedReadVertices.size()-1) {
                continue;
            }

            // Get the information we need about the next vertex.
            const auto& p1 = rankEnforcedOrientedReadVertices[i+1];
            const uint32_t ordinal1 = p1.first;
            const vertex_descriptor v1 = p1.second;
            const LocalMarkerGraphVertex& vertex1 = graph[v1];
            const size_t rank1 = vertex1.rank;

            // Write the cell in between these two vertices.
            const MarkerId markerId1 = getMarkerId(orientedReadId, ordinal1);
            const CompressedMarker& marker1 = markers.begin()[markerId1];
            const uint32_t position0 = marker0.position;
            const uint32_t position1 = marker1.position;
            html << "<td colspan=" << 2*(rank1-rank0)-1;
            if(position1-1 >= position0+k) {
                html <<
                    " title='Oriented read " << orientedReadId << " between markers " <<
                    ordinal0 << "-" << ordinal1 <<
                    ", positions " << position0+k << "-" << position1-1 << "'";
            }
            html <<
                ">"
                "<a style='text-decoration:none;color:black' href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&amp;strand=" << orientedReadId.getStrand() <<
                "&amp;highlightMarker=" << ordinal0 <<
                "&amp;highlightMarker=" << ordinal1 <<
                "'>";

            // Write the sequence in between the two vertices.
            if(position1 > position0+k) {
                if(position1-position0-k > 100) {
                    html << "Too long";
                } else {
                    for(uint32_t position=position0+uint32_t(k); position<position1; position++) {
                        html << getOrientedReadBase(orientedReadId, position);
                    }
                }
                skipBaseCount = 0;
            } else {
                skipBaseCount = position0+k-position1;
            }
            html << "</a>";
        }


        // If necessary, also add an empty cell at the end
        const size_t highestRank = graph[rankEnforcedOrientedReadVertices.back().second].rank;
        if(highestRank < graph.topologicallySortedVertices.size()-1) {
            html << "<td colspan=" << 2*(graph.topologicallySortedVertices.size()-1 - lowestRank) << ">";
        }
    }



    // Finish the table
    html << "</table>";
}



void Assembler::exploreMarkerGraphVertex(const vector<string>& request, ostream& html)
{
    // Get the vertex id.
    GlobalMarkerGraphVertexId vertexId = 0;
    const bool vertexIdIsPresent = getParameterValue(request, "vertexId", vertexId);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Show details for marker graph vertex'> "
        "<input type=text name=vertexId required" <<
        (vertexIdIsPresent ? (" value=" + to_string(vertexId)) : "") <<
        " size=8 title='Enter a vertex id between 0 and " << globalMarkerGraphVertices.size()-1 << "'>";
    html << "</form>";

    // If the vertex id missing or invalid, stop here.
    if(!vertexIdIsPresent || !vertexIdIsPresent) {
        return;
    }
    if(vertexId >= globalMarkerGraphVertices.size()) {
        html << "<p>Invalid vertex id. Must be less than " << globalMarkerGraphVertices.size() << ".";
        return;
    }

    // Access the markers of this vertex.
    MemoryAsContainer<MarkerId> markerIds = globalMarkerGraphVertices[vertexId];
    const size_t markerCount = markerIds.size();
    CZI_ASSERT(markerCount > 0);

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
            CZI_ASSERT(base == kmer[i]);
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
    for(size_t i=0; i<k; i++) {

        Coverage coverage;
        for(size_t j=0; j<markerCount; j++) {
            coverage.addRead(
                AlignedBase(kmer[i]),
                orientedReadIds[j].getStrand(),
                repeatCounts[j][i]);
        }

        const Consensus consensus = (*consensusCaller)(coverage);
        CZI_ASSERT(Base(consensus.base) == kmer[i]);
        consensusRepeatCounts[i] = consensus.repeatCount;
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
    html << "<h1>Marker graph vertex "<< vertexId << "</h1>";
    html << "<p>Vertex coverage (number of markers) is " << markerCount << ".";



    // Write a table with one row for each marker.
    html <<
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
    GlobalMarkerGraphEdgeId edgeId = 0;
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
    array<GlobalMarkerGraphVertexId, 2> vertexIds = {edge.source, edge.target};
    const size_t markerCount = edge.coverage;

    // The marker intervals of this edge.
    const MemoryAsContainer<MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
    CZI_ASSERT(markerIntervals.size() == markerCount);

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



    // Initialize a spoa multiple sequence alignment.
    const spoa::AlignmentType alignmentType = spoa::AlignmentType::kSW;
    const int8_t match = 1;
    const int8_t mismatch = -1;
    const int8_t gap = -1;
    auto alignmentEngine = spoa::createAlignmentEngine(alignmentType, match, mismatch, gap);
    auto alignmentGraph = spoa::createGraph();


    // Add all the sequences to the alignment,
    // including the flanking markers.
    string sequenceString;
    vector<string> msa;
    for(size_t j=0; j!=markerCount; j++) {

        // Get the sequence.
        sequenceString.clear();
        for(const Base base: sequences[j]) {
            sequenceString.push_back(base.character());
        }

        // Add it to the alignment.
        auto alignment = alignmentEngine->align_sequence_with_graph(sequenceString, alignmentGraph);
        alignmentGraph->add_alignment(alignment, sequenceString);
    }

    // Use spoa to compute the multiple sequence alignment.
    msa.clear();
    alignmentGraph->generate_multiple_sequence_alignment(msa);

    // The length of the alignment, including gaps.
    const size_t alignmentLength = msa.front().size();



    // Compute coverage for each base (including "-")
    // at each position of the alignment.
    vector< array<size_t, 5> > baseCoverage(alignmentLength);
    for(size_t i=0; i<alignmentLength; i++) {
        fill(baseCoverage[i].begin(), baseCoverage[i].end(), 0);
        for(size_t j=0; j<markerCount; j++) {
            const AlignedBase base = AlignedBase::fromCharacter(msa[j][i]);
            ++(baseCoverage[i][base.value]);
        }
    }



    // Use the consensus caller to compute the consensus base
    // and repeat count at each position of the alignment.
    vector<AlignedBase> consensusBase(alignmentLength);
    vector<size_t> consensusRepeatCount(alignmentLength);
    vector<size_t> positions(markerCount, 0);
    for(size_t i=0; i<alignmentLength; i++) {
        Coverage coverage;
        for(size_t j=0; j<markerCount; j++) {
            const Strand strand = markerIntervals[j].orientedReadId.getStrand();
            const AlignedBase base = AlignedBase::fromCharacter(msa[j][i]);
            if(base.isGap()) {
                coverage.addRead(base, strand, 0);  // Gap always gets repeat count of 0.
            } else {
                const AlignedBase base = AlignedBase(sequences[j][positions[j]]);
                const size_t repeatCount = repeatCounts[j][positions[j]];
                ++(positions[j]);
                coverage.addRead(base, strand, repeatCount);
            }

        }
        const Consensus consensus = (*consensusCaller)(coverage);
        consensusBase[i] = consensus.base;
        consensusRepeatCount[i] = consensus.repeatCount;
    }



    // Count the number of gaps in the consensus bases.
    const size_t consensusGapCount =
        count(consensusBase.begin(), consensusBase.end(), AlignedBase::gap());

    // Compute concordant and discordant base coverage at each position.
    vector<size_t> concordantBaseCoverage(alignmentLength, 0);
    vector<size_t> discordantBaseCoverage(alignmentLength, 0);
    for(size_t i=0; i<alignmentLength; i++) {
        for(size_t b=0; b<5; b++) {
            const size_t coverage = baseCoverage[i][b];
            if(consensusBase[i] == AlignedBase::fromInteger(b)) {
                concordantBaseCoverage[i] += coverage;
            } else {
                discordantBaseCoverage[i] += coverage;
            }
        }
    }



    // Page title.
    html << "<h1>Marker graph edge "<< edgeId << "</h1>";

    // Basic information about this edge.
    html << "<p>Source vertex " << vertexIds[0];
    html << ", target vertex " << vertexIds[1] << ".";
    html << "<p>Edge coverage (number of markers) is " << markerCount << ".";

    // Begin the table.
    html <<
        "<table>"
        "<tr>"
        "<th>Oriented<br>read"
        "<th>Ordinal0"
        "<th>Ordinal1"
        "<th>Sequence<br>(includes flanking markers)"
        "<th>Repeat<br>counts";




    // Write one row for each oriented read.
    for(size_t j=0; j<markerCount; j++) {
        const MarkerInterval& markerInterval = markerIntervals[j];
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const ReadId readId = orientedReadId.getReadId();
        const Strand strand = orientedReadId.getStrand();

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

        // Sequence, written aligned.
        html << "<td class=centered style='font-family:monospace'>";
        size_t position = 0;
        for(size_t i=0; i<alignmentLength; i++) {
            const char alignmentCharacter = msa[j][i];
            if(alignmentCharacter == '-') {
                html << "-";
            } else {
                if(sequences[j].size() > 2*k && position==k) {
                    html << "<span style='background-color:LightGreen'>";
                }
                const Base base = sequences[j][position];
                ++position;
                CZI_ASSERT(base == Base::fromCharacter(alignmentCharacter));
                html << base;
                if(sequences[j].size() > 2*k && position==sequences[j].size()-k) {
                    html << "</span>";
                }
            }
        }
        CZI_ASSERT(position == sequences[j].size());
        CZI_ASSERT(position == repeatCounts[j].size());



        // Repeat counts, also written aligned.
        html << "<td class=centered style='font-family:monospace'>";
        position = 0;
        for(size_t i=0; i<alignmentLength; i++) {
            const char alignmentCharacter = msa[j][i];
            if(alignmentCharacter == '-') {
                html << "-";
            } else {
                if(sequences[j].size() > 2*k && position==k) {
                    html << "<span style='background-color:LightGreen'>";
                }
                const Base base = sequences[j][position];
                const uint8_t repeatCount = repeatCounts[j][position];
                ++position;
                CZI_ASSERT(base == Base::fromCharacter(alignmentCharacter));
                if(repeatCount < 10) {
                    html << int(repeatCount);
                } else {
                    html << "*";
                }
                if(sequences[j].size() > 2*k && position==sequences[j].size()-k) {
                    html << "</span>";
                }
            }
        }
        CZI_ASSERT(position == sequences[j].size());
        CZI_ASSERT(position == repeatCounts[j].size());
    }



    // Consensus bases and repeat counts.
    html <<
        "<tr>"
        "<th colspan=3 class=left>Consensus (run-length)"
        "<td class=centered style='font-family:monospace'>";
    size_t position = 0;
    for(size_t i=0; i<alignmentLength; i++) {
        if(consensusBase[i].isGap()) {
            html << "-";
        } else {
            if(alignmentLength-consensusGapCount > 2*k && position==k) {
                html << "<span style='background-color:LightGreen'>";
            }
            html << consensusBase[i];
            ++position;
            if(alignmentLength-consensusGapCount > 2*k && position==alignmentLength-consensusGapCount-k) {
                html << "</span>";
            }
        }
    }
    CZI_ASSERT(position == alignmentLength-consensusGapCount);
    html << "<td class=centered style='font-family:monospace'>";
    position = 0;
    for(size_t i=0; i<alignmentLength; i++) {
        if(consensusBase[i].isGap()) {
            html << "-";
        } else {
            if(alignmentLength-consensusGapCount > 2*k && position==k) {
                html << "<span style='background-color:LightGreen'>";
            }
            const size_t repeatCount = consensusRepeatCount[i];
            if(repeatCount < 10) {
                html << repeatCount;
            } else {
                html << "*";
            }
            ++position;
            if(alignmentLength-consensusGapCount > 2*k && position==alignmentLength-consensusGapCount-k) {
                html << "</span>";
            }
        }
    }
    CZI_ASSERT(position == alignmentLength-consensusGapCount);



    // Concordant base coverage.
    html <<
        "<tr>"
        "<th colspan=3 class=left>Concordant base coverage"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<alignmentLength; i++) {
        const size_t coverage = concordantBaseCoverage[i];
        if(coverage == 0) {
            html << ".";
        } else if(coverage < 10) {
            html << coverage;
        } else {
            html << "*";
        }
    }
    html << "<td>";



    // Discordant base coverage.
    html <<
        "<tr>"
        "<th colspan=3 class=left>Discordant base coverage"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<alignmentLength; i++) {
        const size_t coverage = discordantBaseCoverage[i];
        if(coverage == 0) {
            html << ".";
        } else if(coverage < 10) {
            html << coverage;
        } else {
            html << "*";
        }
    }
    html << "<td>";



    // Base coverage.
    for(size_t b=0; b<5; b++) {
        html <<
            "<tr>"
            "<th colspan=3 class=left>Coverage for " << AlignedBase::fromInteger(b) <<
            "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<alignmentLength; i++) {
            const size_t coverage = baseCoverage[i][b];
            if(coverage == 0) {
                html << ".";
            } else if(coverage < 10) {
                html << coverage;
            } else {
                html << "*";
            }
        }
        html << "<td>";
    }



    // Raw consensus sequence.
    html <<
        "<tr>"
        "<th colspan=3 class=left>Consensus (raw)"
        "<td class=centered style='font-family:monospace'>";
    position = 0;
    for(size_t i=0; i<alignmentLength; i++) {
        const AlignedBase base = consensusBase[i];

        if(base.isGap()) {
            continue;
        }

        if(alignmentLength-consensusGapCount > 2*k && position==k) {
            html << "<span style='background-color:LightGreen'>";
        }

        const size_t repeatCount = consensusRepeatCount[i];
        for(size_t k=0; k<repeatCount; k++) {
            html << base;
        }
        ++position;

        if(alignmentLength-consensusGapCount > 2*k && position==alignmentLength-consensusGapCount-k) {
            html << "</span>";
        }
    }
    html << "<td>";



    // End the table.
    html << "</table>";
}
