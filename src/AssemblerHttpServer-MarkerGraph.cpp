// Shasta.
#include "Assembler.hpp"
#include "LocalMarkerGraph2.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

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
    html <<
        "<h3>Display a local subgraph of the global marker graph</h3>"
        "<form>"

        "<table>"

        "<tr title='Read id between 0 and " << reads.size()-1 << "'>"
        "<td>Read id"
        "<td><input type=text required name=readId size=8 style='text-align:center'"
        << (requestParameters.readIdIsPresent ? ("value='"+to_string(requestParameters.readId)+"'") : "") <<
        ">"

        "<tr title='Choose 0 (+) for the input read or 1 (-) for its reverse complement'>"
        "<td>Strand"
        "<td class=centered>";
    writeStrandSelection(html, "strand",
        requestParameters.strandIsPresent && requestParameters.strand==0,
        requestParameters.strandIsPresent && requestParameters.strand==1);

    html <<
        "<tr title='Ordinal for the desired marker in the specified oriented read.'>"
        "<td>Marker ordinal"
        "<td><input type=text required name=ordinal size=8 style='text-align:center'"
        << (requestParameters.ordinalIsPresent ? ("value='"+to_string(requestParameters.ordinal)+"'") : "") <<
        ">"

        "<tr title='Maximum distance from start vertex (number of edges)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        << (requestParameters.maxDistanceIsPresent ? ("value='" + to_string(requestParameters.maxDistance)+"'") : " value='6'") <<
        ">"

        "<tr title='Check for detailed graph with labels'>"
        "<td>Detailed"
        "<td class=centered><input type=checkbox name=detailed"
        << (requestParameters.detailed ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to show vertex ids (only useful for debugging)'>"
        "<td>Show vertex ids"
        "<td class=centered><input type=checkbox name=showVertexId"
        << (requestParameters.showVertexId ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to show in purple the optimal spanning tree of the local marker graph"
        " (always shown if assembled sequence is also selected)'>"
        "<td>Show optimal spanning tree"
        "<td class=centered><input type=checkbox name=showOptimalSpanningTree"
        << (requestParameters.showOptimalSpanningTree ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to show sequence assembled from this local marker graph'>"
        "<td>Show assembled sequence"
        "<td class=centered><input type=checkbox name=showAssembledSequence"
        << (requestParameters.showAssembledSequence ? " checked=checked" : "") <<
        ">"


        "<tr title='Minimum coverage (number of markers) for a vertex or edge to be considered strong. "
        "Affects the coloring of vertices and edges.'>"
        "<td>Coverage threshold"
        "<td><input type=text required name=minCoverage size=8 style='text-align:center'"
        << (requestParameters.minCoverageIsPresent ? (" value='" + to_string(requestParameters.minCoverage)+"'") : " value='3'") <<
        ">"

        "<tr title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "<td>Graphics size in pixels"
        "<td><input type=text required name=sizePixels size=8 style='text-align:center'"
        << (requestParameters.sizePixelsIsPresent ? (" value='" + to_string(requestParameters.sizePixels)+"'") : " value='1600'") <<
        ">"

        "<tr>"
        "<td>Timeout (seconds) for graph creation and layout"
        "<td><input type=text required name=timeout size=8 style='text-align:center'"
        << (requestParameters.timeoutIsPresent ? (" value='" + to_string(requestParameters.timeout)+"'") : " value='30'") <<
        ">"
        "</table>"



        // Radio buttons to choose what portion of the graph to display.
        "Portion of the graph to display:"

        "<br><input type=radio name=portionToDisplay value=all" <<
        (requestParameters.portionToDisplay=="all" ? " checked=on" : "") <<
        ">Entire graph"

        "<br><input type=radio name=portionToDisplay value=spanningTree" <<
        (requestParameters.portionToDisplay=="spanningTree" ? " checked=on" : "") <<
        ">Optimal spanning tree only"

        "<br><input type=radio name=portionToDisplay value=optimalPath" <<
        (requestParameters.portionToDisplay=="optimalPath" ? " checked=on" : "") <<
        ">Optimal spanning tree best path only"

        "<br><input type=radio name=portionToDisplay value=clippedOptimalPath" <<
        (requestParameters.portionToDisplay=="clippedOptimalPath" ? " checked=on" : "") <<
        "><span "
        " title='Longest subpath of the best path in the optimal spanning tree "
        "that does not contain any vertices at maximum distance. "
        "This subpath is used to assemble sequence.'" <<
        ">Path used to assemble sequence</span>"

        "<br><input type=radio name=portionToDisplay value=none" <<
        (requestParameters.portionToDisplay=="none" ? " checked=on" : "") <<
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



    // If any required values are missing, stop here.
    if(!requestParameters.readIdIsPresent || !requestParameters.strandIsPresent || !requestParameters.ordinalIsPresent
        || !requestParameters.maxDistanceIsPresent || !requestParameters.minCoverageIsPresent
        || !requestParameters.timeoutIsPresent) {
        return;
    }



    // Validity checks.
    if(requestParameters.readId > reads.size()) {
        html << "<p>Invalid read id " << requestParameters.readId;
        html << ". Must be between 0 and " << reads.size()-1 << ".";
        return;
    }
    if(requestParameters.strand>1) {
        html << "<p>Invalid strand " << requestParameters.strand;
        html << ". Must be 0 or 1.";
        return;
    }
    const OrientedReadId orientedReadId(requestParameters.readId, requestParameters.strand);
    const auto orientedReadMarkerCount = markers.size(orientedReadId.getValue());
    if(requestParameters.ordinal >= orientedReadMarkerCount) {
        html <<
            "<p>Invalid marker ordinal. "
            "Oriented read " << orientedReadId <<
            " has "  << orientedReadMarkerCount <<
            " markers, and therefore the ordinal must be"
            " between 0 and " << orientedReadMarkerCount-1 << ".";
        return;
    }



    // Create the local marker graph.
    LocalMarkerGraph2 graph(uint32_t(assemblerInfo->k), reads, markers, globalMarkerGraphVertex);
    const auto createStartTime = steady_clock::now();
    if(!extractLocalMarkerGraph(
        orientedReadId,
        requestParameters.ordinal,
        requestParameters.maxDistance,
        requestParameters.timeout,
        graph)) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    if(num_vertices(graph) == 0) {
        html << "<p>The specified marker does not correspond to a vertex of the marker graph.";
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



    // Write the sequence, if requested.
    if(requestParameters.showAssembledSequence) {

        const string fastaSequenceName =
            "AssembledSequence-" + orientedReadId.getString() +
            "-" + to_string(requestParameters.ordinal) + "-" + to_string(requestParameters.maxDistance);
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

        html <<
            "<h2>Sequence assembled from this local marker graph</h2>"
            "<br><span title='Coverage on second line, blank if &ge; 10."
            " Base is lower case if coverage &lt; 10.'>"
            "Assembled " << sequence.size() << " bases. </span>"
            "<a id=fastaDownload>Download in FASTA format</a>"
            "<script>"
            "var element = document.getElementById('fastaDownload');"
            "element.setAttribute('href', 'data:text/plain;charset=utf-8,' +"
            "encodeURIComponent('" << fastaString << "'));"
            "element.setAttribute('download', '" << fastaFileName << "');"
            "</script>"
            "<pre>";

        // Labels.
        for(size_t position=0; position<sequence.size(); position+=10) {
            const string label = to_string(position);
            html << label;
            for(size_t i=0; i<10-label.size(); i++) {
                html << " ";
            }
        }
        html << "\n";

        // Scale.
        for(size_t position=0; position<sequence.size(); position++) {
            if((position%10)==0) {
                html << "|";
            } else if((position%5)==0) {
                html << "+";
            } else {
                html << ".";
            }
        }
        html << "\n";

        // Sequence.
        for(const auto& p: sequence) {
            const auto coverage = p.second;
            if(coverage < 10) {
                html << char(std::tolower(p.first.character()));
            } else {
                html << p.first;
            }
        }
        html << "\n";

        // Coverage.
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
            "<h2>Marker graph near marker " << requestParameters.ordinal <<
            " of oriented read " << orientedReadId <<
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
        BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph2) {
            const LocalMarkerGraph2Vertex& vertex = graph[v];
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
                (requestParameters.showOptimalSpanningTree ? "&showOptimalSpanningTree=on" : "") +
                (requestParameters.showAssembledSequence ? "&showAssembledSequence=on" : "");
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
                    (requestParameters.showOptimalSpanningTree ? "&showOptimalSpanningTree=on" : "") +
                    (requestParameters.showAssembledSequence ? "&showAssembledSequence=on" : "");
                html <<
                    "document.getElementById('vertex" << vertex.vertexId <<
                    "').oncontextmenu = function() {location.href='" << detailUrl << "';"
                    "return false;};\n";
            }
        }
        html << "</script>\n";



        // Position the start vertex at the center of the window.
        const GlobalMarkerGraphVertexId startVertexId =
            getGlobalMarkerGraphVertex(orientedReadId, requestParameters.ordinal);
        html << "<script>\n";
        if(requestParameters.detailed) {
            html <<
                "var element = document.getElementById('a_vertexDistance" << startVertexId << "');\n";
        } else {
            html <<
                "var element = document.getElementById('vertex" << startVertexId << "');\n";
        }
        html <<
            "var r = element.getBoundingClientRect();\n"
            "window.scrollBy((r.left + r.right - window.innerWidth) / 2, (r.top + r.bottom - window.innerHeight) / 2);\n"
            "</script>\n";
    }
}



// Extract  from the request the parameters for the display
// of the local marker graph.
void Assembler::getLocalMarkerGraphRequestParameters(
    const vector<string>& request,
    LocalMarkerGraphRequestParameters& parameters) const
{
    parameters.readId = 0;
    parameters.readIdIsPresent = getParameterValue(
        request, "readId", parameters.readId);

    parameters.strand = 0;
    parameters.strandIsPresent = getParameterValue(
        request, "strand", parameters.strand);

    parameters.ordinal = 0;
    parameters.ordinalIsPresent = getParameterValue(
        request, "ordinal", parameters.ordinal);

    parameters.maxDistance = 0;
    parameters.maxDistanceIsPresent = getParameterValue(
        request, "maxDistance", parameters.maxDistance);

    string detailedString;
    parameters.detailed = getParameterValue(
        request, "detailed", detailedString);

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
