#ifdef SHASTA_HTTP_SERVER

// Shasta.
#include "Assembler.hpp"
#include "AssembledSegment.hpp"
#include "LocalAssemblyGraph.hpp"
#include "platformDependent.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "chrono.hpp"



void Assembler::exploreAssemblyGraph(
    const vector<string>& request,
    ostream& html)
{

    AssemblyGraph& assemblyGraph = *assemblyGraphPointer;

    // Get the request parameters.
    LocalAssemblyGraphRequestParameters requestParameters;
    getLocalAssemblyGraphRequestParameters(request, requestParameters);

    // Write the form.
    const bool allowHighlighting = phasingData.assemblyGraphEdges.isOpen();
    requestParameters.writeForm(html, assemblyGraph.edges.size(), allowHighlighting);

    // If any required values are missing, stop here.
    if(requestParameters.hasMissingRequiredParameters()) {
        return;
    }

    // Validity check.
    if(requestParameters.edgeId > assemblyGraph.edges.size()) {
        html << "<p>Invalid edge id " << requestParameters.edgeId;
        html << ". Must be between 0 and " << assemblyGraph.edges.size()-1 << " inclusive.";
        return;
    }



    // Create the local assembly graph.
    html << "<h1>Local assembly graph</h1>";
    LocalAssemblyGraph graph(assemblyGraph);
    const auto createStartTime = steady_clock::now();
    if(!extractLocalAssemblyGraph(
        requestParameters.edgeId,
        requestParameters.maxDistance,
        requestParameters.timeout,
        graph)) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    html << "<p>The local assembly graph has " << num_vertices(graph);
    html << " vertices and " << num_edges(graph) << " edges.";

    const auto createFinishTime = steady_clock::now();
    if(seconds(createFinishTime - createStartTime) > requestParameters.timeout) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }


#if 0
    // Highlight edges containing the specified oriented read.
    if(allowHighlighting && requestParameters.highlightedReadIdIsPresent) {

        // Create a map of local edges, keyed by global edge id.
        std::map<AssemblyGraph::EdgeId, LocalAssemblyGraph::edge_descriptor> edgeMap;
        BGL_FORALL_EDGES(e, graph, LocalAssemblyGraph) {
            edgeMap.insert(make_pair(graph[e].edgeId, e));
        }

        // Get the edges we need to highlight.
        const OrientedReadId orientedReadId(
            requestParameters.highlightedReadId,
            requestParameters.highlightedStrand);
        const span<AssemblyGraph::EdgeId> highlightedEdges =
            phasingGraph.assemblyGraphEdges[orientedReadId.getValue()];

        // Highlight those edges, if present in the local assembly graph.
        for(const AssemblyGraph::EdgeId edgeId: highlightedEdges) {
            const auto it = edgeMap.find(edgeId);
            if(it == edgeMap.end()) {
                continue;
            }
            const LocalAssemblyGraph::edge_descriptor e = it->second;
            graph[e].isHighlighted = true;
        }

    }
#endif


    // Write it out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    graph.write(
        dotFileName,
        requestParameters.maxDistance,
        requestParameters.useDotLayout,
        requestParameters.showVertexLabels,
        requestParameters.showEdgeLabels);



    // Compute graph layout in svg format.
    const string command =
        timeoutCommand() + " " + to_string(requestParameters.timeout - seconds(createFinishTime - createStartTime)) +
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
    // filesystem::remove(dotFileName);

    // Copy the svg to html.
    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << "<br>" << svgFile.rdbuf();
    svgFile.close();

    // Remove the .svg file.
    filesystem::remove(svgFileName);

}



// Extract  from the request the parameters for the display
// of the local assembly graph.
void Assembler::getLocalAssemblyGraphRequestParameters(
    const vector<string>& request,
    LocalAssemblyGraphRequestParameters& parameters) const
{
    parameters.edgeId = 0;
    parameters.edgeIdIsPresent = getParameterValue(
        request, "edgeId", parameters.edgeId);

    parameters.maxDistance = 0;
    parameters.maxDistanceIsPresent = getParameterValue(
        request, "maxDistance", parameters.maxDistance);

    string useDotLayoutString;
    parameters.useDotLayout = getParameterValue(
        request, "useDotLayout", useDotLayoutString);

    string showVertexLabelsString;
    parameters.showVertexLabels = getParameterValue(
        request, "showVertexLabels", showVertexLabelsString);

    string showEdgeLabelsString;
    parameters.showEdgeLabels = getParameterValue(
        request, "showEdgeLabels", showEdgeLabelsString);

    parameters.sizePixels = 600;
    parameters.sizePixelsIsPresent = getParameterValue(
        request, "sizePixels", parameters.sizePixels);

    parameters.timeout = 30;
    parameters.timeoutIsPresent = getParameterValue(
        request, "timeout", parameters.timeout);

    parameters.highlightedReadIdIsPresent = getParameterValue(
        request, "highlightedReadId", parameters.highlightedReadId);

    parameters.highlightedStrand = 0;
    getParameterValue(request, "highlightedStrand", parameters.highlightedStrand);

}



void Assembler::LocalAssemblyGraphRequestParameters::writeForm(
    ostream& html,
    AssemblyGraph::EdgeId edgeCount,
    bool allowHighlighting) const
{
    html <<
        "<h3>Display a local subgraph of the global assembly graph</h3>"
        "<form>"

        "<table>"

        "<tr title='Edge id between 0 and " << edgeCount << "'>"
        "<td>Edge id"
        "<td><input type=text required name=edgeId size=8 style='text-align:center'"
        << (edgeIdIsPresent ? ("value='"+to_string(edgeId)+"'") : "") <<
        ">"

        "<tr title='Maximum distance from start edge (number of edges)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        << (maxDistanceIsPresent ? ("value='" + to_string(maxDistance)+"'") : " value='6'") <<
        ">"

        "<tr title='Check to use graphviz dot layout. Uncheck to use graphviz sfdp layout.'>"
        "<td>Use dot layout"
        "<td class=centered><input type=checkbox name=useDotLayout"
        << (useDotLayout ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to show vertex labels'>"
        "<td>Label vertices"
        "<td class=centered><input type=checkbox name=showVertexLabels"
        << (showVertexLabels ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to show edge labels'>"
        "<td>Label edges"
        "<td class=centered><input type=checkbox name=showEdgeLabels"
        << (showEdgeLabels ? " checked=checked" : "") <<
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
        ">";

    if(allowHighlighting) {
        html <<
            "<tr>"
            "<td>Oriented read to be highlighted"
            "<td><input type=text name=highlightedReadId size=8 style='text-align:center'";

        if(highlightedReadIdIsPresent) {
            html << " value='" << highlightedReadId << "'";
        }
        html << ">";

        html << "<tr><td>Strand of oriented read to be highlighted><td class=centered>";
        writeStrandSelection(html, "highlightedStrand", highlightedStrand==0, highlightedStrand==1);


    }

    html <<
        "</table>"

        "<br><input type=submit value='Display'>"
        "</form>";
}



bool Assembler::LocalAssemblyGraphRequestParameters::hasMissingRequiredParameters() const
{
    return
        !edgeIdIsPresent ||
        !maxDistanceIsPresent ||
        !timeoutIsPresent;
}



void Assembler::exploreAssemblyGraphEdge(const vector<string>& request, ostream& html)
{
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    html << "<h2>Show information about an edge of the assembly graph</h2>";



    // Get the request parameters.

    AssemblyGraph::EdgeId edgeId = 0;
    const bool edgeIdIsPresent = getParameterValue(request, "edgeId", edgeId);

    string showSequenceString;
    getParameterValue(request, "showSequence", showSequenceString);
    const bool showSequence = (showSequenceString == "on");

    string showDetailsString;
    getParameterValue(request, "showDetails", showDetailsString);
    const bool showDetails = (showDetailsString == "on");

    uint32_t begin = 0;
    const uint32_t beginIsPresent = getParameterValue(request, "begin", begin);

    uint32_t end = 0;
    const uint32_t endIsPresent = getParameterValue(request, "end", end);



    // Write the form.
    html <<
        "<form><table>"

        "<tr><td>Assembly graph edge id<td class=centered><input type=text name=edgeId required "
        "style='text-align:center'" <<
        (edgeIdIsPresent ? (" value='" + to_string(edgeId)) + "'" : "") <<
        " title='Enter an assembly graph edge id between 0 and " << assemblyGraph.edges.size()-1 << " inclusive'"
        ">"

        "<tr><td>Show sequence<td class=centered><input type=checkbox name=showSequence" <<
        (showSequence ? " checked=checked" : "") << ">"

        "<tr><td>Show assembly details<td class=centered><input type=checkbox name=showDetails" <<
        (showDetails ? " checked=checked" : "") << ">"

        "<tr><td>Begin position in raw sequence<td class=centered><input type=text name=begin "
        "style='text-align:center'" <<
        (beginIsPresent ? (" value='" + to_string(begin)) + "'" : "") << ">"

        "<tr><td>End position in raw sequence<td class=centered><input type=text name=end "
        "style='text-align:center'" <<
        (endIsPresent ? (" value='" + to_string(end)) + "'" : "") << ">"

        "</table><br><input type=submit value='Go'>"
        "</form>";



    // If the edge id is missing or invalid, don't do anything.
    if(!edgeIdIsPresent) {
        return;
    }
    if(edgeId >= assemblyGraph.edges.size()) {
        html <<
            "<p>Invalid edge id " << edgeId <<
            ". Enter an assembly graph edge id between 0 and " <<
            assemblyGraph.edges.size()-1 << " inclusive." << endl;
        return;
    }



    // Phasing information.
    if (phasingData.orientedReads.isOpen()) {
        const span<OrientedReadId> orientedReadIds =
            phasingData.orientedReads[edgeId];
        html << "<p>The following oriented reads are internal to the this "
            "assembly graph edge:<br>";
        for(const OrientedReadId orientedReadId: orientedReadIds) {
            html << orientedReadId << " ";
        }
    }




    // If this edge was not assembled, point to the reverse
    // complemented edge, which was assembled.
    if(!assemblyGraph.isAssembledEdge(edgeId)) {
        const AssemblyGraph::EdgeId edgeIdRc = assemblyGraph.reverseComplementEdge[edgeId];
        SHASTA_ASSERT(edgeIdRc != edgeId);
        SHASTA_ASSERT(assemblyGraph.isAssembledEdge(edgeIdRc));
        html <<
            "<h1>Assembly graph edge <a href="
            "'exploreAssemblyGraph?edgeId=" << edgeId <<
            "&maxDistance=6&detailed=on&sizePixels=1600&timeout=30'>" <<
            edgeId << "</a></h1>"
            "<p>This edge was not assembled. Its reverse complement is"
            " assembly graph edge "
            "<a href='exploreAssemblyGraphEdge?edgeId=" << edgeIdRc << "'>" << edgeIdRc <<
            "</a>, which was assembled.";
        return;
    }



    // If showSequence or showDetails are selected,
    // begin and end are required.
    if(showSequence || showDetails) {
        if(!(beginIsPresent && endIsPresent)) {
            html << "<p>Specify begin and end position in raw sequence in the form above.";
            return;
        }
    }



    // Assemble the sequence and output the requested information to html.
    AssembledSegment assembledSegment;
    assembleAssemblyGraphEdge(edgeId, false, assembledSegment);
    assembledSegment.writeHtml(html, showSequence, showDetails, begin, end);
}



void Assembler::exploreAssemblyGraphEdgesSupport(
    const vector<string>& request,
    ostream& html)
{
    AssemblyGraph& assemblyGraph = *assemblyGraphPointer;

    html << "<h2>Display read support for assembly graph edges</h2>";

    // Get the request parameters.
    string edgesString;
    getParameterValue(request, "edges", edgesString);
    uint32_t beginEndMarkerCount = 100;
    getParameterValue(request, "beginEndMarkerCount", beginEndMarkerCount);



    // Parse the edge ids into tokens.
    vector<string> edgesTokens;
    boost::algorithm::split(edgesTokens, edgesString,
        boost::algorithm::is_any_of(" "),
        boost::algorithm::token_compress_on);




    // Write the form.
    html <<
        "<form><table style='width:600px'>"

        "<tr>"
        "<td>Assembly graph edges, "
        "space separated and optionally followed by B or E for begin/end only display, "
        "for example \"45 27E 101B\""
        "<td class=centered>"
        "<input type=text name=edges "
        "style='text-align:center' value='" << edgesString << "'>"

        "<tr><td>Marker count for begin/end "
        "(number of markers for begin/end only display)"
        "<td class=centered><input type=text name=beginEndMarkerCount " <<
        "style='text-align:center' value=" << beginEndMarkerCount << ">"


        "</table><br><input type=submit value='Go'>"
        "</form>";



    // Parse the requested edges to create a list of
    // (edge id, markerBegin, marker end).
    // Here, (markerBegin, marker end) are indexes into
    // marker graph edges for each assembly graph edge.
    vector< tuple<AssemblyGraph::EdgeId, uint32_t, uint32_t> > edges;
    vector<string> edgesStrings;
    for(const string& token: edgesTokens) {
        if(token.size() ==0) {
            continue;
        }

        // See if begin/end only was requested.
        bool beginOnly = false;
        bool endOnly = false;
        const char lastCharacter = token[token.size()-1];
        string edgeString = token;
        if(lastCharacter == 'B') {
            beginOnly = true;
            edgeString = edgeString.substr(0, token.size()-1);
        }
        if(lastCharacter == 'E') {
            endOnly = true;
            edgeString = edgeString.substr(0, token.size()-1);
        }

        // Extract the edge id.
        AssemblyGraph::EdgeId edgeId;
        try {
            edgeId = boost::lexical_cast<AssemblyGraph::EdgeId>(edgeString);
        } catch(const boost::bad_lexical_cast&) {
            html << "<br>Invalid assembly graph edge id " << token;
            return;
        }

        // Check that it is a valid assembly graph edge id.
        if(edgeId >= assemblyGraph.edges.size()) {
            html << "<br>Invalid assembly graph edge id " << edgeId <<
                ". Valid assembly graph edge ids are between 0 and " <<
                assemblyGraph.edges.size()-1 << " included.";
            return;
        }
        edgesStrings.push_back(token);


        // Access the marker graph edges corresponding to the assembly graph edge.
        const span<MarkerGraph::EdgeId> markerGraphEdges =
            assemblyGraph.edgeLists[edgeId];

        // Construct the requested marker interval.
        uint32_t begin = 0;
        uint32_t end = uint32_t(markerGraphEdges.size());
        if(beginOnly) {
            end = min(end, beginEndMarkerCount);
        } else if(endOnly) {
            if(end > beginEndMarkerCount) {
                begin = end - beginEndMarkerCount;
            }
        }
        edges.push_back(make_tuple(edgeId, begin, end));
    }
    SHASTA_ASSERT(edges.size() == edgesStrings.size());



    // Gather marker graph vertex ids.
    vector< vector<MarkerGraph::VertexId> > markerGraphVertexIds(edges.size());
    for(size_t i=0; i<edges.size(); i++)  {
        const auto& t = edges[i];
        const AssemblyGraph::EdgeId assemblyGraphEdgeId = std::get<0>(t);
        const span<MarkerGraph::EdgeId> markerGraphEdgeIds =
            assemblyGraph.edgeLists[assemblyGraphEdgeId];
        const uint32_t begin = std::get<1>(t);
        const uint32_t end = std::get<2>(t);
        for(size_t j=begin; j!=end; j++) {
            const MarkerGraph::EdgeId& markerGraphEdgeId = markerGraphEdgeIds[j];
            const MarkerGraph::Edge& markerGraphEdge = markerGraph.edges[markerGraphEdgeId];
            if(j==begin && j!=0) {  // Never include the very first vertex!
                markerGraphVertexIds[i].push_back(markerGraphEdge.source);
            }
            if(j != markerGraphEdgeIds.size()-1) {    // Never include the very last vertex!
                markerGraphVertexIds[i].push_back(markerGraphEdge.target);
            }
        }
    }



    // Gather the oriented read ids represented
    // in these vertices.
    std::set<OrientedReadId> orientedReadIdsSet;

    // Loop over the requested assembly graph edges.
    for(const vector<MarkerGraph::VertexId>& v: markerGraphVertexIds) {

        // Loop over marker graph vertices in the requested interval
        // for this assembly graph edge.
        for(MarkerGraph::VertexId vertexId: v) {

            // Access the marker ids on this vertex.
            const span<MarkerId> markerIds = markerGraph.vertices[vertexId];

            // Loop over these markers.
            for(const MarkerId markerId: markerIds) {
                OrientedReadId orientedReadId;
                tie(orientedReadId, ignore) = findMarkerId(markerId);
                orientedReadIdsSet.insert(orientedReadId);
            }

        }
    }
    vector<OrientedReadId> orientedReadIds(
        orientedReadIdsSet.begin(),
        orientedReadIdsSet.end());



    // Find out at what positions (in markerGraphVertexIds) each
    // oriented read id appears in each of the requested assembly graph edges.
    // The table is indexed by [oriented read id index][assembly graph edge index],
    // where:
    // - Oriented read id index = index into orientedReadids vector.
    // - Assembly graph edge index = index into markerGraphVertexIds vector
    vector< vector< vector<uint32_t> > > table(
        orientedReadIds.size(),
        vector< vector<uint32_t> >(edges.size()));

    // Loop over requested assembly graph edges.
    for(size_t assemblyGraphEdgeIndex=0;
        assemblyGraphEdgeIndex<markerGraphVertexIds.size();
        assemblyGraphEdgeIndex++) {

        // Loop over the marker grah vertices.
        const vector<MarkerGraph::VertexId>& v = markerGraphVertexIds[assemblyGraphEdgeIndex];
        for(size_t iv=0; iv<v.size(); iv++) {
            const MarkerGraph::VertexId vertexId = v[iv];

            // Access the marker ids on this vertex.
            const span<MarkerId> markerIds = markerGraph.vertices[vertexId];

            // Loop over these markers.
            for(const MarkerId markerId: markerIds) {
                OrientedReadId orientedReadId;
                tie(orientedReadId, ignore) = findMarkerId(markerId);
                const auto it = std::lower_bound(
                    orientedReadIds.begin(), orientedReadIds.end(), orientedReadId);
                SHASTA_ASSERT(it != orientedReadIds.end());
                SHASTA_ASSERT(*it == orientedReadId);
                const size_t orientedReadIdIndex = it - orientedReadIds.begin();

                table[orientedReadIdIndex][assemblyGraphEdgeIndex].push_back(uint32_t(iv));
            }
        }

    }



    // Write a table with the results.
    html << "<p><table>";

    // Assembly graph edge ids.
    html << "<tr><th class=left>Assembly graph edge id";
    for(const string& edgeString : edgesStrings) {
        html << "<th class=centered>" << edgeString;
    }

    // Total number of marker graph edges.
    html << "<tr><th class=left>Total marker graph edge count";
    for(const auto& t: edges) {
        html << "<td class=centered>" << assemblyGraph.edgeLists.size(std::get<0>(t));
    }

    // Marker graph edge requested begin/end.
    html << "<tr><th class=left>Requested marker graph edge begin";
    for(const auto& t: edges) {
        html << "<td class=centered>" << std::get<1>(t);
    }
    html << "<tr><th class=left>Requested marker graph edge end";
    for(const auto& t: edges) {
        html << "<td class=centered>" << std::get<2>(t);
    }
    html << "<tr><th class=left>Requested marker graph edge count";
    for(const auto& t: edges) {
        html << "<td class=centered>" <<
            std::get<2>(t) - std::get<1>(t);
    }



    // One row for each oriented read id.
    for(size_t orientedReadIdIndex=0;
        orientedReadIdIndex< orientedReadIds.size();
        orientedReadIdIndex++) {
        const OrientedReadId orientedReadId = orientedReadIds[orientedReadIdIndex];
        html << "<tr><td>" << orientedReadId;

        // Loop over assembly graph edges.
        for(size_t i=0; i<edges.size(); i++) {
            const uint32_t markerGraphVertexCount = uint32_t(markerGraphVertexIds[i].size());
            const string canvasId = to_string(orientedReadId.getValue()) + "-" + to_string(i);
            html <<
                "<td class=centered>"
                "<canvas id='" << canvasId << "'"
                " width=" << markerGraphVertexCount << " height=10>"
                "</canvas>"
                "<script>"
                "var c = document.getElementById('" << canvasId << "');"
                "var ctx = c.getContext('2d');"
                "ctx.fillStyle = '#e0e0e0';"
                "ctx.fillRect(0, 0, " << markerGraphVertexCount << ", 10);"
                "ctx.fillStyle = '#ff0000';";

            for(const uint32_t iv: table[orientedReadIdIndex][i]) {
                html << "ctx.fillRect(" << iv << ", 0, 1, 10);";
            }

            html << "</script>";
        }
    }


    html << "</table>";

}



#endif
