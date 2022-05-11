// Shasta.
#include "Assembler.hpp"
#include "CompressedAssemblyGraph.hpp"
#include "platformDependent.hpp"
#include "runCommandWithTimeout.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <map>

void Assembler::exploreCompressedAssemblyGraph(
    const vector<string>& request,
    ostream& html)
{
    using vertex_descriptor = CompressedAssemblyGraph::vertex_descriptor;
    using edge_descriptor = CompressedAssemblyGraph::edge_descriptor;

    html << "<h1>Compressed assembly graph</h1>"
        "<p>Each edge of the compressed assembly graph corresponds to "
        " a linear sequence of bubbles in the assembly graph, "
        "which can degenerate into a single edge.";

    // Create the CompressedAssemblyGraph, if necessary.
    if(not compressedAssemblyGraph) {
        html << "<pre>" << timestamp << "Creating the compressed assembly graph.\n";
        createCompressedAssemblyGraph();
        html << timestamp << "Done creating the compressed assembly graph.</pre>";
    }



    // Get the parameters.
    string startEdgeGfaId;
    getParameterValue(request, "startEdgeGfaId", startEdgeGfaId);

    uint32_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    uint32_t sizePixels = 600;
    getParameterValue(request, "sizePixels", sizePixels);

    double vertexScalingFactor = 0.1;
    getParameterValue(request, "vertexScalingFactor", vertexScalingFactor);

    double edgeLengthPower = 0.25;
    getParameterValue(request, "edgeLengthPower", edgeLengthPower);

    double edgeLengthScalingFactor = 1.;
    getParameterValue(request, "edgeLengthScalingFactor", edgeLengthScalingFactor);

    double edgeThicknessScalingFactor = 3.;
    getParameterValue(request, "edgeThicknessScalingFactor", edgeThicknessScalingFactor);

    double edgeArrowScalingFactor = 1.;
    getParameterValue(request, "edgeArrowScalingFactor", edgeArrowScalingFactor);

    double timeout= 30;
    getParameterValue(request, "timeout", timeout);



    // Write the form.
    html <<
        "<form><table>" <<

        "<tr><td>GFA id of the start edge"
        "<td><input type=text required name=startEdgeGfaId size=8 style='text-align:center'" <<
        (startEdgeGfaId.empty() ? "" :  (" value='" + startEdgeGfaId + "'")) << ">"

        "<tr><td>Maximum distance for local subgraph (in edges)"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance << "'>"

        "<tr><td>Graphics size in pixels"
        "<td><input type=text required name=sizePixels size=8 style='text-align:center'"
        " value='" << sizePixels << "'>"

        "<tr><td>Vertex scaling factor"
        "<td><input type=text required name=vertexScalingFactor size=8 style='text-align:center'"
        " value='" << vertexScalingFactor << "'>"

        "<tr><td>Edge length power"
        "<td><input type=text required name=edgeLengthPower size=8 style='text-align:center'"
        " value='" << edgeLengthPower << "'>"

        "<tr><td>Edge length scaling factor"
        "<td><input type=text required name=edgeLengthScalingFactor size=8 style='text-align:center'"
        " value='" << edgeLengthScalingFactor << "'>"

        "<tr><td>Edge thickness scaling factor"
        "<td><input type=text required name=edgeThicknessScalingFactor size=8 style='text-align:center'"
        " value='" << edgeThicknessScalingFactor << "'>"

        "<tr><td>Edge arrow scaling factor"
        "<td><input type=text required name=edgeArrowScalingFactor size=8 style='text-align:center'"
        " value='" << edgeArrowScalingFactor << "'>"

        "<tr><td>Timeout for graph layout (second)"
        "<td><input type=text required name=timeout size=8 style='text-align:center'"
        " value='" << timeout << "'>"

        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";




    // Locate the start edge in the compressed graph.
    if(startEdgeGfaId.empty()) {
        return;
    }
    bool edgeWasFound = false;
    edge_descriptor eStart;
    const CompressedAssemblyGraph& graph = *compressedAssemblyGraph;
    tie(eStart, edgeWasFound) =  graph.getEdgeFromGfaId(startEdgeGfaId);
    if(not edgeWasFound) {
        html << "<p>This edge was not found. See CompressedAssemblyGraph.gfa for valid GFA ids.";
        return;
    }



    // Create the requested local subgraph.
    vector<vertex_descriptor> startVertices;
    startVertices.push_back(source(eStart, graph));
    startVertices.push_back(target(eStart, graph));
    boost::bimap<vertex_descriptor, vertex_descriptor> vertexMap;
    boost::bimap<edge_descriptor, edge_descriptor> edgeMap;
    std::map<vertex_descriptor, uint64_t> distanceMap;
    CompressedAssemblyGraph subgraph(
        graph,
        *this,
        startVertices,
        maxDistance,
        vertexMap,
        edgeMap,
        distanceMap);

    // Compute vertex layout.
    std::map<CompressedAssemblyGraph::vertex_descriptor, array<double, 2 > > vertexPositions;
    subgraph.computeVertexLayout(
        sizePixels,
        vertexScalingFactor,
        edgeLengthPower,
        edgeLengthScalingFactor,
        timeout,
        vertexPositions);

    // Write it in Graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    subgraph.writeGraphviz(
        dotFileName,
        sizePixels,
        vertexScalingFactor,
        edgeLengthScalingFactor,
        edgeThicknessScalingFactor,
        edgeArrowScalingFactor,
        vertexPositions
        );



    // Compute graph layout and write it in svg format.
    const string command = "neato -n2 -O -T svg " + dotFileName;
    bool timeoutTriggered = false;
    bool signalOccurred = false;
    int returnCode;
    runCommandWithTimeout(command, timeout,
        timeoutTriggered, signalOccurred, returnCode);
    if(signalOccurred) {
        html << "<p>Unable to compute graph layout: terminated by a signal. "
            "The failing Command was: <code>" << command << "</code>";
        return;
    }
    if(timeoutTriggered) {
        html << "<p>Timeout exceeded during graph layout computation. "
            "Increase the timeout or decrease the maximum distance to simplify the graph";
        return;
    }
    if(returnCode!=0 ) {
        html << "<p>Unable to compute graph layout: return code " << returnCode <<
            ". The failing Command was: <code>" << command << "</code>";
        return;
    }
    filesystem::remove(dotFileName);


    // Display the graph.
    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << svgFile.rdbuf();
    svgFile.close();
    filesystem::remove(svgFileName);

}
