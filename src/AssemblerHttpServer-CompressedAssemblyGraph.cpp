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

    // Get the GFA if of the start edge for the local subgraph.
    string startEdgeGfaId;
    getParameterValue(request, "startEdgeGfaId", startEdgeGfaId);

    // Get the maximum distance for the local subgraph.
    uint32_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);



    // Write the form.
    html <<
        "<p>Select the start edge and maximum distance for the local subgraph to display."
        "<form><table>" <<

        "<tr><td>GFA id of the start edge"
        "<td><input type=text required name=startEdgeGfaId size=8 style='text-align:center'" <<
        (startEdgeGfaId.empty() ? "" :  (" value='" + startEdgeGfaId + "'")) << ">"

        "<tr><td>Maximum distance for local subgraph (in edges)"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance << "'>"

        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";



    // Locate the start edge in the compressed graph.
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

    // Write it in Graphviz format.
    const double edgeLengthScalingFactor = 1.e-3;
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    subgraph.writeGraphviz(dotFileName, edgeLengthScalingFactor);



    // Compute graph layout and write it in svg format.
    const string command = "dot -O -T svg " + dotFileName;
    const double timeout = 30.;
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
            "The failing Command was: <code>" << command << "</code>";
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
