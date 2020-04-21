#include "Assembler.hpp"
#include "CompressedAssemblyGraph.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>

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
    cout << "The local subgraph has " << num_vertices(subgraph) <<
        " vertices and " << num_edges(subgraph) << " edges." << endl;

    cout << "Related edges:" << endl;
    BGL_FORALL_EDGES(e0, subgraph, CompressedAssemblyGraph) {
        const auto& edge0 = subgraph[e0];
        cout << edge0.gfaId() << ":";
        for(const edge_descriptor e1: edge0.relatedEdges) {
            cout << " " << subgraph[e1].gfaId();
        }
        cout << endl;
    }


}
