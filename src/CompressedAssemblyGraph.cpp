#include "CompressedAssemblyGraph.hpp"
#include "findLinearChains.hpp"
#include "html.hpp"
using namespace shasta;

#include "fstream.hpp"
#include "vector.hpp"



// Create the CompressedAssemblyGraph from the AssemblyGraph.
CompressedAssemblyGraph::CompressedAssemblyGraph(
    const AssemblyGraph& assemblyGraph)
{
    CompressedAssemblyGraph& graph = *this;

    cout << "The assembly graph has " << assemblyGraph.vertices.size() <<
        " vertices and " << assemblyGraph.edges.size() << " edges." << endl;

    // Create a vertex for each vertex of the assembly graph.
    vector<vertex_descriptor> vertexTable;
    createVertices(assemblyGraph.vertices.size(), vertexTable);

    // Create an edge for each set of parallel edges of the assembly graph.
    createEdges(assemblyGraph, vertexTable);

    // Merge linear chains of edges.
    mergeLinearChains();

    cout << "The compressed assembly graph has " <<
        num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;

    // Assign an id to each edge.
    assignEdgeIds();

    // Fill in the assembly graph edges that go into each
    // edge of the compressed assembly graph.
    fillContributingEdges(assemblyGraph);

    // Fill in minimum and maximum marker counts for each edge.
    fillMarkerCounts(assemblyGraph);
}



// Create a vertex for each vertex of the assembly graph.
void CompressedAssemblyGraph::createVertices(
    uint64_t vertexCount,
    vector<vertex_descriptor>& vertexTable)
{
    CompressedAssemblyGraph& graph = *this;

    vertexTable.resize(vertexCount, null_vertex());

    for(VertexId vertexId=0; vertexId<vertexCount; vertexId++) {
        const vertex_descriptor v = add_vertex(CompressedAssemblyGraphVertex(vertexId), graph);
        vertexTable[vertexId] = v;
    }
}



// Create an edge for each set of parallel edges of the assembly graph.
void CompressedAssemblyGraph::createEdges(
    const AssemblyGraph& assemblyGraph,
    const vector<vertex_descriptor>& vertexTable)
{
    CompressedAssemblyGraph& graph = *this;

    // Loop over assembly graph edges.
    for(const AssemblyGraph::Edge& edge: assemblyGraph.edges) {
        const vertex_descriptor v0 = vertexTable[edge.source];
        const vertex_descriptor v1 = vertexTable[edge.target];

        // Do we already have an edge between these two vertices?
        bool edgeExists = false;
        tie(ignore, edgeExists) = boost::edge(v0, v1, graph);

        // If we don't already have an edge between these two vertices, add it.
        // We only create one of each set of parallel edges.
        if(not edgeExists) {
            edge_descriptor e;
            bool edgeWasAdded = false;
            tie(e, edgeWasAdded) = add_edge(v0, v1, graph);
            SHASTA_ASSERT(edgeWasAdded);

            // Store the assembly graph vertices.
            CompressedAssemblyGraphEdge& edge = graph[e];
            edge.vertices.push_back(graph[v0].vertexId);
            edge.vertices.push_back(graph[v1].vertexId);
        }
    }

}



// Merge linear chains of edges.
void CompressedAssemblyGraph::mergeLinearChains()
{
    CompressedAssemblyGraph& graph = *this;

    // Find linear chains.
    vector< std::list<edge_descriptor> > chains;
    findLinearChains(graph, chains);



    // Replace each chain with a single edge.
    for(const std::list<edge_descriptor>& chain: chains) {

        // If the chain has length 1, leave it alone.
        if(chain.size() == 1) {
            continue;
        }

        // Add the new edge.
        const vertex_descriptor v0 = source(chain.front(), graph);
        const vertex_descriptor v1 = target(chain.back(), graph);
        edge_descriptor eNew;
        bool edgeWasAdded = false;
        tie(eNew, edgeWasAdded) = add_edge(v0, v1, graph);
        SHASTA_ASSERT(edgeWasAdded);
        CompressedAssemblyGraphEdge& newEdge = graph[eNew];

        // Fill in the assembly graph vertices corresponding to this new edge.
        newEdge.vertices.push_back(graph[v0].vertexId);
        for(const edge_descriptor e: chain) {
            const vertex_descriptor v = target(e, graph);
            newEdge.vertices.push_back(graph[v].vertexId);
        }

        // Remove the edges of the chain.
        // We will remove the vertices later.
        for(const edge_descriptor e: chain) {
            boost::remove_edge(e, graph);
        }
    }



    // Remove the vertices that have become isolated.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, CompressedAssemblyGraph) {
        if(in_degree(v, graph) == 0 and out_degree(v, graph) == 0) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        remove_vertex(v, graph);
    }

}



// Assign an id to each edge.
void CompressedAssemblyGraph::assignEdgeIds()
{
    CompressedAssemblyGraph& graph = *this;

    uint64_t edgeId = 0;
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        graph[e].id = edgeId++;
    }

}



// Fill in the assembly graph edges that go into each
// edge of the compressed assembly graph.
void CompressedAssemblyGraph::fillContributingEdges(
    const AssemblyGraph& assemblyGraph)
{
    CompressedAssemblyGraph& graph = *this;

    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        CompressedAssemblyGraphEdge& edge = graph[e];
        edge.edges.resize(edge.vertices.size() - 1);
        for(uint64_t i=0; i<edge.edges.size(); i++) {
            const VertexId vertexId0 = edge.vertices[i];
            const VertexId vertexId1 = edge.vertices[i+1];
            const span<const EdgeId> edges0 = assemblyGraph.edgesBySource[vertexId0];
            for(const EdgeId edge01: edges0) {
                if(assemblyGraph.edges[edge01].target == vertexId1) {
                    edge.edges[i].push_back(edge01);
                }
            }
        }

    }
}



string CompressedAssemblyGraphEdge::gfaId() const
{
    if(edges.size()==1 and edges.front().size()==1) {
        // Return the one and only assembly graph edge associated with this
        // compressed assembly graph edge.
        return to_string(edges.front().front());
    } else {
        // Return the id of this compressed assembly graph edge,
        // prefixed with "C".
        return "C" + to_string(id);
    }
}



// GFA output (without sequence).
void CompressedAssemblyGraph::writeGfa(const string& fileName, double basesPerMarker) const
{
    ofstream gfa(fileName);
    writeGfa(gfa, basesPerMarker);
}
void CompressedAssemblyGraph::writeGfa(ostream& gfa, double basesPerMarker) const
{
    const CompressedAssemblyGraph& graph = *this;

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";


    // Write a segment record for each edge.
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        const CompressedAssemblyGraphEdge& edge = graph[e];

        gfa <<
            "S\t" <<
            edge.gfaId() << "\t" <<
            "*\t" <<
            "LN:i:" << uint64_t(basesPerMarker * 0.5 *
                double(edge.minMarkerCount + edge.maxMarkerCount)) <<
            "\n";
    }


    // Write GFA links.
    // For each vertex in the compressed assembly graph there is a link for
    // each combination of in-edges and out-edges.
    // Therefore each vertex generates a number of
    // links equal to the product of its in-degree and out-degree.
    BGL_FORALL_VERTICES(v, graph, CompressedAssemblyGraph) {
        BGL_FORALL_INEDGES(v, eIn, graph, CompressedAssemblyGraph) {
            BGL_FORALL_OUTEDGES(v, eOut, graph, CompressedAssemblyGraph) {
                gfa <<
                    "L\t" <<
                    graph[eIn].gfaId() << "\t" <<
                    "+\t" <<
                    graph[eOut].gfaId() << "\t" <<
                    "+\t" <<
                    "*\n";
            }
        }
    }
}


// HTML output.
void CompressedAssemblyGraph::writeHtml(const string& fileName) const
{
    ofstream html(fileName);
    writeHtml(html);
}
void CompressedAssemblyGraph::writeHtml(ostream& html) const
{
    writeHtmlBegin(html, "Compressed assembly graph");
    html <<
        "<body><h1>Compressed assembly graph</h1>"
        "<p>Each edge of the compressed assembly graph corresponds to either "
        "a single edge of uncompressed assembly graph, "
        "or a chain of bubbles in the uncompressed assembly graph. "
        "The following table summarizes that uncompressed assembly graph edges "
        "that contribute to each edge of the compressed assembly graph."
        "<table>"
        "<tr><th>Compressed<br>edge<th>Position<th>Uncompressed<br>edges\n";

    const CompressedAssemblyGraph& graph = *this;

    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        const CompressedAssemblyGraphEdge& edge = graph[e];


        for(uint64_t position=0; position<edge.edges.size(); position++) {
            const vector<AssemblyGraph::EdgeId>& edgesAtPosition = edge.edges[position];
            html << "<tr><td class=centered title='Compressed edge'>" << edge.gfaId();
            html << "<td class=centered title='Position in compressed edge'>" <<
                position << "<td class=centered title='Uncompressed edges'>";
            for(const AssemblyGraph::EdgeId edgeId: edgesAtPosition) {
                html << " " << edgeId;
            }
            html << "\n";
        }
    }

    html << "</table></body>";
    writeHtmlEnd(html);
}


// Fill in minimum and maximum marker counts for each edge.
void CompressedAssemblyGraph::fillMarkerCounts(const AssemblyGraph& assemblyGraph)
{
    CompressedAssemblyGraph& graph = *this;

    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        graph[e].fillMarkerCounts(assemblyGraph);
    }
}
void CompressedAssemblyGraphEdge::fillMarkerCounts(const AssemblyGraph& assemblyGraph)
{
    minMarkerCount = 0;
    maxMarkerCount = 0;
    for(const vector<AssemblyGraph::EdgeId>& parallelEdges: edges) {
        SHASTA_ASSERT(not parallelEdges.empty());

        // Compute the minimum and maximum number of markers
        // over this set of parallel edges.
        uint64_t minMarkerCountHere = std::numeric_limits<uint64_t>::max();
        uint64_t maxMarkerCountHere = 0;
        for(const AssemblyGraph::EdgeId edgeId: parallelEdges) {
            const uint64_t markerCount = assemblyGraph.edgeLists.size(edgeId);
            minMarkerCountHere = min(minMarkerCountHere, markerCount);
            maxMarkerCountHere = max(maxMarkerCountHere, markerCount);
        }

        // Update the totals.
        minMarkerCount += minMarkerCountHere;
        maxMarkerCount += maxMarkerCountHere;
    }
}
