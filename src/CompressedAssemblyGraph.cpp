#include "CompressedAssemblyGraph.hpp"
using namespace shasta;

#include "vector.hpp"



// Create the CompressedAssemblyGraph from the AssemblyGraph.
CompressedAssemblyGraph::CompressedAssemblyGraph(
    const AssemblyGraph& assemblyGraph)
{
    CompressedAssemblyGraph& graph = *this;

    cout << "The assembly graph has " << assemblyGraph.vertices.size() <<
        " vertices and " << assemblyGraph.edges.size() << " edges." << endl;


    // Create a vertex for each vertex of the Assembly graph.
    vector<vertex_descriptor> vertexTable(assemblyGraph.vertices.size(), null_vertex());
    for(VertexId vertexId=0; vertexId<assemblyGraph.vertices.size(); vertexId++) {
        vertexTable[vertexId] = add_vertex(CompressedAssemblyGraphVertex(vertexId), graph);
    }

    // Create an edge for each set of parallel edges of the assembly graph.
    for(const AssemblyGraph::Edge& edge: assemblyGraph.edges) {
        const vertex_descriptor v0 = vertexTable[edge.source];
        const vertex_descriptor v1 = vertexTable[edge.target];
        bool edgeExists = false;
        edge_descriptor e;
        tie(e, edgeExists) = boost::edge(v0, v1, graph);
        if(not edgeExists) {
            add_edge(v0, v1, graph);
        }
    }
    cout << "Before compression, the compressed assembly graph has " <<
        num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;

    SHASTA_ASSERT(0);
}
