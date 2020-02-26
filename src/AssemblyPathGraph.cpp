// Shasta.
#include "AssemblyPathGraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard library.
#include "fstream.hpp"
#include <set>



AssemblyPathGraph::AssemblyPathGraph(const AssemblyGraph& assemblyGraph)
{
    AssemblyPathGraph& graph = *this;

    // Create a vertex for each assembly graph vertex.
    vector<vertex_descriptor> vertices;
    for(AssemblyGraph::VertexId vertexId=0; vertexId<assemblyGraph.vertices.size(); vertexId++) {
        const vertex_descriptor v = add_vertex(AssemblyPathGraphVertex(vertexId), graph);
        vertices.push_back(v);
    }



    // Create an edge for each assembly graph edge.
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[edgeId];
        const AssemblyGraph::VertexId vertexId0 = edge.source;
        const AssemblyGraph::VertexId vertexId1 = edge.target;
        const vertex_descriptor v0 = vertices[vertexId0];
        const vertex_descriptor v1 = vertices[vertexId1];
        edge_descriptor e;
        tie(e, ignore) = add_edge(v0, v1, AssemblyPathGraphEdge(edgeId), graph);
    }

}



void AssemblyPathGraph::writeGraphviz(const string& fileName) const
{
    ofstream file(fileName);
    writeGraphviz(file);
}



void AssemblyPathGraph::writeGraphviz(ostream& s) const
{
    const AssemblyPathGraph& graph = *this;

    s << "digraph G {\n";

    // Default attributes.
    s << "layout=sfdp;\n";
    s << "K=10;\n";
    s << "overlap=false;\n";
    s << "splines=true;\n";
    s << "smoothing=triangle;\n";
    s << "node [shape=point];\n";

    // This turns off the tooltip on the graph and the edges.
    s << "tooltip = \" \";\n";



    // Vertices.
    BGL_FORALL_VERTICES(v, graph, AssemblyPathGraph) {
        const AssemblyPathGraphVertex& vertex = graph[v];
        s << vertex.vertexId;

        s << " [";

        s << "tooltip=\"" << vertex.vertexId << "\"";

        s << "]";

        s << "\n";
    }



    // Edges. We write each edge as an additional pseudovertex.
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph) {
        const AssemblyPathGraphEdge& edge = graph[e];

        // Get the vertices.
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const AssemblyPathGraphVertex& vertex0 = graph[v0];
        const AssemblyPathGraphVertex& vertex1 = graph[v1];

        const bool isTangle = in_degree(v0, graph) > 1 and out_degree(v1, graph) > 1;

        // Write is as a pseudo vertex.
        const string pseudoVertexName =
            "\"" +
            to_string(vertex0.vertexId) +
            "to" +
            to_string(vertex1.vertexId) +
            "\"";
        s << pseudoVertexName << " [";
        s << "shape=rectangle label=\"";
        for(uint64_t i=0; i<edge.path.size(); i++) {
            AssemblyGraph::EdgeId edgeId = edge.path[i];
            s << edgeId;
            if(i != edge.path.size()-1) {
                s << " ";
            }
            s << "\\n" << edge.pathLength;
        }
        s << "\"";
        if(isTangle) {
            s << " style=filled fillcolor=pink";
        }
        s << "];\n";

            // Write the arrows to/from the pseudovertex.
        s << vertex0.vertexId << "->" << pseudoVertexName << ";\n";
        s << pseudoVertexName << "->" << vertex1.vertexId << ";\n";
    }





    s << "}\n";
}



