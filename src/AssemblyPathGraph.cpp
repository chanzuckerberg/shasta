#include "AssemblyPathGraph.hpp"
using namespace shasta;

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
