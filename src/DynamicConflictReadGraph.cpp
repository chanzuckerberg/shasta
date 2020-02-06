#include "DynamicConflictReadGraph.hpp"
using namespace shasta;


DynamicConflictReadGraph::DynamicConflictReadGraph(const ConflictReadGraph& conflictReadGraph)
{
    DynamicConflictReadGraph& graph = *this;

    // Create the vertices.
    vector<vertex_descriptor> vertexTable;
    for(VertexId vertexId=0; vertexId<conflictReadGraph.vertices.size(); vertexId++) {
        const vertex_descriptor v = boost::add_vertex(DynamicConflictReadGraphVertex(vertexId), graph);
        vertexTable.push_back(v);
    }

    // Create the edges.
    for(EdgeId edgeId=0; edgeId<conflictReadGraph.edges.size(); edgeId++) {
        const VertexId vertexId0 = conflictReadGraph.v0(edgeId);
        const VertexId vertexId1 = conflictReadGraph.v1(edgeId);
        const vertex_descriptor v0 = vertexTable[vertexId0];
        const vertex_descriptor v1 = vertexTable[vertexId1];
        boost::add_edge(v0, v1, DynamicConflictReadGraphEdge(edgeId), graph);
    }
}

