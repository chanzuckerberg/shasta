// Shasta.
#include "DynamicConflictReadGraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"



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



void DynamicConflictReadGraph::writeGraphviz(const string& fileName) const
{
    ofstream file(fileName);
    writeGraphviz(file);
}
void DynamicConflictReadGraph::writeGraphviz(ostream& s) const
{
    s << "graph G{\n";

    const auto& graph = *this;
    BGL_FORALL_VERTICES(v, graph, DynamicConflictReadGraph) {
        writeGraphviz(s, v);
    }
    BGL_FORALL_EDGES(E, graph, DynamicConflictReadGraph) {
        writeGraphviz(s, E);
    }

    s << "}\n";
}



void DynamicConflictReadGraph::writeGraphviz(ostream& s, vertex_descriptor v) const
{
    const auto& graph = *this;
    const auto& vertex = graph[v];
    const OrientedReadId orientedReadId = vertex.getOrientedReadId();

    s << "\"" << orientedReadId << "\"";
    s << "[";
    s << "tooltip=\"" << orientedReadId << "\"";
    s << "];\n";
}



void DynamicConflictReadGraph::writeGraphviz(ostream& s, edge_descriptor e) const
{
    const auto& graph = *this;
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);

    s << "\"" << graph[v0].getOrientedReadId() << "\"--\"" <<
        graph[v1].getOrientedReadId() << "\";\n";
}

