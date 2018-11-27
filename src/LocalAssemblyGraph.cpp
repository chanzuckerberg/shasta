// Shasta.
#include "LocalAssemblyGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard libraries.
#include "fstream.hpp"



LocalAssemblyGraph::LocalAssemblyGraph(const AssemblyGraph& globalAssemblyGraph) :
    globalAssemblyGraph(globalAssemblyGraph)
{

}



// Add a vertex with the given VertexId
// and return its vertex descriptor.
// A vertex with this VertexId must not exist.
LocalAssemblyGraph::vertex_descriptor LocalAssemblyGraph::addVertex(
    AssemblyGraph::VertexId vertexId,
    int distance)
{
    // Check that the vertex does not already exist.
    CZI_ASSERT(vertexMap.find(vertexId) == vertexMap.end());

    // Add the vertex and store it in the vertex map.
    const vertex_descriptor v = add_vertex(LocalAssemblyGraphVertex(vertexId, distance), *this);
    vertexMap.insert(make_pair(vertexId, v));

    return v;
}



// Find out if a vertex with the given VertexId exists.
// If it exists, return make_pair(true, v).
// Otherwise, return make_pair(false, null_vertex());
std::pair<bool, LocalAssemblyGraph::vertex_descriptor>
    LocalAssemblyGraph::findVertex(AssemblyGraph::VertexId vertexId) const
{
    const auto it = vertexMap.find(vertexId);
    if(it == vertexMap.end()) {
        return make_pair(false, null_vertex());
    } else {
        const vertex_descriptor v = it->second;
        return make_pair(true, v);
    }
}



// Write the graph in Graphviz format.
void LocalAssemblyGraph::write(
    const string& fileName,
    int maxDistance,
    bool detailed) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, maxDistance, detailed);
}
void LocalAssemblyGraph::write(
    ostream& s,
    int maxDistance,
    bool detailed) const
{
    Writer writer(*this, maxDistance, detailed);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&LocalAssemblyGraphVertex::vertexId, *this));
}



LocalAssemblyGraph::Writer::Writer(
    const LocalAssemblyGraph& graph,
    int maxDistance,
    bool detailed) :
    graph(graph),
    maxDistance(maxDistance),
    detailed(detailed)
{
}



void LocalAssemblyGraph::Writer::operator()(std::ostream& s) const
{
}



void LocalAssemblyGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
}



void LocalAssemblyGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
}
