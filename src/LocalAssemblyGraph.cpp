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



// Return the "length" of a vertex, that is, the
// number of marker graph edges that it corresponds to.
size_t LocalAssemblyGraph::vertexLength(vertex_descriptor v) const
{
    const LocalAssemblyGraph& graph = *this;
    const VertexId vertexId = graph[v].vertexId;
    return globalAssemblyGraph.vertices.size(vertexId);
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
    // This turns off the tooltip on the graph.
    s << "tooltip = \" \";\n";

    if(detailed) {
        s << "layout=dot;\n";
        s << "rankdir=LR;\n";
        s << "ratio=expand;\n";
        s << "node [fontname=\"Courier New\" shape=rectangle style=filled];\n";
        s << "edge [fontname=\"Courier New\"];\n";
    } else {
        s << "layout=sfdp;\n";
        s << "smoothing=triangle;\n";
        s << "ratio=expand;\n";
        s << "node [shape=point];\n";
    }
}



// Write a vertex in graphviz format.
void LocalAssemblyGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const LocalAssemblyGraphVertex& vertex = graph[v];
    const size_t length = graph.vertexLength(v);

    // Begin vertex attributes.
    s << "[";



    // Color.
    string color;
    if(vertex.distance == maxDistance) {
        // Vertices at maximum distance.
        color = "cyan";
    } else if(vertex.distance == 0) {
        // Start vertex.
        color = "#90ee90";
    } else {
        // All other vertices.
        if(detailed) {
            color = "pink";
        } else {
            color = "black";
        }
    }
    s << " fillcolor=\"" << color << "\"";
    if(!detailed) {
        s << " color=\"" << color << "\"";
    }



    // Size.
    // This could be problematic for the compressed assembly graph.
    /*
    if(detailed) {
        s << " width=" << 1. * double(length);
    } else {
        s << " width=" << 1.e-1 * double(length);
    }
    */



    // Toolip.
    if(!detailed) {
        s << " tooltip=\"Id " << vertex.vertexId;
        s << ", length " << length << "\"";
    } else {
        s << " tooltip=\" \"";
    }



    // Label.
    if(detailed) {
        s << " label=\"Id " << vertex.vertexId << "\\n";
        s << "Length " << length << "\"";
    }


    // Link to detailed information for this vertex.
    s << " URL=\"exploreAssemblyGraphVertex?vertexId=" << vertex.vertexId << "\"";


    // End vertex attributes.
    s << "]";
}



void LocalAssemblyGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
}
