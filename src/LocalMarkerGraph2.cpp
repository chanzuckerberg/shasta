// Nanopore2
#include "LocalMarkerGraph2.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard libraries.
#include "fstream.hpp"
#include "stdexcept.hpp"
#include "utility.hpp"


// Find out if a vertex with the given GlobalMarkerGraphVertexId exists.
// If it exists, return make_pair(true, v).
// Otherwise, return make_pair(false, null_vertex());
std::pair<bool, LocalMarkerGraph2::vertex_descriptor>
    LocalMarkerGraph2::findVertex(GlobalMarkerGraphVertexId vertexId) const
{
    const auto it = vertexMap.find(vertexId);
    if(it == vertexMap.end()) {
        return make_pair(false, null_vertex());
    } else {
        const vertex_descriptor v = it->second;
        return make_pair(true, v);
    }
}


// Add a vertex with the given GlobalMarkerGraphVertexId
// and return its vertex descriptor.
// A vertex with this GlobalMarkerGraphVertexId must not exist.
LocalMarkerGraph2::vertex_descriptor
    LocalMarkerGraph2::addVertex(
    GlobalMarkerGraphVertexId vertexId,
    int distance,
    MemoryAsContainer<OrientedMarkerId> markers)
{
    CZI_ASSERT(vertexMap.find(vertexId) == vertexMap.end());
    const vertex_descriptor v = add_vertex(LocalMarkerGraph2Vertex(vertexId, distance, markers), *this);
    vertexMap.insert(make_pair(vertexId, v));
    return v;
}



// Write the graph in Graphviz format.
void LocalMarkerGraph2::write(
    const string& fileName,
    size_t minCoverage,
    size_t minConsensus,
    int maxDistance,
    bool detailed) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, minCoverage, minConsensus, maxDistance, detailed);
}
void LocalMarkerGraph2::write(
    ostream& s,
    size_t minCoverage,
    size_t minConsensus,
    int maxDistance,
    bool detailed) const
{
    Writer writer(*this, minCoverage, minConsensus, maxDistance, detailed);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&LocalMarkerGraph2Vertex::vertexId, *this));
}

LocalMarkerGraph2::Writer::Writer(
    const LocalMarkerGraph2& graph,
    size_t minCoverage,
    size_t minConsensus,
    int maxDistance,
    bool detailed) :
    graph(graph),
    minCoverage(minCoverage),
    minConsensus(minConsensus),
    maxDistance(maxDistance),
    detailed(detailed)
{
}



void LocalMarkerGraph2::Writer::operator()(std::ostream& s) const
{
    if(detailed) {
        s << "layout=dot;\n";
        s << "ratio=expand;\n";
        s << "node [fontname = \"Courier New\" shape=rectangle];\n";
        s << "edge [fontname = \"Courier New\"];\n";
    } else {
        s << "layout=sfdp;\n";
        s << "ratio=expand;\n";
        s << "node [shape=point];\n";
    }
}


void LocalMarkerGraph2::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const LocalMarkerGraph2Vertex& vertex = graph[v];
    const auto coverage = vertex.markers.size();


    // For compact output, the node shape is already defaulted to point,
    // and we don't write a label. The tooltip contains the vertex id,
    // which can be used to create a local subgraph to be looked at
    // in detailed format (use scripts/CreateLocalSubgraph.py).
    if(!detailed) {

        // Compact output.

        // Begin vertex attributes.
        s << "[";

        // Tooltip.
        s << "tooltip=\"Marker " << vertex.vertexId << ", coverage " << vertex.markers.size() << ", distance " << vertex.distance << "\"";

        // Vertex size.
        s << " width=\"";
        const auto oldPrecision = s.precision(4);
        s << 0.05 * sqrt(double(coverage));
        s.precision(oldPrecision);
        s << "\"";

        // Color.
        string color;
        if(vertex.distance == maxDistance) {
            color = "cyan";
        } else if(vertex.distance == 0) {
            color = "lightGreen";
        } else  if(coverage >= minCoverage) {
            color = "black";
        } else if(coverage == 1) {
            color = "#ff000080";  // Red, half way transparent
        } else if(coverage == 2) {
            color = "#ff800080";  // Orange, half way transparent
        } else {
            color = "#ff80ff80";  // Purple, half way transparent
        }
        s << " fillcolor=\"" << color << "\" color=\"" << color << "\"";

        // End vertex attributes.
        s << "]";

    } else {

        // Detailed output.

        // Begin vertex attributes.
        s << "[";

        // Color.
        string color;
        if(vertex.distance == maxDistance) {
            color = "cyan";
        } else if(vertex.distance == 0) {
            color = "lightGreen";
        } else if(coverage >= minCoverage) {
            color = "green";
        } else if(coverage == 1) {
            color = "#ff0000";  // Red
        } else if(coverage == 2) {
            color = "#ff8000";  // Orange
        } else {
            color = "#ff80ff";  // Purple
        }
        s << " style=filled";
        s << " fillcolor=\"" << color << "\"";

        // Label.
        s << "label=\"Marker " << vertex.vertexId;
        s << "\\nCoverage " << coverage;
        s << "\\nDistance " << vertex.distance << "\"";

        // End vertex attributes.
        s << "]";
    }
}



void LocalMarkerGraph2::Writer::operator()(std::ostream& s, edge_descriptor e) const
{

#if 0
    const LocalMarkerGraph2Edge& edge = graph[e];

    if(!detailed) {

        // Compact output.

        // Begin edge attributes.
        s << "[";

        // End edge attributes.
        s << "]";
    } else {

        // Detailed output.

        // If getting here, we are doing detailed output.

        // Begin edge attributes.
        s << "[";

        // End edge attributes.
        s << "]";
    }
#endif

}

