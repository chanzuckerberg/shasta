// Shasta.
#include "LocalReadGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard library.
#include "fstream.hpp"
#include "stdexcept.hpp"
#include "utility.hpp"



void LocalReadGraph::addVertex(
    ReadId readId,
    uint32_t baseCount,
    uint32_t distance)
{
    // Check that we don't already have a vertex with this ReadId.
    CZI_ASSERT(vertexMap.find(readId) == vertexMap.end());

    // Create the vertex.
    const vertex_descriptor v = add_vertex(LocalReadGraphVertex(readId, baseCount, distance), *this);

    // Store it in the vertex map.
    vertexMap.insert(make_pair(readId, v));
}



void LocalReadGraph::addEdge(
    ReadId readId0,
    ReadId readId1,
    size_t globalEdgeId,
    uint8_t direction0,
    uint8_t direction1)
{
    // Find the vertices corresponding to these two ReadId.
    const auto it0 = vertexMap.find(readId0);
    CZI_ASSERT(it0 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const auto it1 = vertexMap.find(readId1);
    CZI_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v1 = it1->second;

    // Add the edge.
    add_edge(v0, v1, LocalReadGraphEdge(globalEdgeId, direction0, direction1), *this);
}



uint32_t LocalReadGraph::getDistance(ReadId readId) const
{
    const auto it = vertexMap.find(readId);
    CZI_ASSERT(it != vertexMap.end());
    const vertex_descriptor v = it->second;
    return (*this)[v].distance;
}



bool LocalReadGraph::vertexExists(ReadId readId) const
{
   return vertexMap.find(readId) != vertexMap.end();
}



// Write the graph in Graphviz format.
void LocalReadGraph::write(const string& fileName, uint32_t maxDistance) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, maxDistance);
}
void LocalReadGraph::write(ostream& s, uint32_t maxDistance) const
{
    Writer writer(*this, maxDistance);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&LocalReadGraphVertex::readId, *this));
}

LocalReadGraph::Writer::Writer(
    const LocalReadGraph& graph,
    uint32_t maxDistance) :
    graph(graph),
    maxDistance(maxDistance)
{
}



void LocalReadGraph::Writer::operator()(std::ostream& s) const
{
    s << "layout=sfdp;\n";
    s << "ratio=expand;\n";
    s << "node [shape=point];\n";
    s << "edge [penwidth=\"0.2\"];\n";

    // This turns off the tooltip on the graph.
    s << "tooltip = \" \";\n";
}


void LocalReadGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const LocalReadGraphVertex& vertex = graph[v];
    const ReadId readId(vertex.readId);

    s <<
        "["
        " tooltip=\"Read " << readId << ", " << vertex.markerCount <<
        " markers, distance " << vertex.distance << vertex.additionalToolTipText << "\"" <<
        " URL=\"exploreRead?readId=" << readId <<
        "&strand=0\"" <<
        " width=" << sqrt(1.e-5 * vertex.markerCount);
    if(vertex.distance == 0) {
        s << " color=lightGreen fillcolor=lightGreen";
    } else if(vertex.distance == maxDistance) {
            s << " color=cyan fillcolor=cyan";
    }
    s << "]";
}



void LocalReadGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
    const LocalReadGraphEdge& edge = graph[e];
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);
    const LocalReadGraphVertex& vertex0 = graph[v0];
    const LocalReadGraphVertex& vertex1 = graph[v1];

    const string direction0 = edge.direction0 ? "normal" : "inv";
    const string direction1 = edge.direction1 ? "normal" : "inv";

    s <<
        "["
        "tooltip=\"" << vertex0.readId << " " <<
        vertex1.readId << "\""
        " dir=both arrowhead=" << direction0 <<
        " arrowtail=" << direction1 <<
        " arrowsize=0.4"
        "]";
}

