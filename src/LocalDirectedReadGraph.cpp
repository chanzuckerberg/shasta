// Shasta.
#include "LocalDirectedReadGraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard library.
#include "fstream.hpp"
#include "stdexcept.hpp"
#include "utility.hpp"



void LocalDirectedReadGraph::addVertex(
    OrientedReadId orientedReadId,
    uint64_t baseCount,
    uint64_t distance)
{
    // Check that we don't already have a vertex with this OrientedReadId.
    SHASTA_ASSERT(vertexMap.find(orientedReadId) == vertexMap.end());

    // Create the vertex.
    const vertex_descriptor v = add_vertex(LocalDirectedReadGraphVertex(
        orientedReadId, baseCount, distance), *this);

    // Store it in the vertex map.
    vertexMap.insert(make_pair(orientedReadId, v));
}



void LocalDirectedReadGraph::addEdge(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int twiceOffsetAtCenter,
    uint64_t markerCount)
{
    // Find the vertices corresponding to these two OrientedReadId.
    const auto it0 = vertexMap.find(orientedReadId0);
    SHASTA_ASSERT(it0 != vertexMap.end());
    const vertex_descriptor v0 = it0->second;
    const auto it1 = vertexMap.find(orientedReadId1);
    SHASTA_ASSERT(it1 != vertexMap.end());
    const vertex_descriptor v1 = it1->second;

    // Add the edge.
    add_edge(v0, v1,
        LocalDirectedReadGraphEdge(twiceOffsetAtCenter, markerCount),
        *this);
}



uint64_t LocalDirectedReadGraph::getDistance(OrientedReadId orientedReadId) const
{
    const auto it = vertexMap.find(orientedReadId);
    SHASTA_ASSERT(it != vertexMap.end());
    const vertex_descriptor v = it->second;
    return (*this)[v].distance;
}



bool LocalDirectedReadGraph::vertexExists(OrientedReadId orientedReadId) const
{
   return vertexMap.find(orientedReadId) != vertexMap.end();
}



// Write the graph in Graphviz format.
void LocalDirectedReadGraph::write(const string& fileName, uint64_t maxDistance) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, maxDistance);
}
void LocalDirectedReadGraph::write(ostream& s, uint64_t maxDistance) const
{
    Writer writer(*this, maxDistance);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&LocalDirectedReadGraphVertex::orientedReadIdValue, *this));
}

LocalDirectedReadGraph::Writer::Writer(
    const LocalDirectedReadGraph& graph,
    uint64_t maxDistance) :
    graph(graph),
    maxDistance(maxDistance)
{
}



void LocalDirectedReadGraph::Writer::operator()(std::ostream& s) const
{
    s << "layout=sfdp;\n";
    s << "ratio=expand;\n";
    s << "node [shape=point];\n";
    s << "edge [penwidth=\"0.2\"];\n";

    // This turns off the tooltip on the graph.
    s << "tooltip = \" \";\n";
}


void LocalDirectedReadGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const LocalDirectedReadGraphVertex& vertex = graph[v];
    const OrientedReadId orientedReadId(vertex.orientedReadId);

    s <<
        "["
        " tooltip=\"Read " << orientedReadId << ", " << vertex.markerCount <<
        " markers, distance " << vertex.distance << vertex.additionalToolTipText << "\"" <<
        " URL=\"exploreRead?readId=" << orientedReadId.getReadId() <<
        "&strand=" << orientedReadId.getStrand() <<
        "\"" <<
        " width=" << sqrt(1.e-5 * double(vertex.markerCount));
    if(vertex.distance == 0) {
        s << " color=green fillcolor=green";
    } else if(vertex.distance == maxDistance) {
            s << " color=cyan fillcolor=cyan";
    }
    s << "]";
}



void LocalDirectedReadGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
    // const LocalDirectedReadGraphEdge& edge = graph[e];
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);
    const LocalDirectedReadGraphVertex& vertex0 = graph[v0];
    const LocalDirectedReadGraphVertex& vertex1 = graph[v1];

    s << "[";

    s <<
        "tooltip=\"" << vertex0.orientedReadId << " " <<
        vertex1.orientedReadId << "\"";

    s << "]";
}

