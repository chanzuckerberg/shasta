// Shasta.
#include "LocalReadGraph.hpp"
#include "Alignment.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard library.
#include "fstream.hpp"
#include "stdexcept.hpp"
#include "utility.hpp"



void LocalReadGraph::addVertex(
    OrientedReadId orientedReadId,
    uint32_t baseCount,
    bool isChimeric,
    uint32_t distance)
{
    // Check that we don't already have a vertex with this OrientedReadId.
    SHASTA_ASSERT(vertexMap.find(orientedReadId) == vertexMap.end());

    // Create the vertex.
    const vertex_descriptor v = add_vertex(LocalReadGraphVertex(
        orientedReadId, baseCount, isChimeric, distance), *this);

    // Store it in the vertex map.
    vertexMap.insert(make_pair(orientedReadId, v));
}



void LocalReadGraph::addEdge(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    uint32_t markerCount,
    AlignmentType alignmentType,
    bool crossesStrands)
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
        LocalReadGraphEdge(markerCount, alignmentType, crossesStrands),
        *this);
}



uint32_t LocalReadGraph::getDistance(OrientedReadId orientedReadId) const
{
    const auto it = vertexMap.find(orientedReadId);
    SHASTA_ASSERT(it != vertexMap.end());
    const vertex_descriptor v = it->second;
    return (*this)[v].distance;
}



bool LocalReadGraph::vertexExists(OrientedReadId orientedReadId) const
{
   return vertexMap.find(orientedReadId) != vertexMap.end();
}



// Write the graph in Graphviz format.
void LocalReadGraph::write(
        const string& fileName,
        const string& layoutMethod,
        uint32_t maxDistance,
        double vertexScalingFactor,
        double edgeThicknessScalingFactor,
        double edgeArrowScalingFactor,
        uint32_t maxTrim
) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, layoutMethod, maxDistance, vertexScalingFactor,
          edgeThicknessScalingFactor, edgeArrowScalingFactor, maxTrim);
}

void LocalReadGraph::write(
        ostream& s,
        const string& layoutMethod,
        uint32_t maxDistance,
        double vertexScalingFactor,
        double edgeThicknessScalingFactor,
        double edgeArrowScalingFactor,
        uint32_t maxTrim) const
{
    Writer writer(*this, layoutMethod, maxDistance, vertexScalingFactor,
                  edgeThicknessScalingFactor, edgeArrowScalingFactor,
                  maxTrim);

    boost::write_graphviz(s, *this, writer, writer, writer,
                          boost::get(&LocalReadGraphVertex::orientedReadIdValue, *this));
}


LocalReadGraph::Writer::Writer(
        const LocalReadGraph& graph,
        const string& layoutMethod,
        uint32_t maxDistance,
        double vertexScalingFactor,
        double edgeThicknessScalingFactor,
        double edgeArrowScalingFactor,
        uint32_t maxTrim) :
        graph(graph),
        layoutMethod(std::move(layoutMethod)),
        maxDistance(maxDistance),
        vertexScalingFactor(vertexScalingFactor),
        edgeThicknessScalingFactor(edgeThicknessScalingFactor),
        edgeArrowScalingFactor(edgeArrowScalingFactor),
        maxTrim(maxTrim)
{
}


void LocalReadGraph::Writer::operator()(std::ostream& s) const
{
    s << "layout=" + layoutMethod + ";\n";
    s << "ratio=expand;\n";
    s << "node [shape=point];\n";
    s << "edge [penwidth=\"0.2\"];\n";

    // This turns off the tooltip on the graph.
    s << "tooltip = \" \";\n";
}


void LocalReadGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const LocalReadGraphVertex& vertex = graph[v];
    const OrientedReadId orientedReadId(vertex.orientedReadId);

    s <<
        "["
        " tooltip=\"Read " << orientedReadId << ", " << vertex.markerCount <<
        " markers, distance " << vertex.distance << vertex.additionalToolTipText << "\"" <<
        " URL=\"exploreRead?readId=" << orientedReadId.getReadId() <<
        "&strand=" << orientedReadId.getStrand() <<
        "\"" <<
        " width=" << vertexScalingFactor * sqrt(1.e-6 * double(vertex.markerCount)) <<
        " height=" << vertexScalingFactor * sqrt(1.e-6 * double(vertex.markerCount)) <<

        // Id, so we can manipulate the vertex in javascript.
        " id=\"Vertex-" << orientedReadId << "\"";

    if(vertex.distance == 0) {
        s << " color=green fillcolor=green";
    } else if(vertex.distance == maxDistance) {
            s << " color=cyan fillcolor=cyan";
    } else if(vertex.isChimeric) {
        s << " color=red fillcolor=red";
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



    s <<
        "["
        "tooltip=\"" << vertex0.orientedReadId << " " <<
        vertex1.orientedReadId <<
        ", " << edge.markerCount << " aligned markers\"";


    // A containment alignment is drawn in red, at default thickness.
    // A non-containment alignment is drawn in black,
    // with thickness determined by the number of aligned markers.
    if( edge.alignmentType == AlignmentType::read0IsContained ||
        edge.alignmentType == AlignmentType::read1IsContained) {
        s << " color=red";
    } else {
        s << " penwidth=\"" << edgeThicknessScalingFactor * (1.e-4 * edge.markerCount) << "\"";
        s << " arrowsize=\"" << edgeArrowScalingFactor * 0.3 << "\"";
    }


    // An edge that crosses strands is drawn dashed.
    if(edge.crossesStrands) {
        s << " style=dashed";
    }


#if 1
    // The AlignmentType determines the edge endings.
    // Note that the Graphviz convention for undirected graphs
    // is that the head is the second vertex and the tail is the first vertex.
    s << " dir=both ";
    switch(edge.alignmentType) {
    case AlignmentType::read0IsContained:
        s << "arrowhead=none arrowtail=none";
        break;
    case AlignmentType::read1IsContained:
        s << "arrowhead=none arrowtail=none";
        break;
    case AlignmentType::read0IsBackward:
        s << "arrowhead=normal arrowtail=none";
        break;
    case AlignmentType::read1IsBackward:
        s << "arrowhead=none arrowtail=normal";
        break;
    case AlignmentType::ambiguous:
    default:
        s << "arrowhead=diamond arrowtail=diamond";
    }
#endif

    s << "]";
}

