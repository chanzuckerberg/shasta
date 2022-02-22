// Shasta.
#include "LocalReadGraph.hpp"
#include "Alignment.hpp"
#include "writeGraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard library.
#include "fstream.hpp"
#include "stdexcept.hpp"



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
    uint64_t globalEdgeId,
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
        LocalReadGraphEdge(markerCount, crossesStrands, globalEdgeId),
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
        double edgeThicknessScalingFactor
) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, layoutMethod, maxDistance, vertexScalingFactor,
          edgeThicknessScalingFactor);
}

void LocalReadGraph::write(
        ostream& s,
        const string& layoutMethod,
        uint32_t maxDistance,
        double vertexScalingFactor,
        double edgeThicknessScalingFactor) const
{
    Writer writer(*this, layoutMethod, maxDistance, vertexScalingFactor,
                  edgeThicknessScalingFactor);

    boost::write_graphviz(s, *this, writer, writer, writer,
                          boost::get(&LocalReadGraphVertex::orientedReadIdValue, *this));
}


LocalReadGraph::Writer::Writer(
        const LocalReadGraph& graph,
        const string& layoutMethod,
        uint32_t maxDistance,
        double vertexScalingFactor,
        double edgeThicknessScalingFactor) :
        graph(graph),
        layoutMethod(layoutMethod),
        maxDistance(maxDistance),
        vertexScalingFactor(vertexScalingFactor),
        edgeThicknessScalingFactor(edgeThicknessScalingFactor)
{
}


void LocalReadGraph::Writer::operator()(std::ostream& s) const
{
    s << "layout=" + layoutMethod + ";\n";
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

    // Edge thickness is determined by the number of aligned markers.
    s << " penwidth=\"" << edgeThicknessScalingFactor * (1.e-4 * edge.markerCount) << "\"";

    // An edge that crosses strands is drawn dashed.
    if(edge.crossesStrands) {
        s << " style=dashed";
    }

    s << "]";
}



// Compute sfdp layout using graphviz and store the results
// in the vertex positions.
ComputeLayoutReturnCode LocalReadGraph::computeLayout(
    const string& layoutMethod,
    double timeout)
{
    LocalReadGraph& graph = *this;

    // Compute the layout.
    std::map<vertex_descriptor, array<double, 2> > positionMap;
    const ComputeLayoutReturnCode returnCode =
        shasta::computeLayoutGraphviz(graph, layoutMethod, timeout, positionMap);
    if(returnCode != ComputeLayoutReturnCode::Success) {
        return returnCode;
    }

    // Store it in the vertices.
    BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
        const auto it = positionMap.find(v);
        SHASTA_ASSERT(it != positionMap.end());
        graph[v].position = it->second;
    }
    return ComputeLayoutReturnCode::Success;
}



// Write directly to svg, without using Graphviz rendering.
// This assumes that the layout was already computed
// and stored in the vertices.
void LocalReadGraph::writeSvg(
    const string& svgId,
    uint64_t width,
    uint64_t height,
    double vertexScalingFactor,
    double edgeThicknessScalingFactor,
    uint64_t maxDistance,
    ostream& svg) const
{
    using Graph = LocalReadGraph;
    using VertexAttributes = WriteGraph::VertexAttributes;
    using EdgeAttributes = WriteGraph::EdgeAttributes;
    const Graph& graph = *this;



    // Fill in vertex attributes.
    std::map<vertex_descriptor, VertexAttributes> vertexAttributes;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const auto& vertex = graph[v];
        const OrientedReadId orientedReadId = vertex.orientedReadId;
        VertexAttributes attributes;

        attributes.radius = vertexScalingFactor * sqrt(3.e-7 * double(vertex.markerCount));

        attributes.id = "Vertex-" + orientedReadId.getString();

        if(vertex.distance == 0) {
            attributes.color = "lime";
        } else if(vertex.distance == maxDistance) {
            attributes.color = "cyan";
        } else if(vertex.isChimeric) {
            attributes.color = "red";
        }

        attributes.tooltip = "Read " + orientedReadId.getString() + ", " + to_string(vertex.markerCount) +
            " markers, distance " + to_string(vertex.distance) + vertex.additionalToolTipText;

        attributes.url = "exploreRead?readId=" + to_string(orientedReadId.getReadId()) +
            "&strand=" + to_string(orientedReadId.getStrand());

        vertexAttributes.insert(make_pair(v, attributes));
    }



    // Fill in edge attributes.
    std::map<edge_descriptor, EdgeAttributes> edgeAttributes;
    BGL_FORALL_EDGES(e, graph, Graph) {
        const auto& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const auto& vertex0 = graph[v0];
        const auto& vertex1 = graph[v1];

        EdgeAttributes attributes;

        attributes.thickness = edgeThicknessScalingFactor * 1.e-6 * double(edge.markerCount);
        if(edge.color.empty()) {
            attributes.color = "midnightblue";
        } else {
            attributes.color = edge.color;
        }

        attributes.tooltip = vertex0.orientedReadId.getString() + " " +
            vertex1.orientedReadId.getString() +
            ", " + to_string(edge.markerCount) + " aligned markers";

        edgeAttributes.insert(make_pair(e, attributes));
    }

    // Write to svg.
    WriteGraph::writeSvg(graph, svgId, width, height,
        vertexAttributes, edgeAttributes, svg);
}
