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
    uint32_t baseCount,
    uint32_t markerCount,
    uint32_t distance,
    bool isContained)
{
    // Check that we don't already have a vertex with this OrientedReadId.
    SHASTA_ASSERT(vertexMap.find(orientedReadId) == vertexMap.end());

    // Create the vertex.
    const vertex_descriptor v = add_vertex(LocalDirectedReadGraphVertex(
        orientedReadId, baseCount, markerCount, distance, isContained), *this);

    // Store it in the vertex map.
    vertexMap.insert(make_pair(orientedReadId, v));
}



void LocalDirectedReadGraph::addEdge(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const AlignmentInfo& alignmentInfo,
    bool involvesTwoContainedVertices,
    bool involvesOneContainedVertex,
    bool keep,
    bool isConflict,
    uint32_t commonNeighborCount)
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
        LocalDirectedReadGraphEdge(alignmentInfo,
            involvesTwoContainedVertices,
            involvesOneContainedVertex,
            keep,
            isConflict,
            commonNeighborCount),
        *this);
}



uint32_t LocalDirectedReadGraph::getDistance(OrientedReadId orientedReadId) const
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
void LocalDirectedReadGraph::write(
    const string& fileName,
    uint32_t maxDistance,
    double vertexScalingFactor,
    double edgeThicknessScalingFactor,
    double edgeArrowScalingFactor,
    bool colorEdgeArrows,
    bool dashedContainmentEdges,
    uint32_t maxTrim,
    VertexColoringMethod vertexColoringMethod
    ) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, maxDistance, vertexScalingFactor,
        edgeThicknessScalingFactor, edgeArrowScalingFactor,
        colorEdgeArrows,
        dashedContainmentEdges, maxTrim,
        vertexColoringMethod);
}
void LocalDirectedReadGraph::write(
    ostream& s,
    uint32_t maxDistance,
    double vertexScalingFactor,
    double edgeThicknessScalingFactor,
    double edgeArrowScalingFactor,
    bool colorEdgeArrows,
    bool dashedContainmentEdges,
    uint32_t maxTrim,
    VertexColoringMethod vertexColoringMethod) const
{
    Writer writer(*this, maxDistance, vertexScalingFactor,
        edgeThicknessScalingFactor, edgeArrowScalingFactor,
        colorEdgeArrows,
        dashedContainmentEdges,
        maxTrim,
        vertexColoringMethod);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&LocalDirectedReadGraphVertex::orientedReadIdValue, *this));
}

LocalDirectedReadGraph::Writer::Writer(
    const LocalDirectedReadGraph& graph,
    uint32_t maxDistance,
    double vertexScalingFactor,
    double edgeThicknessScalingFactor,
    double edgeArrowScalingFactor,
    bool colorEdgeArrows,
    bool dashedContainmentEdges,
    uint32_t maxTrim,
    VertexColoringMethod vertexColoringMethod) :
    graph(graph),
    maxDistance(maxDistance),
    vertexScalingFactor(vertexScalingFactor),
    edgeThicknessScalingFactor(edgeThicknessScalingFactor),
    edgeArrowScalingFactor(edgeArrowScalingFactor),
    colorEdgeArrows(colorEdgeArrows),
    dashedContainmentEdges(dashedContainmentEdges),
    maxTrim(maxTrim),
    vertexColoringMethod(vertexColoringMethod)
{
}



void LocalDirectedReadGraph::Writer::operator()(std::ostream& s) const
{
    s << "layout=sfdp;\n";
    s << "ratio=expand;\n";
    s << "smoothing=triangle;\n";
    s << "node [shape=point];\n";


    s << "edge [dir=both arrowtail=inv];\n";
    if(colorEdgeArrows) {
        s << "edge [color=\"green:black;0.9:red\"];\n";
    }

    // This turns off the tooltip on the graph.
    s << "tooltip = \" \";\n";
}


void LocalDirectedReadGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const LocalDirectedReadGraphVertex& vertex = graph[v];
    const OrientedReadId orientedReadId(vertex.orientedReadId);

    const bool hasClusterInformation =
        vertex.color != std::numeric_limits<uint32_t>::max();

    // Tooltip.
    s <<
        "["
        " tooltip=\"Read " << orientedReadId << ", " <<
        vertex.baseCount << " bases, " << vertex.markerCount <<
        " markers, distance " << vertex.distance;
    if(vertex.conflictCount) {
        s << ", " << vertex.conflictCount << " conflicting vertices" ;
    }
    if(hasClusterInformation) {
        s << ", conflict read graph component " << vertex.componentId <<
            " color " << vertex.color;
     }
    s << vertex.additionalToolTipText << "\"" <<
        " URL=\"exploreRead?readId=" << orientedReadId.getReadId() <<
        "&strand=" << orientedReadId.getStrand() <<
        "\"" <<
        " width=" << vertexScalingFactor * sqrt(1.e-6 * double(vertex.markerCount)) <<
        " height=" << vertexScalingFactor * sqrt(1.e-6 * double(vertex.markerCount));

    // Id, so we can manipulate the vertex in javascript.
    s << " id=\"Vertex-" << orientedReadId << "\"";




    // Color.
    if(vertex.distance == maxDistance) {
        s << " color=cyan";
    } else if(vertex.isConflictingGreen) {
        s << " color=green";
    } else if(vertex.isConflictingRed) {
        s << " color=red";
    } else {

        switch(vertexColoringMethod) {

        case VertexColoringMethod::ByConflictCount:
            if(vertex.wasRemoved) {
                s << "color=orange";
            } else  if(vertex.conflictCount == 0) {
                s << "color=black";
            } else {
                const double hue = 0.67;
                const double saturation = 0.3;
                const double value = min(1.0, 0.5 + 0.05*double(vertex.conflictCount));
                s << " color=\"" << hue << ","<< saturation << "," << value << "\"";
            }
            break;

        case VertexColoringMethod::ByCluster:
            if(vertex.wasRemoved) {
                s << "color=orange";
            } else if(vertex.color != std::numeric_limits<uint32_t>::max()) {
                s << " color=\"/set18/" << (vertex.color % 8) + 1 << "\"";
            }
            break;

        default:
            break;
        }
    }



    // Shape.
    if(not vertex.additionalToolTipText.empty()) {
        s << " shape=diamond style=filled label=\"\"";
    }

    s << "]";
}



void LocalDirectedReadGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
    const LocalDirectedReadGraphEdge& edge = graph[e];
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);
    const LocalDirectedReadGraphVertex& vertex0 = graph[v0];
    const LocalDirectedReadGraphVertex& vertex1 = graph[v1];
    const OrientedReadId orientedReadId0 = vertex0.orientedReadId;
    const OrientedReadId orientedReadId1 = vertex1.orientedReadId;

    const bool isContainment = edge.alignmentInfo.isContaining(maxTrim);

    s << "[";

    s <<
        "tooltip=\"" << orientedReadId0 << "->" <<
        orientedReadId1 <<
        ", " << edge.alignmentInfo.markerCount << " aligned markers, centers offset " <<
        std::setprecision(6) << edge.alignmentInfo.offsetAtCenter() <<
        " aligned fraction " <<
        std::setprecision(3) <<
        edge.alignmentInfo.alignedFraction(0) << " " <<
        edge.alignmentInfo.alignedFraction(1) <<
        ", common neighbors " << edge.commonNeighborCount <<
        "\"";

    s << " penwidth=\"" << edgeThicknessScalingFactor * (1.e-4 * edge.alignmentInfo.markerCount) << "\"";
    s << " arrowsize=\"" << edgeArrowScalingFactor * 0.3 << "\"";

    if(dashedContainmentEdges and isContainment) {
        s << " style=dashed";
    }


    // Hyperlink to the alignment corresponding to this edge.
    s << " URL=\""
        "exploreAlignment"
        "?readId0=" << orientedReadId0.getReadId() <<
        "&strand0=" << orientedReadId0.getStrand() <<
        "&readId1=" << orientedReadId1.getReadId() <<
        "&strand1=" << orientedReadId0.getStrand() <<
        "\"";



    if(not edge.keep) {
        if(colorEdgeArrows) {
            s << " color=\"green:#0000ff7f;0.9:red\""; // Partially transparent blue with green/red ends
            s << " dir=both arrowtail=inv";
        } else {
            s << " color=\"#0000ff7f\""; // Partially transparent blue.
        }
    } else {

        if(edge.isConflict) {
            s << " color=\"#ff00007f\""; // Partially transparent red.
        }

    }



    s << "]";
}

