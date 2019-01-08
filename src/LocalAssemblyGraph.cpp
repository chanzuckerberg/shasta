// Shasta.
#include "LocalAssemblyGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard libraries.
#include "fstream.hpp"



LocalAssemblyGraph::LocalAssemblyGraph(AssemblyGraph& globalAssemblyGraph) :
    globalAssemblyGraph(globalAssemblyGraph)
{

}



// Add a vertex with the given vertex ids
// and return its vertex descriptor.
// A vertex with this VertexId must not exist.
LocalAssemblyGraph::vertex_descriptor LocalAssemblyGraph::addVertex(
    AssemblyGraph::VertexId assemblyGraphVertexId,
    GlobalMarkerGraphVertexId markerGraphVertexId,
    int distance)
{
    // Check that the vertex does not already exist.
    CZI_ASSERT(vertexMap.find(assemblyGraphVertexId) == vertexMap.end());

    // Add the vertex and store it in the vertex map.
    const vertex_descriptor v = add_vertex(LocalAssemblyGraphVertex(
        assemblyGraphVertexId, markerGraphVertexId, distance), *this);
    vertexMap.insert(make_pair(assemblyGraphVertexId, v));

    return v;
}



// Find out if a vertex with the given assembly graph vertex id exists.
// If it exists, return make_pair(true, v).
// Otherwise, return make_pair(false, null_vertex());
std::pair<bool, LocalAssemblyGraph::vertex_descriptor>
    LocalAssemblyGraph::findVertex(AssemblyGraph::VertexId assemblyGraphVertexId) const
{
    const auto it = vertexMap.find(assemblyGraphVertexId);
    if(it == vertexMap.end()) {
        return make_pair(false, null_vertex());
    } else {
        const vertex_descriptor v = it->second;
        return make_pair(true, v);
    }
}



// Return the number of marker graph edges that an edge corresponds to.
size_t LocalAssemblyGraph::edgeLength(edge_descriptor e) const
{
    const LocalAssemblyGraph& graph = *this;
    const EdgeId edgeId = graph[e].edgeId;
    return globalAssemblyGraph.edgeLists.size(edgeId);
}


// Return the number of bases in the raw assembled sequence of an edge,
// or -1 if not available.
int LocalAssemblyGraph::baseCount(edge_descriptor e) const
{
    if(!globalAssemblyGraph.repeatCounts.isOpen()) {
        return -1;
    }

    // Get the repeat counts for this vertex.
    const LocalAssemblyGraph& graph = *this;
    const EdgeId edgeId = graph[e].edgeId;
    const MemoryAsContainer<uint8_t> repeatCounts = globalAssemblyGraph.repeatCounts[edgeId];


    // Add them up.
    int count = 0;
    for(const uint8_t repeatCount: repeatCounts) {
        count += repeatCount;
    }
    return count;
}




// Write the graph in Graphviz format.
void LocalAssemblyGraph::write(
    const string& fileName,
    int maxDistance,
    bool useDotLayout,
    bool showVertexLabels,
    bool showEdgeLabels)
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, maxDistance, useDotLayout, showVertexLabels, showEdgeLabels);
}
void LocalAssemblyGraph::write(
    ostream& s,
    int maxDistance,
    bool useDotLayout,
    bool showVertexLabels,
    bool showEdgeLabels)
{
    Writer writer(*this, maxDistance, useDotLayout, showVertexLabels, showEdgeLabels);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&LocalAssemblyGraphVertex::assemblyGraphVertexId, *this));
}



LocalAssemblyGraph::Writer::Writer(
    LocalAssemblyGraph& graph,
    int maxDistance,
    bool useDotLayout,
    bool showVertexLabels,
    bool showEdgeLabels) :
    graph(graph),
    maxDistance(maxDistance),
    useDotLayout(useDotLayout),
    showVertexLabels(showVertexLabels),
    showEdgeLabels(showEdgeLabels)
{
}



void LocalAssemblyGraph::Writer::operator()(std::ostream& s) const
{
    // Turn off the tooltip on the graph.
    s << "tooltip = \" \";\n";

    // Define how to use extra space.
    s << "ratio=expand;\n";

    // Graph layout.
    if(useDotLayout) {
        s << "layout=dot;\n";
        s << "rankdir=LR;\n";
    } else {
        s << "layout=sfdp;\n";
        s << "smoothing=triangle;\n";
    }

    // Vertex shape.
    if(showVertexLabels) {
        s << "node [shape=rectangle];\n";
    } else {
        s << "node [shape=point];\n";
    }

}



// Write a vertex in graphviz format.
void LocalAssemblyGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const LocalAssemblyGraphVertex& vertex = graph[v];

    // Begin vertex attributes.
    s << "[";

    // Tooltip.
    s << "tooltip=\""
        "Assembly graph vertex " << vertex.assemblyGraphVertexId <<
        ", marker graph vertex " << vertex.markerGraphVertexId <<
        "\"";

    // URL.
    s << " URL=\"exploreMarkerGraph?"
        "?vertexId=" << vertex.markerGraphVertexId <<
        "&useStoredConnectivity=on"
        "&maxDistance=10"
        "&minCoverage=0"
        "&timeout=30"
        "&useBubbleReplacementEdges=on"
        "&showVertexId=on"
        "\"";

    // Label.
    if(showVertexLabels) {
        s << " label=\""
            "AGVID " << vertex.assemblyGraphVertexId <<
            "\\nMGVID " << vertex.markerGraphVertexId <<
            "\"";
    }

    // End vertex attributes.
    s << "]";


    #if 0
    const LocalAssemblyGraphVertex& vertex = graph[v];
    const size_t length = graph.vertexLength(v);
    const int baseCount = graph.baseCount(v);

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
        s << ", length " << length;
        if(baseCount >= 0) {
            s << ", " << baseCount << " bases";
        }
        s << "\"";
    } else {
        s << " tooltip=\" \"";
    }



    // Label.
    if(detailed) {
        s << " label=\"Vertex " << vertex.vertexId << "\\n";
        s << length << " edges";
        if(baseCount >= 0) {
            s << "\\n" << baseCount << " bases";
        }
        s << "\"";
    }


    // Link to detailed information for this vertex.
    s << " URL=\"exploreAssemblyGraphVertex?vertexId=" << vertex.vertexId << "\"";


    // End vertex attributes.
    s << "]";
#endif
}



void LocalAssemblyGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
    const LocalAssemblyGraphEdge& edge = graph[e];
    const size_t length = graph.edgeLength(e);
    const int baseCount = graph.baseCount(e);

    // Begin edge attributes.
    s << "[";

    // Tooltip.
    s << "tooltip=\"Edge " << edge.edgeId << ", " << baseCount << " bases, " << length << " marker graph edges\"";
    s << "labeltooltip=\"Edge " << edge.edgeId << ", " << baseCount << " bases, " << length << " marker graph edges\"";

    // URL
    s << " URL=\"exploreAssemblyGraphEdge?edgeId=" << edge.edgeId << "\"";

    // Label.
    if(showEdgeLabels) {
        s << "label=\"" << baseCount << "\"";
    }


    // End edge attributes.
    s << "]";
}
