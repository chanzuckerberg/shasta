// Shasta.
#include "LocalAssemblyGraph.hpp"
#include "approximateTopologicalSort.hpp"
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
    MarkerGraph::VertexId markerGraphVertexId,
    int distance)
{
    // Check that the vertex does not already exist.
    SHASTA_ASSERT(vertexMap.find(assemblyGraphVertexId) == vertexMap.end());

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


    // Get the global edge id for this edge or its reverse complemented
    // (the one that was assembled).
    const LocalAssemblyGraph& graph = *this;
    const LocalAssemblyGraphEdge& edge = graph[e];
    EdgeId edgeId = edge.edgeId;
    if(!globalAssemblyGraph.isAssembledEdge(edgeId)) {
        edgeId = globalAssemblyGraph.reverseComplementEdge[edgeId];
    }
    SHASTA_ASSERT(globalAssemblyGraph.isAssembledEdge(edgeId));

    // Get the repeat counts for this edge.
    const MemoryAsContainer<uint8_t> repeatCounts = globalAssemblyGraph.repeatCounts[edgeId];


    // Add them up.
    int count = 0;
    for(const uint8_t repeatCount: repeatCounts) {
        count += repeatCount;
    }
    return count;
}



// Approximate topological sort.
void LocalAssemblyGraph::approximateTopologicalSort()
{
    LocalAssemblyGraph& graph = *this;

    // Create a table containing, for each edge, the number
    // of corresponding marker graph edges.
    vector<pair<uint64_t, edge_descriptor> > edgeTable;
    BGL_FORALL_EDGES(e, graph, LocalAssemblyGraph) {
        const EdgeId edgeId = graph[e].edgeId;
        const size_t markerGraphEdgeCount = globalAssemblyGraph.edgeLists.size(edgeId);
        edgeTable.push_back(make_pair(markerGraphEdgeCount, e));
    }

    // Sort by decreasing corresponding marker graph edges.
    sort(edgeTable.begin(), edgeTable.end(),
        std::greater< pair<uint32_t, edge_descriptor> >());

    // Gather the edge descriptors in this order.
    vector<edge_descriptor> sortedEdges;
    for(const auto& p: edgeTable) {
        sortedEdges.push_back(p.second);
    }

    // Do the approximate topological sort using this order.
    shasta::approximateTopologicalSort(graph, sortedEdges);

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
        "&maxDistance=10"
        "&timeout=30"
        "\"";

    // Label.
    if(showVertexLabels) {
        s << " label=\"" <<
            vertex.assemblyGraphVertexId <<
            "\\n" << vertex.markerGraphVertexId <<
            "\" style=filled fillcolor=\"#999999\"";
    }

    // End vertex attributes.
    s << "]";
}



void LocalAssemblyGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
    // Get the information we need.
    const LocalAssemblyGraphEdge& edge = graph[e];
    const EdgeId edgeId = edge.edgeId;
    const EdgeId edgeIdRc = graph.globalAssemblyGraph.reverseComplementEdge[edgeId];
    const size_t length = graph.edgeLength(e);
    const int baseCount = graph.baseCount(e);
    // const bool wasAssembled = graph.globalAssemblyGraph.isAssembledEdge(edgeId);
    const AssemblyGraph::Edge& globalEdge = graph.globalAssemblyGraph.edges[edgeId];

    // Begin edge attributes.
    s << "[";

    // Tooltip.
    const string tooltipText =
        "Assembly graph edge " + to_string(edgeId) + ", " +
        " reverse complement of assembly graph edge " + to_string(edgeIdRc) + ", " +
        to_string(length) + " marker graph edges, " +
        to_string(baseCount) + " bases" +
        ", vertex coverage min/ave/max " +
        to_string(globalEdge.minVertexCoverage) + "/" +
        to_string(globalEdge.averageVertexCoverage) + "/" +
        to_string(globalEdge.maxVertexCoverage) +
        ", edge coverage min/ave/max " +
        to_string(globalEdge.minEdgeCoverage) + "/" +
        to_string(globalEdge.averageEdgeCoverage) + "/" +
        to_string(globalEdge.maxEdgeCoverage);

    s << " tooltip=\"" << tooltipText << "\"";
    s << " labeltooltip=\"" << tooltipText << "\"";

    // URL
    s << " URL=\"exploreAssemblyGraphEdge?edgeId=" << edge.edgeId << "\"";

    // Color. This controls the color of the edge arrow.
    // const string color = wasAssembled ? "green" : "red";
    // s << " color=\"" << color << "\"";

    // Label.
    if(showEdgeLabels) {
        const string labelColor = "pink";
        s <<
            " label=<<table"
            " color=\"black\""
            " bgcolor=\"" << labelColor << "\""
            " border=\"1\""
            " cellborder=\"0\""
            " cellspacing=\"0\""
            ">"
            "<tr><td>" << edgeId << "</td></tr>"
            "<tr><td>" << length << "</td></tr>"
            "<tr><td>" << baseCount << "</td></tr>"
            "<tr><td>" <<
            globalEdge.minVertexCoverage << "/" <<
            globalEdge.averageVertexCoverage << "/" <<
            globalEdge.maxVertexCoverage <<
            "</td></tr>"
            "<tr><td>" <<
            globalEdge.minEdgeCoverage << "/" <<
            globalEdge.averageEdgeCoverage << "/" <<
            globalEdge.maxEdgeCoverage <<
            "</td></tr>"
            "</table>> decorate=true";
    }

    // If the edge was not marked as a DAG edge during approximate topological sort,
    // tell graphviz not to use it in constraint assignment.
    // This results in better graph layouts when using dot,
    // because back-edges tend to be low coverage edges.
    if(useDotLayout && !edge.isDagEdge) {
        s << " constraint=false";
    }

    // Thickness.
    const double thickness = 0.05 * std::pow(double(baseCount), 0.333333);
    s << " penwidth=";
    const auto oldPrecision = s.precision(4);
    s <<  thickness;
    s.precision(oldPrecision);


    // End edge attributes.
    s << "]";
}
