#ifndef SHASTA_ASSEMBLY_GRAPH2_HPP
#define SHASTA_ASSEMBLY_GRAPH2_HPP

// Assembly graph for assembly mode2.



// Shasta.
#include "MarkerGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>
#include "string.hpp"
#include "vector.hpp"



namespace shasta {
    class AssemblyGraph2;
    class AssemblyGraph2Vertex;
    class AssemblyGraph2Edge;
    class MarkerGraph;

    using AssemblyGraph2BaseClass =
        boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
        AssemblyGraph2Vertex, AssemblyGraph2Edge>;

    using MarkerGraphPath = vector<MarkerGraph::EdgeId>;
}



class shasta::AssemblyGraph2Vertex {
public:
    MarkerGraph::VertexId markerGraphVertexId;

    AssemblyGraph2Vertex(MarkerGraph::VertexId markerGraphVertexId) :
        markerGraphVertexId(markerGraphVertexId) {}
};



class shasta::AssemblyGraph2Edge {
public:

    // Each assembly graph edge corresponds to
    // a set of paths in the marker graph.
    // This way it can describe a bubble in the marker graph.
    vector<MarkerGraphPath> markerGraphPaths;

    // The default constructor creates an edge without any paths.
    AssemblyGraph2Edge() {}

    // This constructor creates an edge with a single path.
    AssemblyGraph2Edge(const MarkerGraphPath& path) :
        markerGraphPaths(1, path) {}

    bool isBubble() const
    {
        return markerGraphPaths.size() > 1;
    }
};



class shasta::AssemblyGraph2 : public AssemblyGraph2BaseClass{
public:

    // The constructor creates an edge for each linear path
    // in the marker graph. Therefore, immediately after construction,
    // each edge has a single MarkerGraphPath (no bubbles).
    AssemblyGraph2(const MarkerGraph&);

    void writeCsv(const string& baseName) const;
    void writeVerticesCsv(const string& baseName) const;
    void writeEdgesCsv(const string& baseName) const;
    void writeEdgeDetailsCsv(const string& baseName) const;

private:

    // A non-owned reference to the MarkerGraph.
    const MarkerGraph& markerGraph;

    // Map that gives us the vertex descriptor corresponding to
    // each marker graph vertex.
    std::map<MarkerGraph::VertexId, vertex_descriptor> vertexMap;

    // Get the vertex descriptor for the vertex corresponding to
    // a given MarkerGraph::VertexId, creating the vertex if necessary.
    vertex_descriptor getVertex(MarkerGraph::VertexId);

    // Create a new edge corresponding to the given path.
    // Also create the vertices if necessary.
    void addEdge(const MarkerGraphPath&);

    // Finds edges that form bubbles, then combine
    // each of them into a single edge with multiple paths.
    void gatherBubbles();

};



#endif

