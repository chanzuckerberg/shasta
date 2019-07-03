#ifndef SHASTA_LOCAL_ASSEMBLY_GRAPH_HPP
#define SHASTA_LOCAL_ASSEMBLY_GRAPH_HPP


/*******************************************************************************

The local marker graph created by class LocalMarkerGraph is a subgraph
of the global assembly graph, created by starting at a given vertex,
and extending out to a specified distance in both directions.
Distance is number of edges on the global assembly graph.

*******************************************************************************/

// Shasta
#include "AssemblyGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

namespace shasta {

    // Forward declaration of types declared in this file.
    class LocalAssemblyGraphVertex;
    class LocalAssemblyGraphEdge;
    class LocalAssemblyGraph;
    using LocalAssemblyGraphBaseClass = boost::adjacency_list<
        boost::listS,   // Allow parallel edges!
        boost::listS,
        boost::bidirectionalS,
        LocalAssemblyGraphVertex,
        LocalAssemblyGraphEdge
        >;

}



class shasta::LocalAssemblyGraphVertex {
public:

    // The vertex id of the vertex of the global assembly
    // graph that corresponds to this vertex.
    AssemblyGraph::VertexId assemblyGraphVertexId;

    // The vertex id of the vertex of the global marker
    // graph that corresponds to this vertex.
    MarkerGraph::VertexId markerGraphVertexId;

    // The distance from the start vertex.
    int distance;

    // Fields used by approximateTopologicalSort.
    uint32_t color = 0;
    size_t rank = 0;

    LocalAssemblyGraphVertex(
        AssemblyGraph::VertexId assemblyGraphVertexId,
        MarkerGraph::VertexId markerGraphVertexId,
        int distance) :
        assemblyGraphVertexId(assemblyGraphVertexId),
        markerGraphVertexId(markerGraphVertexId),
        distance(distance)
        {}

};



class shasta::LocalAssemblyGraphEdge {
public:
    // The global edge id of the edge of the global assembly
    // graph that corresponds to this edge.
    AssemblyGraph::EdgeId edgeId;

    // Field used by approximateTopologicalSort.
    bool isDagEdge = true;
};



class shasta::LocalAssemblyGraph :
    public LocalAssemblyGraphBaseClass {
public:

    LocalAssemblyGraph(
        AssemblyGraph&
        );

    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::EdgeId;

    // Add a vertex with the given vertex ids
    // and return its vertex descriptor.
    // A vertex with this vertex id must not exist.
    vertex_descriptor addVertex(
        VertexId,
        MarkerGraph::VertexId,
        int distance);

    // Find out if a vertex with the given assembly graph vertex id exists.
    // If it exists, return make_pair(true, v).
    // Otherwise, return make_pair(false, null_vertex());
    pair<bool, vertex_descriptor> findVertex(VertexId) const;

    // Return the number of marker graph edges that an edge corresponds to.
    size_t edgeLength(edge_descriptor) const;

    // Return the number of bases in the raw assembled sequence of an edge,
    // or -1 if not available.
    int baseCount(edge_descriptor) const;

    // Write in Graphviz format.
    void write(
        ostream&,
        int maxDistance,
        bool useDotLayout,
        bool showVertexLabels,
        bool showEdgeLabels);
    void write(
        const string& fileName,
        int maxDistance,
        bool useDotLayout,
        bool showVertexLabels,
        bool showEdgeLabels);

    // Approximate topological sort.
    void approximateTopologicalSort();

private:

    // Map a global assembly graph vertex id to a vertex descriptor for the local graph.
    std::map<VertexId, vertex_descriptor> vertexMap;

    // Reference to the global assembly graph.
    AssemblyGraph& globalAssemblyGraph;

    // Writer class used for Graphviz output.
    class Writer {
    public:
        Writer(
            LocalAssemblyGraph&,
            int maxDistance,
            bool useDotLayout,
            bool showVertexLabels,
            bool showEdgeLabels);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const LocalAssemblyGraph& graph;
        int maxDistance;
        bool useDotLayout;
        bool showVertexLabels;
        bool showEdgeLabels;
    };
    friend class Writer;
};


#endif
