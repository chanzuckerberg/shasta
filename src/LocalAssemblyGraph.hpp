#ifndef CZI_SHASTA_LOCAL_ASSEMBLY_GRAPH_HPP
#define CZI_SHASTA_LOCAL_ASSEMBLY_GRAPH_HPP


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

namespace ChanZuckerberg {
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
}



class ChanZuckerberg::shasta::LocalAssemblyGraphVertex {
public:

    // The global vertex id of the vertex of the global assembly
    // graph that corresponds to this vertex.
    AssemblyGraph::VertexId vertexId;

    // The distance from the start vertex.
    int distance;

    LocalAssemblyGraphVertex(
        AssemblyGraph::VertexId vertexId,
        int distance) :
        vertexId(vertexId),
        distance(distance)
        {}

};



class ChanZuckerberg::shasta::LocalAssemblyGraphEdge {
public:
};



class ChanZuckerberg::shasta::LocalAssemblyGraph :
    public LocalAssemblyGraphBaseClass {
public:

    LocalAssemblyGraph(
        AssemblyGraph&
        );

    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::EdgeId;

    // Add a vertex with the given vertex id
    // and return its vertex descriptor.
    // A vertex with this vertex id must not exist.
    vertex_descriptor addVertex(
        VertexId,
        int distance);

    // Find out if a vertex with the given VertexId exists.
    // If it exists, return make_pair(true, v).
    // Otherwise, return make_pair(false, null_vertex());
    pair<bool, vertex_descriptor> findVertex(VertexId) const;

    // Return the number of marker graph edges that the vertex corresponds to.
    size_t vertexLength(vertex_descriptor) const;

    // Return the number of bases in the raw assembled sequence of a vertex,
    // or -1 if not available.
    int baseCount(vertex_descriptor) const;

    // Write in Graphviz format.
    void write(
        ostream&,
        int maxDistance,
        bool detailed);
    void write(
        const string& fileName,
        int maxDistance,
        bool detailed);


private:

    // Map a global vertex id to a vertex descriptor for the local graph.
    std::map<VertexId, vertex_descriptor> vertexMap;

    // Reference to the global assembly graph.
    AssemblyGraph& globalAssemblyGraph;

    // Writer class used for Graphviz output.
    class Writer {
    public:
        Writer(
            LocalAssemblyGraph&,
            int maxDistance,
            bool detailed);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const LocalAssemblyGraph& graph;
        int maxDistance;
        bool detailed;
    };
    friend class Writer;
};


#endif
