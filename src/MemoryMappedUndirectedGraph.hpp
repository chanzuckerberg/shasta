#ifndef SHASTA_MEMORY_MAPPED_UNDIRECTED_GRAPH_HPP
#define SHASTA_MEMORY_MAPPED_UNDIRECTED_GRAPH_HPP

// Shasta
#include "deduplicate.hpp"
#include "MemoryMappedVector.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

// Standard library.
#include "chrono.hpp"
#include "iterator.hpp"
#include <limits>
#include <map>
#include <queue>
#include "string.hpp"

// An undirected graph stored on memory mapped data structures.
// It allows parallel edges (more than one edge
// between a pair of vertices), but not self-edges.
// It has the following restrictions:
// - Vertices and edges can only be added, not removed.
// - When an edge is added, the edgesBySource
//   and edgesByTarget data structures become invalid
//   and must be recreated with a call to computeConnectivity.

namespace shasta {
    namespace MemoryMapped {
        template<class Vertex, class Edge> class UndirectedGraph;
    }
}



template<class Vertex, class Edge> class shasta::MemoryMapped::UndirectedGraph {
public:

    // The types used to identify vertices and edges.
    using VertexId = uint64_t;
    using EdgeId = uint64_t;
    static const VertexId invalidVertexId = std::numeric_limits<VertexId>::max();
    static const EdgeId invalidEdgeId = std::numeric_limits<EdgeId>::max();

    // The graph vertices.
    // Vertices are identified by their index in this vector
    // (a VertexId).
    Vector<Vertex> vertices;
    VertexId addVertex(const Vertex& vertex)
    {
        const VertexId vertexId = vertices.size();
        vertices.push_back(vertex);
        return vertexId;
    }

    // Accessor for graph vertices.
    Vertex& getVertex(VertexId vertexId)
    {
        return vertices[vertexId];
    }
    const Vertex& getVertex(VertexId vertexId) const
    {
        return vertices[vertexId];
    }

    // The graph edges of the graph.
    // Vertices are identified by their index in this vector
    // (an EdgeId).
    // Edges are always stored with v0<v1.
    class EdgeInformation {
    public:
        VertexId v0;
        VertexId v1;
        Edge edge;
        EdgeInformation(VertexId v0, VertexId v1, const Edge& edge) :
            v0(min(v0, v1)), v1(max(v0, v1)), edge(edge)
        {
            SHASTA_ASSERT(v0 < v1);
        }
        EdgeInformation() : v0(invalidVertexId), v1(invalidVertexId) {}
    };
    Vector<EdgeInformation> edges;
    EdgeId addEdge(
        VertexId v0,
        VertexId v1,
        const Edge& edge)
    {
        const EdgeId edgeId = edges.size();
        edges.push_back(EdgeInformation(v0, v1, edge));
        return edgeId;
    }

    // Accessors for graph edges.
    Edge& getEdge(EdgeId edgeId) {
        return edges[edgeId].edge;
    }
    const Edge& getEdge(EdgeId edgeId) const {
        return edges[edgeId].edge;
    }

    // Return the lower-numbered vertex.
    VertexId v0(EdgeId edgeId) const {
        return edges[edgeId].v0;
    }

    // Return the higher-numbered vertex.
    VertexId target(EdgeId edgeId) const {
        return edges[edgeId].v1;
    }

    // The list of edges incident to a vertex.
    // Indexed by a VertexId.
    // That is, edgesByVertex[vertexId] contains the EdgeId's
    // of all edges incident to vertexId.
    // This vector becomes invalid whenever am edge is added
    // and it has to be recreated with a call to
    // computeConnectivity.
    VectorOfVectors<VertexId, uint64_t> edgesByVertex;

    // Accessors for edges by source and by target.
    MemoryAsContainer<EdgeId> incidentEdges(VertexId vertexId)
    {
        return edgesByVertex[vertexId];
    }
    uint64_t degree(VertexId vertexId)
    {
        return incidentEdges.size(vertexId);
    }


    // A call to this function recomputes from scratch
    // the edgesByVertex data structure.
    // This could be made multithreaded.
    void computeConnectivity()
    {
        // Clear out edgesByVertex data structure.
        edgesByVertex.clear();

        // Pass 1: count the number of edges that each vertex
        // is adjacent to.
        edgesByVertex.beginPass1(vertices.size());
        for(const EdgeInformation& edge: edges) {
            edgesByVertex.incrementCount(edge.v0);
            edgesByVertex.incrementCount(edge.v1);
        }

        // Pass 2: store the edges that each vertex
        // is the source/target of.
        edgesByVertex.beginPass2();
        for(EdgeId edgeId=0; edgeId<edges.size(); edgeId++) {
            const EdgeInformation& edge = edges[edgeId];
            edgesByVertex.store(edge.v0, edgeId);
            edgesByVertex.store(edge.v1, edgeId);
        }
        edgesByVertex.endPass2();

        // Pass 3: sort the egde ids for each vertex.
        for(VertexId vertexId=0; vertexId<vertices.size(); vertexId++) {
            MemoryAsContainer<EdgeId> s = incidentEdges(vertexId);
            sort(s.begin(), s.end());
        }

    }






    // The default constructor puts the graph in an uninitialized state.
    // The graph can only be used after calling createNew or accessExisting.
    UndirectedGraph() {}

    // Create a new UndirectedGraph.
    // If the name is empty, this is created in anonymous memory.
    // Otherwise, it is created in memory mapped files with the
    // specified base name.
    void createNew(
        const string& baseNameArgument,
        size_t pageSize)
    {
        baseName = baseNameArgument;
        vertices.createNew(fileName("Vertices"), pageSize);
        edges.createNew(fileName("Edges"), pageSize);
        edgesByVertex.createNew(fileName("EdgesByVertex"), pageSize);
    }

    // Access an existing undirected graph.
    void accessExisting(
        const string& baseNameArgument,
        bool readWriteAccess)
    {
        baseName = baseNameArgument;
        vertices.accessExisting(fileName("Vertices"), readWriteAccess);
        edges.accessExisting(fileName("Edges"), readWriteAccess);
        edgesByVertex.accessExisting(fileName("EdgesByVertex"), readWriteAccess);
    }
    void accessExistingReadOnly(const string& baseName)
    {
        accessExisting(baseName, false);
    }
    void accessExistingReadWrite(const string& baseName)
    {
        accessExisting(baseName, true);
    }

    bool isOpen()
    {
        return
            vertices.isOpen and
            edges.isOpen and
            edgesByVertex.isOpen();
    }

    bool isOpenWithWriteAccess()
    {
        return
            vertices.isOpenWithWriteAccess and
            edges.isOpenWithWriteAccess and
            edgesByVertex.isOpenWithWriteAccess();
    }



private:

    // The base name used for memory mapped files.
    // If this is empty, all memory mapped files are created
    // with an empty name, which means they are created
    // on anonymous memory.
    string baseName;

    // Given a data structure name, return the file name to be used.
    string fileName(const string& name) const
    {
        if(baseName.empty()) {
            return "";
        } else {
            return baseName + "-" + name;
        }
    }
};

#endif


