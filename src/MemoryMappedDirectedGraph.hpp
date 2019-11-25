#ifndef SHASTA_MEMORY_MAPPED_DIRECTED_GRAPH_HPP
#define SHASTA_MEMORY_MAPPED_DIRECTED_GRAPH_HPP

// Shasta
#include "deduplicate.hpp"
#include "MemoryMappedVector.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

// Standard library.
#include "chrono.hpp"
#include <limits>
#include <map>
#include <queue>
#include "string.hpp"

// A directed graph stored on memory mapped data structures.
// It allows parallel edges (more than one edge
// between a pair of vertices).
// It has the following restrictions:
// - Vertices and edges can only be added, not removed.
// - When an edge is added, the edgesBySource
//   and edgesByTarget data structure, become invalid
//   and must be recreated with a call to

namespace shasta {
    namespace MemoryMapped {
        template<class Vertex, class Edge> class DirectedGraph;
    }
}



template<class Vertex, class Edge> class shasta::MemoryMapped::DirectedGraph {
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
    class EdgeInformation {
    public:
        VertexId v0;
        VertexId v1;
        Edge edge;
        EdgeInformation(VertexId v0, VertexId v1, const Edge& edge) :
            v0(v0), v1(v1), edge(edge) {}
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
    VertexId source(EdgeId edgeId) const {
        return edges[edgeId].v0;
    }
    VertexId target(EdgeId edgeId) const {
        return edges[edgeId].v1;
    }

    // The list of edges that have each vertex as their source
    // or their target. Indexed by a VertexId.
    // That is, edgeBySource[vertexId] contains the EdgeId's
    // of all edges with source vertexId,
    // and edgeByTarget[vertexId] contains the EdgeId's
    // of all edges with target vertexId.
    // These vectors become invalid whenever am edge is added
    // And they have to be recreated with a call to
    // computeConnectivity.
    VectorOfVectors<VertexId, uint64_t> edgesBySource;
    VectorOfVectors<VertexId, uint64_t> edgesByTarget;

    // Accessors for edges by source and by target.
    MemoryAsContainer<EdgeId> outEdges(VertexId vertexId)
    {
        return edgesBySource[vertexId];
    }
    MemoryAsContainer<EdgeId> inEdges(VertexId vertexId)
    {
        return edgesByTarget[vertexId];
    }
    uint64_t outDegree(VertexId vertexId)
    {
        return edgesBySource.size(vertexId);
    }
    uint64_t inDegree(VertexId vertexId)
    {
        return edgesByTarget.size(vertexId);
    }
    uint64_t totalDegree(VertexId vertexId)
    {
        return inDegree(vertexId) + outDegree(vertexId);
    }


    // A call to this function recomputes from scratch
    // the edgesBySource and edgesByTarget data structures.
    // This could be made multithreaded.
    void computeConnectivity()
    {
        // Clear out edgesBySource and edgesByTarget data structures.
        edgesBySource.clear();
        edgesByTarget.clear();

        // Pass 1: count the number of edges that each vertex
        // is the source/target of.
        edgesBySource.beginPass1(vertices.size());
        edgesByTarget.beginPass1(vertices.size());
        for(const EdgeInformation& edge: edges) {
            edgesBySource.incrementCount(edge.v0);
            edgesByTarget.incrementCount(edge.v1);
        }

        // Pass 2: store the edges that each vertex
        // is the source/target of.
        edgesBySource.beginPass2();
        edgesByTarget.beginPass2();
        for(EdgeId edgeId=0; edgeId<edges.size(); edgeId++) {
            const EdgeInformation& edge = edges[edgeId];
            edgesBySource.store(edge.v0, edgeId);
            edgesByTarget.store(edge.v1, edgeId);
        }
        edgesBySource.endPass2();
        edgesByTarget.endPass2();

        // Pass 3: sort the egde ids for each vertex.
        for(VertexId vertexId=0; vertexId<vertices.size(); vertexId++) {
            MemoryAsContainer<EdgeId> s = outEdges(vertexId);
            sort(s.begin(), s.end());
            MemoryAsContainer<EdgeId> t = inEdges(vertexId);
            sort(t.begin(), t.end());
        }

    }



    // The default constructor puts the graph in an uninitialized state.
    // The graph can only be used after calling createNew or accessExisting.
    DirectedGraph() {}

    // Create a new DirectedGraph.
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
        edgesBySource.createNew(fileName("EdgesBySource"), pageSize);
        edgesByTarget.createNew(fileName("EdgesByTarget"), pageSize);
    }

    // Access an existing Directed graph.
    void accessExisting(
        const string& baseNameArgument,
        bool readWriteAccess)
    {
        baseName = baseNameArgument;
        vertices.accessExisting(fileName("Vertices"), readWriteAccess);
        edges.accessExisting(fileName("Edges"), readWriteAccess);
        edgesBySource.accessExisting(fileName("EdgesBySource"), readWriteAccess);
        edgesByTarget.accessExisting(fileName("EdgesByTarget"), readWriteAccess);
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
            edgesBySource.isOpen() and
            edgesByTarget.isOpen();
    }

    bool isOpenWithWriteAccess()
    {
        return
            vertices.isOpenWithWriteAccess and
            edges.isOpenWithWriteAccess and
            edgesBySource.isOpenWithWriteAccess() and
            edgesByTarget.isOpenWithWriteAccess();
    }



    // Find the neighborhood of a vertex, allowing only edges
    // permitted by the edge filter object.
    class AbstractEdgeFilter {
    public:
        virtual bool allowEdge(EdgeId, const Edge&) const = 0;
    };
    bool findNeighborhood(
        VertexId vStart,
        uint64_t maxDistance,
        const AbstractEdgeFilter& edgeFilter,
        bool forward,   // True if allowed to move forward.
        bool backward,  // True if allowed to move backward.
        double timeout,
        std::map<VertexId, uint64_t>& neighbors
        )
    {
        const auto startTime = steady_clock::now();
        const bool debug = false;

        // Initialize the BFS.
        std::queue<VertexId> q;
        neighbors.clear();
        neighbors.insert(make_pair(vStart, 0));
        q.push(vStart);

        if(debug) {
            cout << "BFS begins at " << OrientedReadId(OrientedReadId::Int(vStart)) <<
                " maximum distance " << maxDistance << endl;
        }

        // Do the BFS.
        while(!q.empty()) {

            // See if we exceeded the timeout.
            if(timeout>0. and (seconds(steady_clock::now() - startTime) > timeout)) {
                neighbors.clear();
                return false;
            }

            // Dequeue a vertex.
            const VertexId v0 = q.front();
            q.pop();
            const uint64_t distance0 = neighbors[v0];
            const uint64_t distance1 = distance0 + 1;

            if(debug) {
                cout << "Dequeued " << OrientedReadId(OrientedReadId::Int(v0)) <<
                    " at distance " << distance0 <<
                    " out-degree " << outDegree(v0) <<
                    " in-degree " << inDegree(v0) << endl;
            }

            // Move forward.
            if(forward) {
                for(const EdgeId& edgeId: outEdges(v0)) {
                    if(not edgeFilter.allowEdge(edgeId, getEdge(edgeId))) {
                        continue;
                    }
                    const VertexId v1 = target(edgeId);
                    if(neighbors.find(v1) != neighbors.end()) {
                        // We already encountered v1.
                        continue;
                    }
                    neighbors.insert(make_pair(v1, distance1));
                    if(distance1 < maxDistance) {
                        q.push(v1);
                    }
                }
            }

            // Move backward.
            if(backward) {
                for(const EdgeId& edgeId: inEdges(v0)) {
                    if(not edgeFilter.allowEdge(edgeId, getEdge(edgeId))) {
                        continue;
                    }
                    const VertexId v1 = source(edgeId);
                    if(neighbors.find(v1) != neighbors.end()) {
                        // We already encountered v1.
                        continue;
                    }
                    neighbors.insert(make_pair(v1, distance1));
                    if(distance1 < maxDistance) {
                        q.push(v1);
                    }
                }
            }
        }

        return true;
    }



    // Compute the immediate (distance 1) neighborhood of two vertices
    // and their intersecion and union.
    void compareNeighborhoods(
        VertexId v0,
        VertexId v1,
        vector<VertexId>& neighbors0,
        vector<VertexId>& neighbors1,
        vector<VertexId>& intersectionVertices,
        vector<VertexId>& unionVertices)
    {
        // Compute neighbors of v0.
        neighbors0.clear();
        for(const EdgeId edgeId: outEdges(v0)) {
            const VertexId v = target(edgeId);
            if(v!=v0 && v!=v1) {
                neighbors0.push_back(target(edgeId));
            }
        }
        for(const EdgeId edgeId: inEdges(v0)) {
            const VertexId v = source(edgeId);
            if(v!=v0 && v!=v1) {
                neighbors0.push_back(target(edgeId));
            }
        }
        deduplicate(neighbors0);

        // Compute neighbors of v1.
        neighbors1.clear();
        for(const EdgeId edgeId: outEdges(v1)) {
            const VertexId v = target(edgeId);
            if(v!=v0 && v!=v1) {
                neighbors1.push_back(target(edgeId));
            }
        }
        for(const EdgeId edgeId: inEdges(v1)) {
            const VertexId v = source(edgeId);
            if(v!=v0 && v!=v1) {
                neighbors1.push_back(target(edgeId));
            }
        }
        deduplicate(neighbors1);


        // Compute intersection and union.
        intersectionVertices.clear();
        std::set_intersection(
            neighbors0.begin(), neighbors0.end(),
            neighbors1.begin(), neighbors1.end(),
            back_inserter(intersectionVertices));
        unionVertices.clear();
        std::set_union(
            neighbors0.begin(), neighbors0.end(),
            neighbors1.begin(), neighbors1.end(),
            back_inserter(unionVertices));
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


