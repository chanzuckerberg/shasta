#ifndef SHASTA_COMPACT_UNDIRECTED_GRAPH_HPP
#define SHASTA_COMPACT_UNDIRECTED_GRAPH_HPP



// Class CompactUndirectedGraph is used to represent
// an undirected graph with minimal memory allocation overhead.

// It implement some but not all API required by the Boost Graph library.
// Therefore it will not necessarily work with all
// graph algorithms of the Boost Graph library.

// This is similar to boost compressed_sparse_row_graph,
// which however does not support undirected graphs.

// Class CompactUndirectedGraph also provides a clear operation
// that allows reusing the same object multiple times
// with minimal memory allocation activity.

// Shasta.
#include "SHASTA_ASSERT.hpp"

// Boost Graph library.
#include <boost/graph/graph_selectors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "array.hpp"
#include "iostream.hpp"
#include  <limits>
#include "utility.hpp"
#include "vector.hpp"

// Forward definitions.
namespace ChanZuckerberg {
    namespace shasta {
        template<class Vertex, class Edge> class CompactUndirectedGraph;
        void testCompactUndirectedGraph1();
        void testCompactUndirectedGraph2();
    }
}



template<class Vertex, class Edge>
    class ChanZuckerberg::shasta::CompactUndirectedGraph {
public:

    // Type used by the vertex_descriptor and edge_descriptor.
    // This determines the maximum number of vertices and edges
    // that can be stored.
    using Int = size_t;



    // Strongly typed vertex descriptor. It simply stores the index
    // of a vertex in the vertices vector.
    // This is safer than just using an Int.
    class vertex_descriptor {
    public:
        // Index into the vertices vector.
        Int v;

        explicit vertex_descriptor(Int v = std::numeric_limits<Int>::max()) : v(v) {}
        bool operator==(const vertex_descriptor& that) const
        {
            return v == that.v;
        }
        bool operator!=(const vertex_descriptor& that) const
        {
            return v != that.v;
        }
    };
    static vertex_descriptor null_vertex()
    {
        return vertex_descriptor();
    }



    // An edge descriptor stores the index of the vertex in the edges array,
    // plus a flag that tells us whether the edge is being seen in
    // the same direction as stored or in the reverse direction.
    class edge_descriptor {
    public:

        // Index info the edges vector.
        Int e;

        // The vertex that the edge is seen as starting from.
        // If this equal null_vertex, the edge is returned as stored.
        vertex_descriptor v;

        explicit edge_descriptor(
            Int e = std::numeric_limits<Int>::max(),
            vertex_descriptor v = null_vertex()) :
            e(e), v(v) {}
    };



    // The state determines what operations are allowed.
    // States are reached in the order presented here.
    // Calling clear removes all vertices and edges and
    // puts the graph back in the AddingVertices state.
    enum class State {
        AddingVertices,
        AddingEdges,
        Processing
    };



    // Operations allowed in all states.
    State getState() const;
    Int vertexCount() const;
    Int edgeCount() const;
    Vertex& operator[](vertex_descriptor);
    const Vertex& operator[](vertex_descriptor) const;
    Edge& operator[](edge_descriptor);
    const Edge& operator[](edge_descriptor) const;

    // Operations allowed only when getState() == AddingVertices.
    vertex_descriptor addVertex(const Vertex& vertex = Vertex());
    void sortVertices();
    void doneAddingVertices();  // Transitions to state = AddingEdges.

    // Operations allowed only when getState() == AddingEdges.
    edge_descriptor addEdge(vertex_descriptor, vertex_descriptor, const Edge& edge = Edge());
    void doneAddingEdges();  // Transitions to state = Processing.

    // All remaining operations are only allowed when getState() == Processing.

    // Clear all vertices and edges and put the graph back in the AddingVertices state.
    void clear();
    CompactUndirectedGraph();

    // Vertices of an edge.
    vertex_descriptor source(edge_descriptor e) const
    {
        SHASTA_ASSERT(e.e < edgeTable.size());
        const EdgeInfo& edgeInfo = edgeTable[e.e];
        if(e.v == null_vertex()) {
            return edgeInfo.vertices[0];    // Return as stored.
        } else {
            SHASTA_ASSERT(e.v==edgeInfo.vertices[0] || e.v==edgeInfo.vertices[1]);
            return e.v;
        }
    }
    vertex_descriptor target(edge_descriptor e) const
    {
        SHASTA_ASSERT(e.e < edgeTable.size());
        const EdgeInfo& edgeInfo = edgeTable[e.e];
        if(e.v == null_vertex()) {
            return edgeInfo.vertices[1];    // Return as stored.
        } else {
            if(e.v==edgeInfo.vertices[0]) {
                return edgeInfo.vertices[1];
            } else {
                SHASTA_ASSERT(e.v == edgeInfo.vertices[1]);
                return edgeInfo.vertices[0];
            }
            return e.v;
        }
    }

    // Dump the data structures.
    void dump(ostream&) const;



    // Iteration over vertices.
    class vertex_iterator {
    public:
        Int v;
        vertex_iterator(Int v= std::numeric_limits<Int>::max()) : v(v) {}
        vertex_descriptor operator*() const
        {
            return vertex_descriptor(v);
        }
        vertex_iterator& operator++()
        {
           v++;
           return *this;
        }
        bool operator==(const vertex_iterator& that) const
        {
                return v == that.v;
        }
        bool operator!=(const vertex_iterator& that) const
        {
                return v != that.v;
        }
        vertex_iterator& operator=(const vertex_iterator& that)
        {
            v = that.v;
            return *this;
        }
        vertex_iterator(const vertex_iterator& that) : v(that.v) {}
    };
    vertex_iterator verticesBegin() const
    {
        return 0;
    }
    vertex_iterator verticesEnd() const
    {
        return vertexCount();
    }
    pair<vertex_iterator, vertex_iterator> allVertices() const
    {
        return make_pair(verticesBegin(), verticesEnd());
    }



    // Iteration over edges.
    // This always returns edges in the same direction as stored.
    class edge_iterator {
    public:
        Int e;
        explicit edge_iterator(Int e= std::numeric_limits<Int>::max()) : e(e) {}
        edge_descriptor operator*() const
        {
            // Always returns edges in the same direction as stored.
            return edge_descriptor(e);
        }
        edge_iterator& operator++()
        {
           e++;
           return *this;
        }
        bool operator==(const edge_iterator& that) const
        {
                return e == that.e;
        }
        bool operator!=(const edge_iterator& that) const
        {
                return e != that.e;
        }
        edge_iterator& operator=(const edge_iterator& that)
        {
            e = that.e;
            return *this;
        }
        edge_iterator(const edge_iterator& that) : e(that.e) {}
    };
    edge_iterator edgesBegin() const
    {
        return edge_iterator(0);
    }
    edge_iterator edgesEnd() const
    {
        return edge_iterator(edgeCount());
    }
    pair<edge_iterator, edge_iterator> allEdges() const
    {
        return make_pair(edgesBegin(), edgesEnd());
    }



    // Iteration over out-edges of a vertex.
    // This always returns edges with the source as the
    // vertex over which we are iterating.
    class out_edge_iterator {
    public:

        // Pointer into the edgeLists vector.
        const Int* edgeListsPointer;

        // The vertex that is seen as the source of the edge.
        // It his is null_vertex(), the edge is seen as stored.
        Int v;

        out_edge_iterator(const Int* edgeListsPointer, Int v) :
            edgeListsPointer(edgeListsPointer), v(v) {}
        edge_descriptor operator*() const
        {
            return edge_descriptor(*edgeListsPointer, vertex_descriptor(v));
        }
        out_edge_iterator& operator++()
        {
           edgeListsPointer++;
           return *this;
        }
        bool operator==(const out_edge_iterator& that) const
        {
                return edgeListsPointer==that.edgeListsPointer && v==that.v;
        }
        bool operator!=(const out_edge_iterator& that) const
        {
                return !(*this == that);
        }
        out_edge_iterator& operator=(const out_edge_iterator& that)
        {
            edgeListsPointer = that.edgeListsPointer;
            v = that.v;
            return *this;
        }
        out_edge_iterator(const out_edge_iterator& that) :
            edgeListsPointer(that.edgeListsPointer), v(that.v) {}
    };
    out_edge_iterator outEdgesBegin(vertex_descriptor v) const
    {
        return out_edge_iterator(edgeLists.data() + vertexTable[v.v].second, v.v);
    }
    out_edge_iterator outEdgesEnd(vertex_descriptor v) const
    {
        return out_edge_iterator(edgeLists.data() + vertexTable[v.v+1].second, v.v);
    }
    pair<out_edge_iterator, out_edge_iterator> allOutEdges(vertex_descriptor v) const
    {
        return make_pair(outEdgesBegin(v), outEdgesEnd(v));
    }


    // Traits required for compatibility with the Boost Graph library.
    using directed_category = boost::undirectedS;
    using edge_parallel_category = boost::allow_parallel_edge_tag;
    using traversal_category = boost::adjacency_graph_tag;

private:

    State state = State::AddingVertices;

    // The vertices, indexed by the vertex_descriptor.
    // If State==Processing, the Int is the index in
    // edgeLists of the first edge descriptor for each vertex.
    // There is also a dummy vertex at the end
    // that points to the one past the end of edgeLists.
    vector< pair<Vertex, Int> > vertexTable;

    // The edges, indexed by the edge descriptor.
    class EdgeInfo {
    public:
        array<vertex_descriptor, 2> vertices;
        Edge edge;
        EdgeInfo(
            vertex_descriptor v0,
            vertex_descriptor v1,
            const Edge& edge) :
            vertices(array<vertex_descriptor, 2>({v0, v1})),
            edge(edge)
        {
        }
    };
    vector<EdgeInfo> edgeTable;

    // The edges descriptors of each vertex
    // (or rather the underlying Int).
    // Each vertex points to the index in this vector
    // for the first edge of the vertex.
    vector<Int> edgeLists;
};



// Implement some of the boost adjacency_list API, but not all of it.
namespace ChanZuckerberg {
    namespace shasta {

        template<class Vertex, class Edge> std::pair<
            typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::vertex_iterator,
            typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::vertex_iterator
            > vertices(const ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>& graph)
        {
            return graph.allVertices();
        }

        template<class Vertex, class Edge> std::pair<
            typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::edge_iterator,
            typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::edge_iterator
            > edges(const ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>& graph)
        {
            return graph.allEdges();
        }

        template<class Vertex, class Edge>
            typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::vertex_descriptor
            source(
            typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::edge_descriptor e,
            const CompactUndirectedGraph<Vertex, Edge>& graph)
        {
            return graph.source(e);
        }

        template<class Vertex, class Edge>
            typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::vertex_descriptor
            target(
            typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::edge_descriptor e,
            const CompactUndirectedGraph<Vertex, Edge>& graph)
        {
            return graph.target(e);
        }
    }

        template<class Vertex, class Edge> std::pair<
            typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::out_edge_iterator,
            typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::out_edge_iterator
            > out_edges(
                typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::vertex_descriptor v,
                const ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>& graph)
        {
            return graph.allOutEdges(v);
        }
}



// Implementation follows.

template<class Vertex, class Edge>
    inline
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    CompactUndirectedGraph()
{
    clear();
}

template<class Vertex, class Edge>
    inline void
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    clear()
{
    state = State::AddingVertices;
    vertexTable.clear();
    edgeTable.clear();
    edgeLists.clear();
}

template<class Vertex, class Edge>
    inline typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::State
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    getState() const
{
    return state;
}

template<class Vertex, class Edge>
    inline typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::Int
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    vertexCount() const
{
    if(state == State::Processing) {
        return Int(vertexTable.size() -  1);
    } else {
        return Int(vertexTable.size());
    }
}

template<class Vertex, class Edge>
    inline typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::Int
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    edgeCount() const
{
    return Int(edgeTable.size());
}

template<class Vertex, class Edge>
    inline Vertex&
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    operator[](vertex_descriptor v)
{
    SHASTA_ASSERT(v.v < vertexTable.size());
    return vertexTable[v.v].first;
}

template<class Vertex, class Edge>
    inline const Vertex&
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    operator[](vertex_descriptor v) const
{
    SHASTA_ASSERT(v.v < vertexTable.size());
    return vertexTable[v.v].first;
}

template<class Vertex, class Edge>
    inline Edge&
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    operator[](edge_descriptor e)
{
    SHASTA_ASSERT(e.e < edgeTable.size());
    return edgeTable[e.e].edge;
}

template<class Vertex, class Edge>
    inline const Edge&
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    operator[](edge_descriptor e) const
{
    SHASTA_ASSERT(e.e < edgeTable.size());
    return edgeTable[e.e].edge;
}

template<class Vertex, class Edge>
    inline typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::vertex_descriptor
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    addVertex(const Vertex& vertex)
{
    SHASTA_ASSERT(state == State::AddingVertices);
    const vertex_descriptor v = vertex_descriptor(Int(vertexTable.size()));
    vertexTable.push_back(make_pair(vertex, Int(0)));
    return v;
}

template<class Vertex, class Edge>
    inline void
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    sortVertices()
{
    SHASTA_ASSERT(state == State::AddingVertices);
    sort(vertexTable.begin(), vertexTable.end());
}

template<class Vertex, class Edge>
    inline void
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    doneAddingVertices()
{
    state = State::AddingEdges;
}

template<class Vertex, class Edge>
    inline typename ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::edge_descriptor
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    addEdge(
    vertex_descriptor v0,
    vertex_descriptor v1,
    const Edge& edge)
{
    SHASTA_ASSERT(state == State::AddingEdges);
    const edge_descriptor e = edge_descriptor(Int(edgeTable.size()));
    edgeTable.push_back(EdgeInfo(v0, v1, edge));
    return e;
}

template<class Vertex, class Edge>
    inline void
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    doneAddingEdges()
{
    // First, we store the degree of each vertex
    // (number of edges to the vertex).
    for(const EdgeInfo& edgeInfo: edgeTable) {
        const vertex_descriptor v0 = edgeInfo.vertices[0];
        const vertex_descriptor v1 = edgeInfo.vertices[1];
        ++(vertexTable[v0.v].second);
        ++(vertexTable[v1.v].second);
    }


    // Now accumulate to make each vertex point to one past the last
    // edge of the vertex.
    Int n = Int(0);
    for(Int v=0; v<Int(vertexTable.size()); v++) {
        Int& i = vertexTable[v].second;
        n += i;
        i = n;
    }
    vertexTable.push_back(make_pair(Vertex(), n));

    // Now store in reverse order.
    edgeLists.resize(n);
    for(Int e=0; e<Int(edgeTable.size()); e++) {
        const EdgeInfo& edgeInfo = edgeTable[e];
        const vertex_descriptor v0 = edgeInfo.vertices[0];
        const vertex_descriptor v1 = edgeInfo.vertices[1];
        edgeLists[--vertexTable[v0.v].second] = e;
        edgeLists[--vertexTable[v1.v].second] = e;
    }
    SHASTA_ASSERT(vertexTable.front().second == Int(0));

    // Reverse the order, for each vertex,
    // so it is sorted.
    for(Int v=0; v<Int(vertexTable.size()-1); v++) {
        std::reverse(
            edgeLists.begin() + vertexTable[v].second,
            edgeLists.begin() + vertexTable[v + Int(1)].second);
    }

    SHASTA_ASSERT(edgeLists.size() == 2*edgeTable.size());
    SHASTA_ASSERT(vertexTable.back().second == Int(edgeLists.size()));


    state = State::Processing;
}

template<class Vertex, class Edge>
    inline void
    ChanZuckerberg::shasta::CompactUndirectedGraph<Vertex, Edge>::
    dump(ostream& s) const
{
    s << vertexCount() << " vertices, ";
    s << edgeCount() << " edges." << endl;

    s << "Edges:\n";
    for(Int e=0; e<Int(edgeTable.size()); e++) {
        s << e << " " << edgeTable[e].vertices[0].v << " ";
        s << edgeTable[e].vertices[1].v << "\n";
    }

    s << "Indices of first/last edge of each vertex:\n";
    for(Int v=0; v<Int(vertexTable.size()); v++) {
        s << v << " " << vertexTable[v].second << "\n";
    }

    s << "Edge lists vector:\n";
    for(size_t i=0; i<edgeLists.size(); i++) {
        s << i << ": " << edgeLists[i] << "\n";
    }
}

#endif
