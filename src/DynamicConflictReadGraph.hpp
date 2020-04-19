#ifndef SHASTA_DYNAMIC_CONFLICT_READ_GRAPH_HPP
#define SHASTA_DYNAMIC_CONFLICT_READ_GRAPH_HPP

// A non-persistent representation of the conflict read graph
// that uses the Boost Graph library.
// Used by Assembler::cleanupConflictReadGraph.

// Shasta.
#include "ConflictReadGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <limits>

namespace shasta {
    class DynamicConflictReadGraph;
    class DynamicConflictReadGraphVertex;
    class DynamicConflictReadGraphEdge;

    using DynamicConflictReadGraphBaseClass =
        boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::undirectedS,
        DynamicConflictReadGraphVertex,
        DynamicConflictReadGraphEdge
        >;
}



class shasta::DynamicConflictReadGraphVertex {
public:
    ConflictReadGraph::VertexId vertexId;
    DynamicConflictReadGraphVertex(
        ConflictReadGraph::VertexId vertexId = std::numeric_limits<ConflictReadGraph::VertexId>::max()
        ) :
        vertexId(vertexId) {}

    OrientedReadId getOrientedReadId() const
    {
        return ConflictReadGraph::getOrientedReadId(vertexId);
    }

    // Coloring information.
    static const uint64_t invalid = std::numeric_limits<uint64_t>::max();
    uint64_t componentId = invalid;
    uint64_t color = invalid;
    bool isColored() const
    {
        return componentId != invalid and color != invalid;
    }
};



class shasta::DynamicConflictReadGraphEdge {
public:
    ConflictReadGraph::EdgeId edgeId;
    DynamicConflictReadGraphEdge(ConflictReadGraph::EdgeId edgeId) :
        edgeId(edgeId) {}
};



class shasta::DynamicConflictReadGraph : public DynamicConflictReadGraphBaseClass {
public:
    using VertexId = ConflictReadGraph::VertexId;
    using EdgeId = ConflictReadGraph::EdgeId;

    DynamicConflictReadGraph(const ConflictReadGraph&);
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

private:
    void writeGraphviz(ostream&, vertex_descriptor) const;
    void writeGraphviz(ostream&, edge_descriptor) const;
};

#endif

