#ifndef SHASTA_DYNAMIC_CONFLICT_READ_GRAPH_HPP
#define SHASTA_DYNAMIC_CONFLICT_READ_GRAPH_HPP

// A non-persistent representation of the conflict read graph
// that uses the Boost Graph library.
// Used by Assembler::cleanupConflictReadGraph.

// Shasta.
#include "ConflictReadGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

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
};

#endif

