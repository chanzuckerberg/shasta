#ifndef SHASTA_BUBBLE_GRAPH_HPP
#define SHASTA_BUBBLE_GRAPH_HPP

// Class BubbleGraph is used in AssemblyGraph::findBubbleChains.
// Each vertex corresponds to an assembly graph vertex.
// Each edge corresponds to either:
// - A bubble.
// - An assembly graph edge v0->v1 that satisfies the following:
//   * Does not belong to a bubble.
//   * v1 is the source vertex of a bubble, or
//     v0 is the target vertex of a bubble.

// Shasta.
#include "AssemblyGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <list>
#include <map>
#include <set>
#include "string.hpp"

namespace shasta {
    class BubbleGraph;
    class BubbleGraphEdge;
    class BubbleGraphVertex;
}



class shasta::BubbleGraphVertex {
public:
    AssemblyGraph::VertexId vertexId;
};



class shasta::BubbleGraphEdge {
public:
    bool isBubble;
    AssemblyGraph::EdgeId edgeId;
    uint64_t bubbleId;  // Index into the bubbles vector.

    // The linear chain this edge belongs to.
    uint64_t linearChainId = std::numeric_limits<uint64_t>::max();
};



class shasta::BubbleGraph :
    public boost::adjacency_list<
    boost::listS,
    boost::listS,
    boost::bidirectionalS,
    BubbleGraphVertex,
    BubbleGraphEdge
    > {
public:
    BubbleGraph(const AssemblyGraph&);

    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::VertexId;
    using Bubble = AssemblyGraph::Bubble;
    using BubbleId = uint64_t;

private:

    void addBubbles(const AssemblyGraph&);
    void addBubble(BubbleId, const Bubble&, const AssemblyGraph&);
    void addEdges(const AssemblyGraph&);

    // Map that gives the vertex_descriptor corresponding to a VertexId.
    std::map<VertexId, vertex_descriptor> vertexMap;

    // The vertices that are sources or targets of a bubble.
    std::set<VertexId> bubbleSources;
    std::set<VertexId> bubbleTargets;

    // Return the vertex_descriptor corresponding to
    // a vertex id, creating the vertex if necessary.
    vertex_descriptor getVertex(VertexId);

    // Assembly graph edges that belong to a bubble.
    std::set<EdgeId> bubbleEdges;

    // Linear chains (paths) in the bubble graph. Each edge
    // belongs to exactly one linear chain.
    class LinearChain {
    public:
        std::list<edge_descriptor> edges;
        bool isCircular = false;
    };
    vector<LinearChain> linearChains;
public:
    void findLinearChains();
    void writeLinearChains(const string& fileName, const AssemblyGraph&);
    void writeLinearChains(ostream&, const AssemblyGraph&) const;
};

#endif
