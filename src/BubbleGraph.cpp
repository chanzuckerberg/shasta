#include "BubbleGraph.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>

#include "fstream.hpp"
#include <numeric>



BubbleGraph::BubbleGraph(const AssemblyGraph& assemblyGraph)
{
    addBubbles(assemblyGraph);
    addEdges(assemblyGraph);
}


BubbleGraph::vertex_descriptor BubbleGraph::getVertex(VertexId vertexId)
{
    const auto it = vertexMap.find(vertexId);
    if(it == vertexMap.end()) {
        const vertex_descriptor v = add_vertex(*this);
        vertexMap.insert(make_pair(vertexId, v));
        return v;
    } else {
        return it->second;
    }
}


void BubbleGraph::addBubbles(const AssemblyGraph& assemblyGraph)
{
    for(BubbleId bubbleId=0; bubbleId<assemblyGraph.bubbles.size(); bubbleId++) {
        const AssemblyGraph::Bubble& bubble = assemblyGraph.bubbles[bubbleId];
        addBubble(bubbleId, bubble, assemblyGraph);
    }
}



void BubbleGraph::addBubble(
    BubbleId bubbleId,
    const Bubble& bubble,
    const AssemblyGraph& assemblyGraph)
{
    BubbleGraph& graph = *this;

    // Locate the vertices of this bubble.
    const VertexId vertexId0 = bubble.v0;
    const VertexId vertexId1 = bubble.v1;

    const vertex_descriptor v0 = getVertex(vertexId0);
    const vertex_descriptor v1 = getVertex(vertexId1);

    bubbleSources.insert(vertexId0);
    bubbleTargets.insert(vertexId1);

    // Add the edge corresponding to this bubble.
    bool edgeWasAdded = false;
    edge_descriptor e;
    tie(e, edgeWasAdded) = add_edge(v0, v1, graph);
    SHASTA_ASSERT(edgeWasAdded);
    BubbleGraphEdge& edge = graph[e];
    edge.bubbleId = bubbleId;
    edge.isBubble = true;
    edge.edgeId = std::numeric_limits<EdgeId>::max();

    // Keep track of the bubble edges.
    for(EdgeId edgeId: assemblyGraph.edgesBySource[vertexId0]) {
        bubbleEdges.insert(edgeId);
    }
}



void BubbleGraph::addEdges(const AssemblyGraph& assemblyGraph)
{
    BubbleGraph& graph = *this;

    // Loop over edges of the assembly graph.
    for(EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {

        // If this edge belongs to a bubble, skip it.
        if(bubbleEdges.find(edgeId) != bubbleEdges.end()) {
            continue;
        }

        const AssemblyGraph::Edge& assemblyGraphEdge = assemblyGraph.edges[edgeId];
        const VertexId vertexId0 = assemblyGraphEdge.source;
        const VertexId vertexId1 = assemblyGraphEdge.target;

        // If vertexId0 is not a bubble target and vertexId1 is not a bubble source,
        // skip it.
        if(
            bubbleTargets.find(vertexId0) == bubbleTargets.end() and
            bubbleSources.find(vertexId1) == bubbleSources.end()) {
            continue;
        }

        // Add the edge.
        bool edgeWasAdded = false;
        edge_descriptor e;
        tie(e, edgeWasAdded) = add_edge(getVertex(vertexId0), getVertex(vertexId1), graph);
        SHASTA_ASSERT(edgeWasAdded);
        BubbleGraphEdge& edge = graph[e];
        edge.bubbleId = std::numeric_limits<BubbleId>::max();
        edge.isBubble = false;
        edge.edgeId = edgeId;
    }

}



void BubbleGraph::findLinearChains()
{
    BubbleGraph& graph = *this;

    // Loop over possible start edges.
    uint64_t chainId = 0;
    BGL_FORALL_EDGES(startEdge, graph, BubbleGraph) {

        // If this edge is already part of a chain, skip it.
        if(graph[startEdge].linearChainId != std::numeric_limits<uint64_t>::max()) {
            continue;
        }

        // Add a new chain consisting of the start edge.
        linearChains.push_back(LinearChain());
        LinearChain& chain = linearChains.back();
        chain.edges.push_back(startEdge);
        graph[startEdge].linearChainId = chainId;

        // Extend forward.
        edge_descriptor e = startEdge;
        while(true) {
            const vertex_descriptor v = target(e, graph);
            if(out_degree(v, graph) != 1) {
                break;
            }
            BGL_FORALL_OUTEDGES(v, eNext, graph, BubbleGraph) {
                e = eNext;
                break;
            }
            if(e == startEdge) {
                chain.isCircular = true;
                break;
            }
            chain.edges.push_back(e);
            graph[e].linearChainId = chainId;
        }


        // Extend backward.
        if(not chain.isCircular) {
            edge_descriptor e = startEdge;
            while(true) {
                const vertex_descriptor v = source(e, graph);
                if(in_degree(v, graph) != 1) {
                    break;
                }
                BGL_FORALL_INEDGES(v, ePrevious, graph, BubbleGraph) {
                    e = ePrevious;
                    break;
                }
                if(e == startEdge) {
                    chain.isCircular = true;
                    break;
                }
                chain.edges.push_front(e);
                graph[e].linearChainId = chainId;
            }

        }

        // Prepare for the next chain.
        ++chainId;

    }


    // Check that we added all the edges.
    BGL_FORALL_EDGES(e, graph, BubbleGraph) {
        SHASTA_ASSERT(graph[e].linearChainId != std::numeric_limits<uint64_t>::max());
    }
}



void BubbleGraph::writeLinearChains(
    const string& fileName,
    const AssemblyGraph& assemblyGraph)
{
    ofstream csv(fileName);
    writeLinearChains(csv, assemblyGraph);
};
void BubbleGraph::writeLinearChains(
    ostream& csv,
    const AssemblyGraph& assemblyGraph) const
{
    const BubbleGraph& graph = *this;

    // Find the maximum ploidy of the bubbles.
    uint64_t maximumPloidy = 0;
    for(const Bubble& bubble: assemblyGraph.bubbles) {
        const VertexId v0 = bubble.v0;
        const uint64_t ploidy = assemblyGraph.edgesBySource[v0].size();
        maximumPloidy = max(maximumPloidy, ploidy);
    }


    csv << "Chain,Circular,Position,";
    for(uint64_t ploidy=0; ploidy<maximumPloidy; ploidy++) {
        csv << "Segment" << ploidy << ",";
    }
    csv << "\n";


    for(uint64_t chainId=0; chainId<linearChains.size(); chainId++) {
        const LinearChain& chain = linearChains[chainId];
        uint64_t position=0;
        for(const edge_descriptor e: chain.edges) {
            const BubbleGraphEdge& edge = graph[e];
            csv << chainId << ",";
            csv << (chain.isCircular ? "Yes," : "No,");
            csv << position++ << ",";
            if(edge.isBubble) {
                const BubbleId bubbleId = edge.bubbleId;
                const Bubble& bubble = assemblyGraph.bubbles[bubbleId];
                const VertexId v0 = bubble.v0;
                for(EdgeId edgeId: assemblyGraph.edgesBySource[v0]) {
                    csv << edgeId << ",";
                }

            } else {
                csv << edge.edgeId << ",";
            }
            csv << "\n";
        }
    }
}
