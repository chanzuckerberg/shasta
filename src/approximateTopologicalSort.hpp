#ifndef SHASTA_APPROXIMATE_TOPOLOGICAL_SORT_HPP
#define SHASTA_APPROXIMATE_TOPOLOGICAL_SORT_HPP

/*******************************************************************************

Approximate topological sort of a directed graph,
non necessarily acyclic.

The computed approximate topological sort of the given graph
is an exact topological sort for an acyclic graph
that uses a subset of the edges of the input graph.
Edges are added to this subset in the order specified
in the call, but edges that would cause a cycle are not added.

The topological sort of the acyclic graph is maintained dynamically
using the PK algorithm described in:
David. J. Pearce and Paul H. J. Kelly,
A Dynamic Topological Sort Algorithm for Directed Acyclic graphs,
ACM Journal of Experimental Algorithmics 11, 1-24 (2006).

Note that a C++ implementation of this algorithm
is available on GitHub:
https://github.com/DavePearce/DynamicTopologicalSort
However, this code was written from scratch using the algorithm
description in the paper.

This function is intrusive and requires the following
fields to exist in the vertices and edges of the input graph:
Vertices:
- color (an integer type), used for DFSs in the affected region.
- rank (an integer type), the computed topological sort order,
  beginning at 0.
Edge:
- isDagEdge (bool), set to true by this function for edges that contributed
  to the topological sort computations. The remaining
  edges have isDagEdge set to false. They were excluded because they caused cycles.

Note that the order in which the edges are processed determines which
edges are classified as causing cycles.
Therefore, the edges should be processed in decreasing order of "importance".
Edges processed first have lower chance of being flagged as causing cycles
and therefore excluded from the topological sort.

As an example, if the input graph consists of a single cycle,
Only the last edge processed will be classified as causing a cycle.

********************************************************************************/

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <stack>
#include "vector.hpp"

namespace shasta {
    template<class Graph> class ApproximateTopologicalSortEdgePredicate;
    template<class Graph> void approximateTopologicalSort(
        Graph&,
        const vector<typename Graph::edge_descriptor>&);
}



// Predicate class used by function approximateTopologicalSort.
template<class Graph> class shasta::ApproximateTopologicalSortEdgePredicate {
public:
    ApproximateTopologicalSortEdgePredicate(const Graph& graph) : graph(&graph) {}
    const Graph* graph;
    bool operator()(const typename Graph::edge_descriptor e) const
    {
        return (*graph)[e].isDagEdge;
    }
};



template<class Graph> void shasta::approximateTopologicalSort(
    Graph& graph,
    const vector<typename Graph::edge_descriptor>& edgesToProcess)
{
    // Type names for readability.
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;

    // Initialize.
    size_t rank = 0;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        auto& vertex = graph[v];
        vertex.color = 0;
        vertex.rank = rank++;
    }
    BGL_FORALL_EDGES_T(e, graph, Graph) {
        graph[e].isDagEdge = false;
    }

    // We will work with a filtered graph that includes only the edges
    // that we have already processed and that did not introduce cycles.
    using EdgePredicate = ApproximateTopologicalSortEdgePredicate<Graph>;
    const EdgePredicate edgePredicate(graph);
    using FilteredGraph = boost::filtered_graph<Graph, EdgePredicate>;
    FilteredGraph filteredGraph(graph, edgePredicate);

    // Vectors to hold deltaF and deltaB
    // (see definition 2.6 in the paper).
    // Stored as pairs (rank, vertex descriptor).
    // Defined here outside the loop to reduce memory allocation activity.
    vector< pair<size_t, vertex_descriptor> > deltaF;
    vector< pair<size_t, vertex_descriptor> > deltaB;

    // Vector to hold the ranks of vertices in deltaF and deltaB.
    vector<size_t> deltaRanks;

    // Stack used for DFS's.
    std::stack<vertex_descriptor> vertexStack;

    // Loop over the given edges in the specified order.
    for(const edge_descriptor e: edgesToProcess) {
        const vertex_descriptor vX = source(e, graph);
        const vertex_descriptor vY = target(e, graph);
        auto& vertexX = graph[vX];
        auto& vertexY = graph[vY];

        // If this edge is consistent with the existing rank,
        // just mark the edge as belonging to the DAG.
        // The topological sort is not affected.
        if(vertexX.rank < vertexY.rank) {
            graph[e].isDagEdge = true;
            continue;
        }

        deltaRanks.clear();

        // Use a forward DFS starting at vY, and limited to the affected region,
        // to compute deltaF (see definition 2.6 in the paper).
        deltaF.clear();
        if(!vertexStack.empty()) {
            throw runtime_error("Vertex stack is not empty.");
        }
        vertexStack.push(vY);
        deltaF.push_back(make_pair(vertexY.rank, vY));
        deltaRanks.push_back(vertexY.rank);
        graph[vY].color = 1;
        while(!vertexStack.empty()) {
            const vertex_descriptor v0 = vertexStack.top();
            vertexStack.pop();
            BGL_FORALL_OUTEDGES_T(v0, e01, filteredGraph, FilteredGraph) {
                const vertex_descriptor v1 = target(e01, graph);
                auto& vertex1 = graph[v1];
                if(vertex1.rank > vertexX.rank) {
                    // Outside the affected region.
                    continue;
                }
                if(vertex1.color == 1) {
                    // We already have this vertex.
                    continue;
                }
                vertex1.color = 1;
                vertexStack.push(v1);
                deltaF.push_back(make_pair(vertex1.rank, v1));
                deltaRanks.push_back(vertex1.rank);
            }
        }
        sort(deltaF.begin(), deltaF.end());
        for(const auto& p: deltaF) {
            graph[p.second].color = 0;
        }



        // Use a backward DFS starting at vX, and limited to the affected region,
        // to compute deltaB (see definition 2.6 in the paper).
        deltaB.clear();
        if(!vertexStack.empty()) {
            throw runtime_error("Vertex stack is not empty.");
        }
        vertexStack.push(vX);
        deltaB.push_back(make_pair(vertexX.rank, vX));
        deltaRanks.push_back(vertexX.rank);
        graph[vX].color = 1;
        while(!vertexStack.empty()) {
            const vertex_descriptor v0 = vertexStack.top();
            vertexStack.pop();
            BGL_FORALL_INEDGES_T(v0, e01, filteredGraph, FilteredGraph) {
                const vertex_descriptor v1 = source(e01, graph);
                auto& vertex1 = graph[v1];
                if(vertex1.rank < vertexY.rank) {
                    // Outside the affected region.
                    continue;
                }
                if(vertex1.color == 1) {
                    // We already have this vertex.
                    continue;
                }
                vertex1.color = 1;
                vertexStack.push(v1);
                deltaB.push_back(make_pair(vertex1.rank, v1));
                deltaRanks.push_back(vertex1.rank);
            }
        }
        sort(deltaB.begin(), deltaB.end());
        for(const auto& p: deltaB) {
            graph[p.second].color = 0;
        }


        // Sort the ranks. If we find any duplicates, deltaF and deltaB
        // intersect, which means this edge introduces a cycle.
        // So we ignore it.
        sort(deltaRanks.begin(), deltaRanks.end());
        bool createsCycle = false;
        for(size_t i=1; i<deltaRanks.size(); i++) {
            if(deltaRanks[i-1] == deltaRanks[i]) {
                createsCycle = true;
                break;
            }
        }
        if(!createsCycle) {

            // This edge does not create a cycle.
            // redistribute the deltaRanks to vertices in deltaB and deltaF.
            size_t i = 0;
            for(const auto& p: deltaB) {
               graph[p.second].rank = deltaRanks[i++];
            }
            for(const auto& p: deltaF) {
               graph[p.second].rank = deltaRanks[i++];
            }
            graph[e].isDagEdge = true;
        }

    }
}

#endif
