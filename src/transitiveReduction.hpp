#ifndef SHASTA_TRANSITIVE_REDUCTION_HPP
#define SHASTA_TRANSITIVE_REDUCTION_HPP

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/topological_sort.hpp>

// Standard library.
#include "iterator.hpp"
using std::back_inserter;
#include <queue>
#include "vector.hpp"

namespace shasta {
    template<class Graph> void transitiveReduction(Graph&);
}



// Transitive reduction of a directed graph without cycles.
// Class Graph must be a boost::adjacency_list with
// the first three template arguments set to <listS, vecS, directedS>.
// If the graph has cycles, this throws boost::not_a_dag.
template<class Graph> void shasta::transitiveReduction(Graph& graph)
{
    using namespace boost;
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;
    
    // Check the Graph type.
    // Use C++20 concepts instead.
    static_assert(
        std::is_same<typename Graph::out_edge_list_selector, listS>::value, 
        "shasta::transitiveReduction requires an adjacency_list "
        "with the first template argument set to boost::listS.");
    static_assert(
        std::is_same<typename Graph::vertex_list_selector, vecS>::value, 
        "shasta::transitiveReduction requires an adjacency_list "
        "with the second template argument set to boost::vecS.");
    static_assert(
        std::is_same<typename Graph::directed_selector, bidirectionalS>::value, 
        "shasta::transitiveReduction requires an adjacency_list "
        "with the third template argument set to boost::bidirectionalS.");
        
    // Use boost topological_sort to get a vector of vertex descriptors
    // in reverse toplogical order.        
    vector<vertex_descriptor> sortedVertices;
    topological_sort(graph, back_inserter(sortedVertices));
    
    // Now construct a vector containg the rank of each vertex in topological order.
    vector<uint64_t> vertexRank(num_vertices(graph));
    uint64_t rank = num_vertices(graph);
    for(const vertex_descriptor v: sortedVertices) {
        vertexRank[v] = --rank;
    }
    
    
    
    // Find the edges that should be removed.
    vector<edge_descriptor> edgesToBeRemoved;
    vector<bool> wasVisited(num_vertices(graph));
    BGL_FORALL_EDGES_T(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        
        // Edge e should be removed if there is a path
        // from v0 to v1 that does not use edge e.
        
        
        
        // Do a forward BFS starting at v0 and ending at v1 but:
        // - Don't use edge e in the BFS.
        // - Don't use any vertices that have topological order 
        //   greater than the topological order of v1,
        //   because there can be no paths ending at v1 
        //   that use these vertices.
        // If the BFS encounters v1, edge e can be removed.
        
        // Initialize the BFS.
        std::queue<vertex_descriptor> q;
        q.push(v0);
        std::fill(wasVisited.begin(), wasVisited.end(), false);
        wasVisited[v0] = true;
        
        // BFS loop.
        while(not q.empty()) {
        
            // Dequeue a vertex.
            const vertex_descriptor vv0 = q.front();
            q.pop();
            
            // Loop over its out-edges.
            BGL_FORALL_OUTEDGES_T(vv0, ee, graph, Graph) {
            
                // Don't use edge e in the BFS.
                if(ee == e) {
                    continue;
                }
                
                // Get the other vertex in edge ee.
                const vertex_descriptor vv1 = target(ee, graph);
                
                // If vv1 was already visited in this BFS, skip it.
                if(wasVisited[vv1]) {
                    continue;
                }
                
                // If vv1 follows v1 in topological order, skip it. 
                if(vertexRank[vv1] > vertexRank[v1]) {
                    continue;
                }
                
                if(vv1 == v1) {
                    // We reached v1. Edge e can be removed and we can stop the BFS.
                    edgesToBeRemoved.push_back(e);
                    q = {};
                } else {
                    // Continue the BFS.
                    wasVisited[vv1] = true;
                    q.push(vv1);
                }
            }
        }
    }
    
    
    
    // Remove the edges.
    for(const edge_descriptor e: edgesToBeRemoved) {
        remove_edge(e, graph);
    }
}



#endif
