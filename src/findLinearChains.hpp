#ifndef SHASTA_FIND_LINEAR_CHAINS_HPP
#define SHASTA_FIND_LINEAR_CHAINS_HPP

// Find linear chains in a directed graph

#include <boost/graph/iteration_macros.hpp>
#include <list>
#include "vector.hpp"

namespace shasta {

    template<class Graph> void findLinearChains(
        const Graph&,
        vector< std::list<typename Graph::edge_descriptor> >&);

}


template<class Graph> inline void shasta::findLinearChains(
    const Graph& graph,
    vector< std::list<typename Graph::edge_descriptor> >& chains)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;

    // The edges we have already encountered.
    std::set<edge_descriptor> edgesFound;


    chains.clear();

    // Consider all possible start edges for the chain.
    BGL_FORALL_EDGES_T(eStart, graph, Graph) {

        // If we already assigned this edge to a chain, skip it.
        if(edgesFound.find(eStart) != edgesFound.end()) {
            continue;
        }

        // Add a new chain consisting of the start edge.
        chains.resize(chains.size() + 1);
        std::list<edge_descriptor>& chain = chains.back();
        chain.push_back(eStart);
        edgesFound.insert(eStart);

        // Extend forward.
        bool isCircular = false;
        edge_descriptor e = eStart;
        while(true) {
            const vertex_descriptor v = target(e, graph);
            if(in_degree(v, graph) != 1) {
                break;
            }
            if(out_degree(v, graph) != 1) {
                break;
            }
            BGL_FORALL_OUTEDGES_T(v, eNext, graph, Graph) {
                e = eNext;
                break;
            }
            if(e == eStart) {
                isCircular = true;
                break;
            }
            chain.push_back(e);
            SHASTA_ASSERT(edgesFound.find(e) == edgesFound.end());
            edgesFound.insert(e);
        }


        // Extend backward.
        if(not isCircular) {
            edge_descriptor e = eStart;
            while(true) {
                const vertex_descriptor v = source(e, graph);
                if(in_degree(v, graph) != 1) {
                    break;
                }
                if(out_degree(v, graph) != 1) {
                    break;
                }
                BGL_FORALL_INEDGES_T(v, ePrevious, graph, Graph) {
                    e = ePrevious;
                    break;
                }
                if(e == eStart) {
                    isCircular = true;
                    break;
                }
                chain.push_front(e);
                SHASTA_ASSERT(edgesFound.find(e) == edgesFound.end());
                edgesFound.insert(e);
            }
        }

    }

    // Check that all edges were found.
    SHASTA_ASSERT(edgesFound.size() == num_edges(graph));
}



#endif

