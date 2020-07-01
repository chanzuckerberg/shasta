#ifndef SHASTA_FIND_LINEAR_CHAINS_HPP
#define SHASTA_FIND_LINEAR_CHAINS_HPP

// Find linear chains in a directed graph

#include <boost/graph/iteration_macros.hpp>
#include <list>
#include "vector.hpp"

namespace shasta {

    // Find linear chains of edges (paths).
    template<class Graph> void findLinearChains(
        const Graph&,
        vector< std::list<typename Graph::edge_descriptor> >&);


    // Find linear chains of vertices.
    template<class Graph> void findLinearVertexChains(
        const Graph&,
        vector< std::list<typename Graph::vertex_descriptor> >&);
    template<class Graph> void findLinearVertexChains(
        const Graph&,
        vector< vector<typename Graph::vertex_descriptor> >&);
}



// Find linear chains of edges (paths).
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



// Find linear chains of vertices.
template<class Graph> void shasta::findLinearVertexChains(
    const Graph& graph,
    vector< std::list<typename Graph::vertex_descriptor> >& chains)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;

    // The vertices we have already encountered.
    std::set<vertex_descriptor> verticesFound;

    chains.clear();

    // Consider all possible start vertices for the chain.
    BGL_FORALL_VERTICES_T(vStart, graph, Graph) {

        // If we already assigned this vertex to a chain, skip it.
        if(verticesFound.find(vStart) != verticesFound.end()) {
            continue;
        }

        // Add a new chain consisting of the start vertex.
        chains.resize(chains.size() + 1);
        std::list<vertex_descriptor>& chain = chains.back();
        chain.push_back(vStart);
        verticesFound.insert(vStart);

        // Extend forward.
        bool isCircular = false;
        vertex_descriptor v = vStart;
        while(true) {
            if(out_degree(v, graph) != 1) {
                break;
            }
            BGL_FORALL_OUTEDGES_T(v, e, graph, Graph) {
                v = target(e, graph);
                break;
            }
            if(v == vStart) {
                isCircular = true;
                break;
            }
            if(in_degree(v, graph) != 1) {
                break;
            }
            chain.push_back(v);
            SHASTA_ASSERT(verticesFound.find(v) == verticesFound.end());
            verticesFound.insert(v);
        }

        // Extend backward.
        if(not isCircular) {
            vertex_descriptor v = vStart;
            while(true) {
                if(in_degree(v, graph) != 1) {
                    break;
                }
                BGL_FORALL_INEDGES_T(v, e, graph, Graph) {
                    v = source(e, graph);
                    break;
                }
                if(out_degree(v, graph) != 1) {
                    break;
                }
                chain.push_front(v);
                SHASTA_ASSERT(verticesFound.find(v) == verticesFound.end());
                verticesFound.insert(v);
            }
        }
    }



    // Check that all vertices were found.
    SHASTA_ASSERT(verticesFound.size() == num_vertices(graph));
}



template<class Graph> void shasta::findLinearVertexChains(
    const Graph& graph,
    vector< vector<typename Graph::vertex_descriptor> >& chains)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;

    // Find the chains.
    vector< std::list<typename Graph::vertex_descriptor> > chainLists;
    findLinearVertexChains(graph, chainLists);

    // Copy lists to vectors.
    chains.clear();
    chains.reserve(chainLists.size());
    for(const auto& chain: chainLists) {
        chains.push_back(vector<vertex_descriptor>());
        copy(chain.begin(), chain.end(), back_inserter(chains.back()));
    }
}


#endif

