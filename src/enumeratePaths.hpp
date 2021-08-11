#ifndef SHASTA_ENUMERATE_PATHS_HPP
#define SHASTA_ENUMERATE_PATHS_HPP

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

#include "algorithm.hpp"
#include <stack>
#include "tuple.hpp"
#include "vector.hpp"

namespace shasta {

template<class G> void enumerateSelfAvoidingPaths(const G&,
    typename G::vertex_descriptor vA, typename G::vertex_descriptor vB,
    vector<vector<typename G::edge_descriptor> > &paths);

}

// Enumerate self-avoiding paths starting at v0 and ending at v1.
// Self-avoiding means that an edge cannot be used twice.
template<class G> void shasta::enumerateSelfAvoidingPaths(const G &g,
    typename G::vertex_descriptor vA, typename G::vertex_descriptor vB,
    vector<vector<typename G::edge_descriptor> > &paths)
{
    using vertex_descriptor = typename G::vertex_descriptor;
    using edge_descriptor = typename G::edge_descriptor;
    using out_edge_iterator = typename G::out_edge_iterator;
    using Path = vector<typename G::edge_descriptor>;
    using std::stack;

    paths.clear();
    stack<Path> partialPaths;

    // For some reason I was not able to get BGL_FORALL_OUTEDGES_T 
    // to work, so using explicit iterators instead.
    out_edge_iterator it, end;

    std::tie(it, end) = boost::out_edges(vA, g);
    for (; it != end; ++it) {
        const edge_descriptor e = *it;
        const vertex_descriptor v1 = boost::target(e, g);
        const Path path(1, e);
        if (v1 == vB) {
            paths.push_back(path);
        } else {
            partialPaths.push(path);
        }
    }

    // Recursively add edges to the partial paths.
    while (not partialPaths.empty()) {
        const Path path = partialPaths.top();
        partialPaths.pop();

        const vertex_descriptor v0 = boost::target(path.back(), g);
        std::tie(it, end) = boost::out_edges(v0, g);
        for (; it != end; ++it) {
            const edge_descriptor e = *it;

            // If e is already on the path, skip it.
            if (find(path.begin(), path.end(), e) != path.end()) {
                continue;
            }

            const vertex_descriptor v1 = target(e, g);
            Path newPath = path;
            newPath.push_back(e);
            if (v1 == vB) {
                paths.push_back(newPath);
            } else {
                partialPaths.push(newPath);
            }
        }
    }
}

#endif

