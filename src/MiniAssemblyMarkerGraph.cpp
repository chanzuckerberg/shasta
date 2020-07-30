#include "MiniAssemblyMarkerGraph.hpp"
#include "findLinearChains.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>



void MiniAssemblyMarkerGraph::removeSelfEdges()
{
    Graph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, Graph) {
        if(source(e,graph) == target(e, graph)) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }

}



void MiniAssemblyMarkerGraph::removeIsolatedVertices()
{
    Graph& graph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph)
    {
        if(in_degree(v, graph)==0 and out_degree(v, graph)==0) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        remove_vertex(v, graph);
    }

}



void MiniAssemblyMarkerGraph::removeLowCoverageEdges(
    uint64_t minTotalEdgeCoverage,
    uint64_t minPerStrandEdgeCoverage)
{
    Graph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, Graph) {

        // Check total coverage.
        if(graph[e].coverage() < minTotalEdgeCoverage) {
            edgesToBeRemoved.push_back(e);
            continue;
        }

        // Check coverage per strand.
        array<uint64_t, 2> perStrandCoverage = {0, 0};
        for(const auto& marker: graph[e].markers) {
            const auto sequenceId = marker.first;
            const OrientedReadId orientedReadId = orientedReadIds[sequenceId];
            const Strand strand = orientedReadId.getStrand();
            ++perStrandCoverage[strand];
        }
        if(perStrandCoverage[0]<minPerStrandEdgeCoverage or
           perStrandCoverage[1]<minPerStrandEdgeCoverage) {
            edgesToBeRemoved.push_back(e);
        }
    }

    // Remove the edges we flagged.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }

}


void MiniAssemblyMarkerGraph::fillBubbleContainment(Bubble& bubble)
{
    Graph& graph = *this;

    // The contains vector is indexed by [sequenceId][branchId].
    bubble.contains.resize(orientedReadIds.size(),
        vector<bool>(bubble.branches.size(), false));

    // Loop over all branches.
    for(uint64_t branchId=0; branchId<bubble.branches.size(); branchId++) {
        const vector<edge_descriptor>& branch = bubble.branches[branchId];

        // Loop over edges of this branch.
        for(const edge_descriptor e: branch) {
            const auto& markers = graph[e].markers;

            // Loop over markers of this edge.
            for(const auto marker: markers) {
                const uint64_t sequenceId = marker.first;
                SHASTA_ASSERT(sequenceId < bubble.contains.size());
                bubble.contains[sequenceId][branchId] = true;
            }
        }

    }


    // Now we can create the branch table.
    bubble.branchTable.resize(orientedReadIds.size(), -1);
    for(uint64_t sequenceId=0; sequenceId<orientedReadIds.size(); sequenceId++) {
        for(uint64_t branchId=0; branchId<bubble.branches.size(); branchId++) {
            if(bubble.contains[sequenceId][branchId]) {
                int64_t& value = bubble.branchTable[sequenceId];
                if(value == -1) {
                    // This is the first time.
                    value = branchId;
                } else if(value == -2) {
                    // Do nothing.
                } else {
                    // This is the second time.
                    value = -2;
                }
            }
        }
    }
}



void MiniAssemblyMarkerGraph::findBubbles()
{
    Graph& graph = *this;

    // Find linear chains.
    using Chain = std::list<edge_descriptor>;
    vector<Chain> chains;
    findLinearChains(graph, chains);

    // Store keyed by start/end vertex.
    std::map< pair<vertex_descriptor, vertex_descriptor>, vector<Chain> > chainMap;
    for(const Chain& chain: chains) {
        SHASTA_ASSERT(not chain.empty());
        const vertex_descriptor v0 = source(chain.front(), graph);
        const vertex_descriptor v1 = target(chain.back(), graph);
        chainMap[make_pair(v0, v1)].push_back(chain);
    }

    // Each non-trivial set of chains with the same start/end vertices
    // generates a bubble.
    bubbles.clear();
    for(const auto& p: chainMap) {
        const vector<Chain>& chains = p.second;
        if(chains.size() < 2) {
            continue;
        }
        bubbles.resize(bubbles.size() + 1);
        Bubble& bubble = bubbles.back();
        bubble.v0 = p.first.first;
        bubble.v1 = p.first.second;

        // Each chain with these start/end vertices generates
        // a branch in the bubble.
        for(const Chain& chain: chains) {
            bubble.branches.resize(bubble.branches.size() + 1);
            copy(chain.begin(), chain.end(), back_inserter(bubble.branches.back()));
        }
        fillBubbleContainment(bubble);
    }
}


