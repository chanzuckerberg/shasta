#include "SegmentGraph.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>

#include <queue>


void SegmentGraph::removeLowCoverageEdges(uint64_t minCoverage)
{
    SegmentGraph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, SegmentGraph) {
        if(graph[e].coverage < minCoverage) {
            edgesToBeRemoved.push_back(e);
        }
    }

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



// Approximate transitive reduction.
// Processes edge in order of increasing coverage.
// If for an edge v0->v1 there is a path from v0 to v1
// of distance up to maxDistance and
// that does not use that edge, remove edge v0->v1.
void SegmentGraph::transitiveReduction(uint64_t maxDistance)
{
    SegmentGraph& graph = *this;

    // Gather edges and their coverage.
    vector< pair<uint64_t, edge_descriptor> > edgeTable;
    BGL_FORALL_EDGES(e, graph, SegmentGraph) {
        edgeTable.push_back(make_pair(graph[e].coverage, e));
    }
    sort(edgeTable.begin(), edgeTable.end());



    // Process the edges in order of increasing coverage.
    for(const auto& p: edgeTable) {
        const edge_descriptor e01 = p.second;
        const vertex_descriptor v0 = source(e01, graph);
        const vertex_descriptor v1 = target(e01, graph);

        // Do a BFS starting at v0, looking for v1, and stopping at the maximum distance.
        std::map<vertex_descriptor, uint64_t> distanceMap;
        std::queue<vertex_descriptor> q;
        distanceMap.insert(make_pair(v0, 0));
        q.push(v0);
        bool found = false;
        while(not q.empty()) {
            const vertex_descriptor u0 = q.front();
            q.pop();
            const uint64_t distance0 = distanceMap[u0];
            const uint64_t distance1 = distance0 + 1;
            BGL_FORALL_OUTEDGES(u0, f01, graph, SegmentGraph) {
                if(f01 == e01) {
                    continue;
                }
                const vertex_descriptor u1 = target(f01, graph);
                if(u1 == v1) {
                    boost::remove_edge(e01, graph);
                    found = true;
                    break;
                }
                if(distanceMap.find(u1) != distanceMap.end()) {
                    continue;
                }
                distanceMap.insert(make_pair(u1, distance1));
                if(distance1 < maxDistance) {
                    q.push(u1);
                }
            }
            if(found) {
                break;
            }
        }

    }


}


void SegmentGraph::writeGraphviz(const string& fileName) const
{
    const SegmentGraph& graph = *this;

    ofstream graphOut(fileName);
    graphOut << "digraph SegmentGraph {\n";

    // Write the vertices.
    BGL_FORALL_VERTICES(v, graph, SegmentGraph) {
        graphOut << graph[v].assemblyGraphEdgeId << ";\n";
    }

    // Write the edges.
    BGL_FORALL_EDGES(e, graph, SegmentGraph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        graphOut << graph[v0].assemblyGraphEdgeId << "->";
        graphOut << graph[v1].assemblyGraphEdgeId;
        graphOut << " [penwidth=\"" << 0.1*double(graph[e].coverage) << "\"];\n";
    }

    graphOut << "}\n";
}

