#include "SegmentGraph.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>

#include <list>
#include <queue>
#include <set>


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




// Find chains in the segment graph.
// This assumes that no vertex has in-degree or out-degree
// greater than one.
void SegmentGraph::findChains()
{
    SegmentGraph& graph = *this;

    chains.clear();

    // Flag all vertices as unprocessed.
    std::set<vertex_descriptor> unprocessedVertices;
    BGL_FORALL_VERTICES(v, graph, SegmentGraph) {
        SHASTA_ASSERT(in_degree(v, graph) < 2);
        SHASTA_ASSERT(out_degree(v, graph) < 2);
        unprocessedVertices.insert(v);
    }



    // Keep going as long as there are unprocessed vertices.
    while(not unprocessedVertices.empty()) {
        const vertex_descriptor v0 = *unprocessedVertices.begin();
        unprocessedVertices.erase(v0);

        std::list<vertex_descriptor> chain;
        chain.push_back(v0);

        // Follow the chain forward.
        bool isCircular = false;
        vertex_descriptor v = v0;
        while(true) {
            out_edge_iterator begin, end;
            tie(begin, end) = out_edges(v, graph);
            if(begin == end) {
                break;  // Dead end
            }
            edge_descriptor e = *begin;
            v = target(e, graph);
            if(v == v0) {
                isCircular = true;
                break;
            } else {
                chain.push_back(v);
                SHASTA_ASSERT(unprocessedVertices.find(v) != unprocessedVertices.end());
                unprocessedVertices.erase(v);
            }
        }

        // Follow the chain backward.
        if(not isCircular) {
            v = v0;
            while(true) {
                in_edge_iterator begin, end;
                tie(begin, end) = in_edges(v, graph);
                if(begin == end) {
                    break;  // Dead end
                }
                edge_descriptor e = *begin;
                v = source(e, graph);
                SHASTA_ASSERT(v != v0);
                chain.push_front(v);
                SHASTA_ASSERT(unprocessedVertices.find(v) != unprocessedVertices.end());
                unprocessedVertices.erase(v);
            }
        }

        // Store this chain.
        chains.push_back(Chain());
        Chain& storedChain = chains.back();
        storedChain.isCircular = isCircular;
        copy(chain.begin(), chain.end(), back_inserter(storedChain.vertices));
    }
    cout << "Found " << chains.size() << " chains." << endl;



    // Write the chains.
    ofstream chainsCsv("Chains.csv");
    chainsCsv << "Chain,Circular,Segment\n";
    for(uint64_t chainId=0; chainId< chains.size(); chainId++) {
        const Chain& chain = chains[chainId];
        for(uint64_t i=0; i<chain.vertices.size(); i++) {
            chainsCsv << chainId << ",";
            chainsCsv << (chain.isCircular ? "Yes" : "No") << ",";
            chainsCsv << graph[chain.vertices[i]].assemblyGraphEdgeId << "\n";
        }
    }
}
