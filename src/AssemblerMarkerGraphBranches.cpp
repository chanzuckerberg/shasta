#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
using namespace shasta;

#include "boost/graph/iteration_macros.hpp"


// Class used only in Assembler::analyzeMarkerGraphBranches.
namespace shasta {
    namespace analyzeMarkerGraphBranches {
        class Branch {
        public:
            using EdgeId = MarkerGraph::EdgeId;
            vector<EdgeId> edges;

            bool operator==(const Branch& that) const
            {
                return edges == that.edges;
            }
            bool operator<(const Branch& that) const
            {
                return edges < that.edges;
            }

            // The oriented read ids in each branch.
            vector< vector<OrientedReadId> > orientedReadIds;


            // Remove OrientedReadId's that appear in more than one edge.
            void removeAmbiguousOrientedReadIds()
            {
                std::map<OrientedReadId, uint64_t> m;
                for(const auto& v: orientedReadIds) {
                    for (const OrientedReadId orientedReadId: v) {
                        ++m[orientedReadId];
                    }
                }
                for(auto& v: orientedReadIds) {
                    vector<OrientedReadId> u;
                    for(const OrientedReadId orientedReadId: v) {
                        if(m[orientedReadId] == 1) {
                            u.push_back(orientedReadId);
                        }
                    }
                    v = u;
                }
            }
        };
    }
}


// Analyze branches in the marker graph.
void Assembler::analyzeMarkerGraphBranches()
{
    using namespace analyzeMarkerGraphBranches;
    using VertexId = MarkerGraph::VertexId;
    using EdgeId = MarkerGraph::EdgeId;

    // Sanity checkts.
    SHASTA_ASSERT(markers.isOpen());
    SHASTA_ASSERT(markerGraph.edges.isOpen);
    SHASTA_ASSERT(markerGraph.edgesBySource.isOpen());
    SHASTA_ASSERT(markerGraph.edgesByTarget.isOpen());
    const VertexId vertexCount = markerGraph.edgesBySource.size();
    SHASTA_ASSERT(markerGraph.edgesByTarget.size() == vertexCount);

    // Find the branches.
    vector<Branch> branches;
    for(VertexId vertexId=0; vertexId<vertexCount; vertexId++) {
        const auto outgoingEdges = markerGraph.edgesBySource[vertexId];
        if(outgoingEdges.size() > 1) {
            Branch branch;
            for(const EdgeId edgeId: outgoingEdges) {
                branch.edges.push_back(edgeId);
            }
            branches.push_back(branch);
        }
    }
    deduplicate(branches);
    cout << "Found " << branches.size() << " branches." << endl;

    // Histogram the branches by number of edges.
    vector<uint64_t> histogram;
    for(const Branch& branch: branches) {
        const uint64_t n = branch.edges.size();
        if(histogram.size() <= n) {
            histogram.resize(n+1, 0);
        }
        ++histogram[n];
    }
    for(uint64_t n=2; n<histogram.size(); n++) {
        const uint64_t count = histogram[n];
        if(count) {
            cout << "Found " << count << " branches with " << n << " edges." << endl;
        }
    }


    // Write out the branch edges.
    {
        ofstream csv("BranchEdges.csv");
        csv << "BranchId,EdgeId0,EdgeId1,EdgeId2\n";
        for(uint64_t branchId=0; branchId<branches.size(); branchId++) {
            const Branch& branch = branches[branchId];
            csv << branchId << ",";
            for(const EdgeId edgeId: branch.edges) {
                csv << edgeId << ",";
            }
            csv << "\n";
        }
    }


    // Fill in the oriented read ids of each branch.
    for(Branch& branch: branches) {
        branch.orientedReadIds.resize(branch.edges.size());
        for(uint64_t i=0; i<branch.edges.size(); i++) {
            const EdgeId edgeId = branch.edges[i];
            for(const MarkerInterval& markerInterval: markerGraph.edgeMarkerIntervals[edgeId]) {
                branch.orientedReadIds[i].push_back(markerInterval.orientedReadId);
            }
            deduplicate(branch.orientedReadIds[i]);
            branch.removeAmbiguousOrientedReadIds();
        }
    }



    // Write the oriented reads of each branch.
    {
        ofstream csv("BranchDetails.csv");
        csv << "BranchId,Edge,OrientedReadId,\n";
        for(uint64_t branchId=0; branchId<branches.size(); branchId++) {
            const Branch& branch = branches[branchId];
            for(uint64_t i=0; i<branch.orientedReadIds.size(); i++) {
                const vector<OrientedReadId>& v = branch.orientedReadIds[i];
                for(const OrientedReadId orientedReadId: v) {
                    csv << branchId << ",";
                    csv << i << ",";
                    csv << orientedReadId << ",\n";
                }
            }
        }

    }



    // Now for each pair of oriented reads count how many times they appear
    // on the same edge or on different edges of a branch.
    // The map is keyed by the pair of OrientedReadIds, with the lowest one first.
    std::map< pair<OrientedReadId, OrientedReadId>, pair<uint64_t, uint64_t> > m;
    for(uint64_t branchId=0; branchId<branches.size(); branchId++) {
        const Branch& branch = branches[branchId];
        for(uint64_t i=0; i<branch.orientedReadIds.size(); i++) {
            const vector<OrientedReadId>& vi = branch.orientedReadIds[i];
            for(uint64_t j=i; j<branch.orientedReadIds.size(); j++) {
                const vector<OrientedReadId>& vj = branch.orientedReadIds[j];

                for(const OrientedReadId oi: vi) {
                    for(const OrientedReadId oj: vj) {
                        pair<OrientedReadId, OrientedReadId> p = make_pair(oi, oj);
                        if(oi == oj) {
                            SHASTA_ASSERT(i == j);
                            continue;
                        }
                        SHASTA_ASSERT(oi != oj);
                        if(oj < oi) {
                            if(i == j) {
                                continue; // Don't count it twice!
                            } else {
                                swap(p.first, p.second);
                            }
                        }
                        /*
                        const bool debug = ((p.first==OrientedReadId(699,0)) and (p.second==OrientedReadId(700,0)));
                        if(debug) {
                            cout << "Before " << branchId << " " << i << " " << j << " " << m[p].first << " " << m[p].second << "\n";
                        }
                        */
                        if(i == j) {
                            ++m[p].first;
                        } else {
                            ++m[p].second;
                        }
                        /*
                        if(debug) {
                            cout << "After " << m[p].first << " " << m[p].second << "\n";
                        }
                        */
                    }
                }

            }
        }
    }



    // Write out the counts we found.
    {
        ofstream csv("Counts.csv");
        for(const auto& p: m) {
            csv << p.first.first << ",";
            csv << p.first.second << ",";
            csv << p.second.first << ",";
            csv << p.second.second << "\n";
        }
    }



    // Each vertex is an OrientedReadId.
    // Each edge contains the number of times the two appear
    // on the same edge or different edges of a branch.
    using Graph = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        OrientedReadId,
        pair<uint64_t, uint64_t>
        >;
    using vertex_descriptor = Graph::vertex_descriptor;
    using edge_descriptor = Graph::edge_descriptor;
    Graph graph;

    std::map<OrientedReadId, vertex_descriptor> vertexMap;

    for(const auto& p: m) {
        const OrientedReadId o0 = p.first.first;
        const OrientedReadId o1 = p.first.second;

        vertex_descriptor v0;
        auto it0 = vertexMap.find(o0);
        if(it0 == vertexMap.end()) {
            v0 = add_vertex(o0, graph);
            vertexMap.insert(make_pair(o0, v0));
        } else {
            v0 = it0->second;
        }

        vertex_descriptor v1;
        auto it1 = vertexMap.find(o1);
        if(it1 == vertexMap.end()) {
            v1 = add_vertex(o1, graph);
            vertexMap.insert(make_pair(o1, v1));
        } else {
            v1 = it1->second;
        }

        add_edge(v0, v1, p.second, graph);
    }
    cout << "The initial graph has " << num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges." << endl;



    // k-NN: for each vertex only keep the k strongest edges.
    const uint64_t k = 6;
    std::set<edge_descriptor> edgesToBeKept;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        vector< pair<edge_descriptor, uint64_t> > x;
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            const auto& p = graph[e];
            const uint64_t total = p.first + p.second;
            x.push_back(make_pair(e, total));
        }
        sort(x.begin(), x.end(), OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());
        if(x.size() > k) {
            x.resize(k);
        }
        for(const auto& p: x) {
            edgesToBeKept.insert(p.first);
        }
    }

    std::set<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, Graph) {
        if(edgesToBeKept.find(e) == edgesToBeKept.end()) {
            edgesToBeRemoved.insert(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        remove_edge(e, graph);
    }
    cout << "The k-NN graph has " << num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges." << endl;



    // Write out the k-nn graph.
    {
        ofstream out("Graph.dot");
        out << "graph G{\n";
        BGL_FORALL_VERTICES(v, graph, Graph) {
            out << "\"" << graph[v] << "\";\n";
        }
        BGL_FORALL_EDGES(e, graph, Graph) {
            const auto& p = graph[e];
            if(p.second > 2) {
                continue;
            }
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            out << "\"" << graph[v0] << "\"--\"" << graph[v1] << "\";\n";
        }
        out << "}\n";
    }
}


