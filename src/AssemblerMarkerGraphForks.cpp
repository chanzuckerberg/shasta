#include "Assembler.hpp"
#include "Forks.hpp"
using namespace shasta;


void Assembler::analyzeMarkerGraphForks()
{
    Forks forks(markers, markerGraph);
    const uint32_t maxDistance = 1000000;

    // while(true) {
        try {
            cout << "Enter a VertexId and direction (0=forward, 1=backward):" << endl;
            MarkerGraph::VertexId vertexId;
            uint64_t direction;
            std::cin >> vertexId >> direction;
            forks.analyze(vertexId, Forks::ForkDirection(direction), maxDistance);
        } catch(exception& e) {
            cout << e.what() << endl;
        }
    // }
}



#if 0



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

    // Sanity checks.
    SHASTA_ASSERT(markers.isOpen());
    SHASTA_ASSERT(markerGraph.edges.isOpen);
    SHASTA_ASSERT(markerGraph.edgesBySource.isOpen());
    SHASTA_ASSERT(markerGraph.edgesByTarget.isOpen());
    const VertexId vertexCount = markerGraph.edgesBySource.size();
    SHASTA_ASSERT(markerGraph.edgesByTarget.size() == vertexCount);

    // Find the branches.
    // For now only keep diploid branches.
    vector<Branch> branches;
    for(VertexId vertexId=0; vertexId<vertexCount; vertexId++) {
        const auto outgoingEdges = markerGraph.edgesBySource[vertexId];
        if(outgoingEdges.size() == 2) {
            Branch branch;
            for(const EdgeId edgeId: outgoingEdges) {
                branch.edges.push_back(edgeId);
            }
            branches.push_back(branch);
        }
        const auto incomingEdges = markerGraph.edgesByTarget[vertexId];
        if(incomingEdges.size() == 2) {
            Branch branch;
            for(const EdgeId edgeId: incomingEdges) {
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


    // Create a matrix with a row for each oriented read and a column for each branch.
    // Rows are indexed by OrientedReadId::getValue() and the matrix is stored
    // in Fortran storage scheme (that is, by columns).
    // An element is set to +1 if the oriented read appears in the first side of the branch
    // and to -1 if the oriented read appears on the second side of the branch.
    // All other matrix elements are setv to 0.
    const uint64_t readCount = getReads().readCount();
    const uint64_t orientedReadCount = 2 * readCount;
    vector<double> A(orientedReadCount * branches.size(), 0.);
    for(uint64_t branchId=0; branchId<branches.size(); branchId++) {
        const Branch& branch = branches[branchId];
        for(uint64_t i=0; i<branch.orientedReadIds.size(); i++) {
            const vector<OrientedReadId>& v = branch.orientedReadIds[i];
            for(const OrientedReadId orientedReadId: v) {
                double& matrixElement = A[branchId * orientedReadCount + uint64_t(orientedReadId.getValue())];
                if(i ==0) {
                    matrixElement = 1.;
                } else if (i==1) {
                    matrixElement = -1.;
                } else {
                    SHASTA_ASSERT(0);
                }
            }
        }
    }



    // Compute the SVD.
    const int M = int(orientedReadCount);
    const int N = int(branches.size());
    const string JOBU = "A";
    const string JOBVT = "A";
    const int LDA = M;
    vector<double> S(min(M, N));
    vector<double> U(M*M);
    const int LDU = M;
    vector<double> VT(N*N);
    const int LDVT = N;
    const int LWORK = 10 * max(M, N);
    vector<double> WORK(LWORK);
    int INFO = 0;
    cout << timestamp << "Computing the svd. Matrix is " << M << " by " << N << endl;
    dgesvd_(
        JOBU.data(), JOBVT.data(),
        M, N,
        &A[0], LDA, &S[0], &U[0], LDU, &VT[0], LDVT, &WORK[0], LWORK, INFO);
    cout << timestamp << "Done computing the svd." << endl;
    if(INFO != 0) {
        throw runtime_error("Error " + to_string(INFO) +
            " computing SVD decomposition of local read graph.");
    }



    // Write out the component of the first k left vector for each oriented read.
    {
        const uint64_t k = 40;
        ofstream csv("U.csv");
        cout << "Singular values: " << endl;
        csv << "OrientedReadId,";
        for(uint64_t i=0; i<k; i++) {
            cout << i << " " << S[i] << endl;
            csv << "U" << i << ",";
        }
        csv << "\n";
        for(ReadId readId=0; readId<readCount; readId++) {
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);
                csv << orientedReadId << ",";
                for(uint64_t i=0; i<k; i++) {
                    csv << U[i*orientedReadCount + orientedReadId.getValue()] << ",";
                }
                csv << "\n";
            }
        }
    }



#if 0
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
#endif
}


#endif
