#include "MetaMarkerGraph.hpp"
#include "chokePoints.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>

#include <list>
#include <queue>
#include <set>
#include <stack>


void MetaMarkerGraph::addVertex(
    SegmentId segmentId,
    uint64_t markerCount,
    const vector< pair<OrientedReadId, uint64_t> >& orientedReads)
{
    add_vertex(
        MetaMarkerGraphVertex(segmentId, markerCount, orientedReads),
        *this);
}



// Generate the sequenceNumber field for each vertex.
void MetaMarkerGraph::generateSequenceNumbers()
{
    MetaMarkerGraph& graph = *this;
    vector<uint64_t> vertexCountBySegmentId;

    BGL_FORALL_VERTICES(v, graph, MetaMarkerGraph) {
        MetaMarkerGraphVertex& vertex = graph[v];
        const SegmentId segmentId = vertex.segmentId;

        if(segmentId >= vertexCountBySegmentId.size()) {
            vertexCountBySegmentId.resize(segmentId + 1, 0);
        }

        graph[v].sequenceNumber = vertexCountBySegmentId[segmentId];
        ++vertexCountBySegmentId[segmentId];
    }

}


void MetaMarkerGraph::createEdges()
{
    using Graph = MetaMarkerGraph;
    Graph& graph = *this;



    // Construct the sequence of vertices encountered by each oriented read.
    std::map<OrientedReadId, vector<vertex_descriptor> > pseudoPaths;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const vector< pair<OrientedReadId, uint64_t> >& orientedReads =
            graph[v].orientedReads;

        // Loop over oriented reads of this vertex and the corresponding
        // metaOrdinals.
        for(const auto& p: orientedReads) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t metaOrdinal = p.second;

            // Access the pseudo-path for this oriented read.
            vector<vertex_descriptor>& pseudoPath = pseudoPaths[orientedReadId];

            // Make sure we have a slot for this metaOrdinal.
            if(pseudoPath.size() <= metaOrdinal) {
                pseudoPath.resize(metaOrdinal+1, null_vertex());
            }

            // Now we can store it.
            pseudoPath[metaOrdinal] = v;
        }
    }


    // Check that the pseudo-path don't have any missing vertices.
    for(const auto& p: pseudoPaths) {
        const vector<vertex_descriptor>& pseudoPath = p.second;
        for(const vertex_descriptor v: pseudoPath) {
            SHASTA_ASSERT(v != null_vertex());
        }
    }




    // Now we can create the edges by looping over the pseudo-path
    // of each orientedRead.
    for(const auto& p: pseudoPaths) {
        const OrientedReadId orientedReadId = p.first;
        const vector<vertex_descriptor>& pseudoPath = p.second;

        // Loop over successive vertices in the pseudo-path of this oriented read.
        for(uint64_t metaOrdinal1=1; metaOrdinal1<pseudoPath.size(); metaOrdinal1++) {
            const uint64_t metaOrdinal0 = metaOrdinal1 - 1;
            const vertex_descriptor v0 = pseudoPath[metaOrdinal0];
            const vertex_descriptor v1 = pseudoPath[metaOrdinal1];

            // Access the edge between these two vertices, creating it if necessary.
            bool edgeExists = false;;
            edge_descriptor e;
            tie(e, edgeExists) = edge(v0, v1, graph);
            if(not edgeExists) {
                tie(e, edgeExists) = add_edge(v0, v1, graph);
            }

            // Store this metaOrdinal in the edge.
            graph[e].orientedReads.push_back(make_pair(orientedReadId, metaOrdinal0));
        }
    }
}


// This does transitive reduction considering paths of arbitrary length.
// We directly test for reachability.
// The resulting transitive reduction is unique only if the graph
// is acyclic. In the general case, multiple transitive reductions
// are possible. For this reason, we process edges in order
// of increasing coverage.

// The boost graph library as an undocumented transitive reduction
// which requires an acyclic graph (because it uses topological_sort)
// and also does not compile for this graph type.
// The transitive closure algorithm is documented and general, but
// does not help here.
void MetaMarkerGraph::transitiveReduction()
{
    MetaMarkerGraph& graph = *this;

    // Sort edges by increasing coverage.
    vector< pair<uint64_t, edge_descriptor> > sortedEdges;
    BGL_FORALL_EDGES(e, graph, MetaMarkerGraph) {
        sortedEdges.push_back(make_pair(graph[e].orientedReads.size(), e));
    }
    sort(sortedEdges.begin(), sortedEdges.end());



    // Process all edges in order of increasing coverage.
    for(const auto& p: sortedEdges) {
        const edge_descriptor e01 = p.second;
        const vertex_descriptor v0 = source(e01, graph);
        const vertex_descriptor v1 = target(e01, graph);

        // Do a BFS starting at v0, without using this edge,
        // and stopping if we encounter v1.
        std::queue<vertex_descriptor> q;
        q.push(v0);
        std::set<vertex_descriptor> verticesFound;
        verticesFound.insert(v0);
        bool done = false;
        while(not q.empty()) {
            const vertex_descriptor u0 = q.front();
            q.pop();
            BGL_FORALL_OUTEDGES(u0, f01, graph, MetaMarkerGraph) {
                if(f01 == e01) {
                    continue;
                }
                const vertex_descriptor u1 = target(f01, graph);
                if(u1 == v1) {
                    boost::remove_edge(e01, graph);
                    done = true;
                    break;
                }
                if(verticesFound.find(u1) == verticesFound.end()) {
                    q.push(u1);
                    verticesFound.insert(u1);
                }
            }
            if(done) {
                break;
            }
        }
    }
}



// Recursively prune all leafs with coverage less than minCoverage.
void MetaMarkerGraph::prune(uint64_t minCoverage)
{
    MetaMarkerGraph& graph = *this;
    const bool debug = false;

    // Gather all leafs with low coverage.
    std::stack<vertex_descriptor> unprocessedLowCoverageLeafs;
    BGL_FORALL_VERTICES(v, graph, MetaMarkerGraph) {
        if(graph[v].coverage() >= minCoverage) {
            continue;
        }
        const bool isLeaf = (out_degree(v, graph)==0) or (in_degree(v, graph) == 0);
        if(not isLeaf) {
            continue;
        }
        unprocessedLowCoverageLeafs.push(v);
    }



    // Main recursive pruning loop.
    while(not unprocessedLowCoverageLeafs.empty()) {
        const vertex_descriptor v0 = unprocessedLowCoverageLeafs.top();
        if(debug) {
            cout << "Found leaf " << v0 << " " << graph[v0].gfaId() << endl;
        }
        unprocessedLowCoverageLeafs.pop();
        SHASTA_ASSERT(graph[v0].coverage() < minCoverage);
        SHASTA_ASSERT((out_degree(v0, graph)==0) or (in_degree(v0, graph) == 0));

        // Forward leaf. We have to check if any parents now become low coverage leafs.
        if((out_degree(v0, graph)==0)) {
            BGL_FORALL_INEDGES(v0, e, graph, MetaMarkerGraph) {
                const vertex_descriptor v1 = source(e, graph);
                if(graph[v1].coverage() >= minCoverage) {
                    continue;   // Has high coverage. Skip.
                }
                if(out_degree(v1, graph) > 1) {
                    continue;   // v1 has other children. Skip.
                }
                if(in_degree(v1, graph) == 0) {
                    continue;   // v1 is already a leaf. No reason to add to the stack. Skip.
                }
                unprocessedLowCoverageLeafs.push(v1);
                if(debug) {
                    cout << "Added to stack " << v1 << " " << graph[v1].gfaId() << endl;
                }
            }
        }

        // Backward leaf. We have to check if any children now become low coverage leafs.
        else if((in_degree(v0, graph)==0)) {
            BGL_FORALL_OUTEDGES(v0, e, graph, MetaMarkerGraph) {
                const vertex_descriptor v1 = target(e, graph);
                if(graph[v1].coverage() >= minCoverage) {
                    continue;   // Has high coverage. Skip.
                }
                if(in_degree(v1, graph) > 1) {
                    continue;   // v1 has other parents. Skip.
                }
                if(out_degree(v1, graph) == 0) {
                    continue;   // v1 is already a leaf. No reason to add to the stack. Skip.
                }
                unprocessedLowCoverageLeafs.push(v1);
                if(debug) {
                    cout << "Added to stack " << v1 << " " << graph[v1].gfaId() << endl;
                }
            }
        }

        else {
            SHASTA_ASSERT(0);
        }

        if(debug) {
            cout << "Removed " << v0 << " "  << graph[v0].gfaId() << endl;
        }
        clear_vertex(v0, graph);
        remove_vertex(v0, graph);
    }

}



// Find the vertex corresponding to a given segment id,
// or null_vertex if not exactly one found.
MetaMarkerGraph::vertex_descriptor MetaMarkerGraph::findVertex(SegmentId segmentId) const
{
    const MetaMarkerGraph& graph = *this;
    vertex_descriptor v = null_vertex();

    BGL_FORALL_VERTICES(u, graph, MetaMarkerGraph) {
        if(graph[u].segmentId == segmentId) {
            if(v == null_vertex()) {
                v = u;
            } else {
                // More than one found.
                return null_vertex();
            }
        }
    }
    return v;
}


// Construct the linear chain (path) that includes a given segment.
void MetaMarkerGraph::findLinearChain(
    SegmentId segmentId0,
    vector<SegmentId>& chain) const
{
    const MetaMarkerGraph& graph = *this;
    chain.clear();

    // Locate the one and only vertex corresponding to the start segment.
    const vertex_descriptor v0 = findVertex(segmentId0);
    if(v0 == null_vertex()) {
        return;
    }


    // If this is a branching vertex, give up and return an empty chain.
    if(in_degree(v0, graph)!=1 or out_degree(v0, graph)!=1) {
        return;
    }


    // Construct the forward portion of the chain.
    std::list<SegmentId> chainList;
    chainList.push_back(segmentId0);
    vertex_descriptor v = v0;
    while(true) {
        out_edge_iterator it;
        if(out_degree(v, graph) == 0) {
            break;
        }
        tie(it, ignore)= out_edges(v, graph);
        const edge_descriptor e = *it;
        v = target(e, graph);
        if((in_degree(v, graph)<2) and (out_degree(v, graph)<2)) {
            chainList.push_back(graph[v].segmentId);
            cout << "Forward " << graph[v].segmentId  << " " <<
                in_degree(v, graph) << " " << out_degree(v, graph)<< endl;
        } else {
            break;
        }
    }


    // Construct the backward portion of the chain.
    v = v0;
    while(true) {
        in_edge_iterator it;
        if(in_degree(v, graph) == 0) {
            break;
        }
        tie(it, ignore)= in_edges(v, graph);
        const edge_descriptor e = *it;
        v = source(e, graph);
        if((in_degree(v, graph)<2) and (out_degree(v, graph)<2)) {
            chainList.push_front(graph[v].segmentId);
            cout << "Backward " << graph[v].segmentId << endl;
        } else {
            break;
        }
    }

    chain.clear();
    copy(chainList.begin(), chainList.end(), back_inserter(chain));
}



void MetaMarkerGraph::findForwardChokePoints(
    SegmentId startSegmentId,
    vector<vertex_descriptor>& chokePoints)
{
    // To restore this, the second template argument must be vecS.
    SHASTA_ASSERT(0);

#if 0
    const MetaMarkerGraph& graph = *this;
    chokePoints.clear();

    // Locate the one and only vertex corresponding to the start segment.
    const vertex_descriptor startVertex = findVertex(startSegmentId);
    if(startVertex == null_vertex()) {
        return;
    }

    // Find the forward choke points.
    shasta::findForwardChokePoints(graph, startVertex, chokePoints);
#endif
}




void MetaMarkerGraph::findBackwardChokePoints(
    SegmentId startSegmentId,
    vector<vertex_descriptor>& chokePoints)
{
    // To restore this, the second template argument must be vecS.
    SHASTA_ASSERT(0);

#if 0
    const MetaMarkerGraph& graph = *this;
    chokePoints.clear();

    // Locate the one and only vertex corresponding to the start segment.
    const vertex_descriptor startVertex = findVertex(startSegmentId);
    if(startVertex == null_vertex()) {
        return;
    }

    // Find the backward choke points.
    shasta::findBackwardChokePoints(graph, startVertex, chokePoints);
#endif
}



void MetaMarkerGraph::writeGraphviz(
    const string& fileName,
    AssemblyGraph::EdgeId startSegmentId) const
{
    using Graph = MetaMarkerGraph;
    const Graph& graph = *this;

    ofstream graphOut(fileName);
    graphOut << "digraph MetaMarkerGraph {\n";


    BGL_FORALL_VERTICES(v, graph, Graph) {
        const MetaMarkerGraphVertex& vertex = graph[v];

        graphOut <<
            vertex.gfaId() <<
            " [label=\"" << vertex.segmentId << "\" tooltip=\"Segment " <<
            vertex.segmentId << ", " <<
            vertex.markerCount << " markers, coverage " <<
            vertex.orientedReads.size() << "\"";

        if(startSegmentId!=std::numeric_limits<SegmentId>::max() and
            vertex.segmentId == startSegmentId) {
            graphOut << " style=filled fillcolor=pink";
        } else {
            const string color = vertex.color();
            graphOut << " style=filled color=" << color << " fillcolor=" << color;
        }

        graphOut << "];\n";
    }



    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const auto coverage = graph[e].orientedReads.size();
        const auto color = graph[e].color();
        graphOut << graph[v0].gfaId() << "->" << graph[v1].gfaId() <<
            " ["
            "tooltip=\"Coverage " << coverage << "\""
            // " color=" << color << "]"
            "];\n";
    }

    graphOut << "}\n";
}


void MetaMarkerGraph::writeGfa(const string& fileName) const
{
    using Graph = MetaMarkerGraph;
    const Graph& graph = *this;

    ofstream gfa(fileName);

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";


    // Write a segment record for each vertex.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        gfa <<
            "S\t" <<
            graph[v].gfaId() << "\t" <<
            "*\t" <<
            "LN:i:" << graph[v].markerCount <<
            "\n";
    }

    // Write link records.
    BGL_FORALL_VERTICES(v0, graph, Graph) {
        BGL_FORALL_OUTEDGES(v0, e01, graph, Graph) {
            const vertex_descriptor v1 = target(e01, graph);
            gfa <<
                "L\t" <<
                graph[v0].gfaId() << "\t" <<
                "+\t" <<
                graph[v1].gfaId() << "\t" <<
                "+\t" <<
                "*\n";
        }
    }
}



void MetaMarkerGraph::writeVerticesCsv(const string& fileName) const
{
    using Graph = MetaMarkerGraph;
    const Graph& graph = *this;

    ofstream csv(fileName);
    csv << "GfaId,Segment id,Sequence number,Marker count,Coverage,Color\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const MetaMarkerGraphVertex& vertex = graph[v];
        csv << vertex.gfaId() << ",";
        csv << vertex.segmentId << ",";
        csv << vertex.sequenceNumber << ",";
        csv << vertex.markerCount << ",";
        csv << vertex.orientedReads.size() << ",";
        csv << vertex.color() << "\n";
    }

}



void MetaMarkerGraph::writeVerticesDetailCsv(const string& fileName) const
{
    using Graph = MetaMarkerGraph;
    const Graph& graph = *this;

    ofstream csv(fileName);
    csv << "GfaId,Segment id,Sequence number,Oriented read id,Position in pseudo-path\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const MetaMarkerGraphVertex& vertex = graph[v];
        for(const auto& p: vertex.orientedReads) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t positionInPseudoPath = p.second;
            csv << vertex.gfaId() << ",";
            csv << vertex.segmentId << ",";
            csv << vertex.sequenceNumber << ",";
            csv << orientedReadId << ",";
            csv << positionInPseudoPath << "\n";
        }
    }

}



void MetaMarkerGraph::writeEdgesCsv(const string& fileName) const
{
    using Graph = MetaMarkerGraph;
    const Graph& graph = *this;

    ofstream csv(fileName);
    csv << "VertexId0,VertexId1,Coverage\n";

    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        csv << graph[v0].gfaId() << ",";
        csv << graph[v1].gfaId() << ",";
        csv << graph[e].orientedReads.size() << "\n";
    }

}

