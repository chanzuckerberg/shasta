#include "MetaMarkerGraph.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>

#include <queue>
#include <set>



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


void MetaMarkerGraph::writeGraphviz(const string& fileName) const
{
    using Graph = MetaMarkerGraph;
    const Graph& graph = *this;

    ofstream graphOut(fileName);
    graphOut << "digraph MetaMarkerGraph {\n";
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        graphOut << graph[v0].vertexId << "->" << graph[v1].vertexId <<
            " [penwidth=" << int(0.3* double(graph[e].orientedReads.size())) << "]"
            << ";\n";
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
            graph[v].vertexId << "\t" <<
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
                graph[v0].vertexId << "\t" <<
                "+\t" <<
                graph[v1].vertexId << "\t" <<
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
    csv << "VertexId,Segment id,Marker count,Coverage,Segment id and coverage\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const MetaMarkerGraphVertex& vertex = graph[v];
        csv << vertex.vertexId << ",";
        csv << vertex.segmentId << ",";
        csv << vertex.markerCount << ",";
        csv << vertex.orientedReads.size() << ",";
        csv << vertex.segmentId << "/";
        csv << vertex.orientedReads.size() << "\n";
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
        csv << graph[v0].vertexId << ",";
        csv << graph[v1].vertexId << ",";
        csv << graph[e].orientedReads.size() << "\n";
    }

}

