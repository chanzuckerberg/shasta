// Shasta.
#include "AssemblyPathGraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard library.
#include "fstream.hpp"
#include <set>



AssemblyPathGraph::AssemblyPathGraph(const AssemblyGraph& assemblyGraph)
{
    AssemblyPathGraph& graph = *this;

    // Create a vertex for each assembly graph vertex.
    vector<vertex_descriptor> vertices;
    for(AssemblyGraph::VertexId vertexId=0; vertexId<assemblyGraph.vertices.size(); vertexId++) {
        const vertex_descriptor v = add_vertex(AssemblyPathGraphVertex(vertexId), graph);
        vertices.push_back(v);
    }



    // Create an edge for each assembly graph edge.
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[edgeId];
        const AssemblyGraph::VertexId vertexId0 = edge.source;
        const AssemblyGraph::VertexId vertexId1 = edge.target;
        const vertex_descriptor v0 = vertices[vertexId0];
        const vertex_descriptor v1 = vertices[vertexId1];
        edge_descriptor e;
        tie(e, ignore) = add_edge(v0, v1, AssemblyPathGraphEdge(edgeId), graph);
    }

}



void AssemblyPathGraph::writeGraphviz(const string& fileName) const
{
    ofstream file(fileName);
    writeGraphviz(file);
}



void AssemblyPathGraph::writeGraphviz(ostream& s) const
{
    const AssemblyPathGraph& graph = *this;

    s << "digraph G {\n";

    // Default attributes.
    s << "layout=sfdp;\n";
    s << "K=10;\n";
    s << "overlap=false;\n";
    s << "splines=true;\n";
    s << "smoothing=triangle;\n";
    s << "node [shape=point];\n";

    // This turns off the tooltip on the graph and the edges.
    s << "tooltip = \" \";\n";



    // Vertices.
    BGL_FORALL_VERTICES(v, graph, AssemblyPathGraph) {
        const AssemblyPathGraphVertex& vertex = graph[v];
        s << vertex.vertexId;

        s << " [";

        s << "tooltip=\"" << vertex.vertexId << "\"";

        s << "]";

        s << "\n";
    }



    // Edges. We write each edge as an additional pseudovertex.
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph) {
        const AssemblyPathGraphEdge& edge = graph[e];

        // Get the vertices.
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const AssemblyPathGraphVertex& vertex0 = graph[v0];
        const AssemblyPathGraphVertex& vertex1 = graph[v1];

        const bool isTangle = in_degree(v0, graph) > 1 and out_degree(v1, graph) > 1;

        // Write is as a pseudo vertex.
        const string pseudoVertexName =
            "\"" +
            to_string(vertex0.vertexId) +
            "to" +
            to_string(vertex1.vertexId) +
            "\"";
        s << pseudoVertexName << " [";
        s << "shape=rectangle label=\"" << edge << "\\n" << edge.pathLength << "\"";
        if(isTangle) {
            s << " style=filled fillcolor=pink";
        }
        s << "];\n";

            // Write the arrows to/from the pseudovertex.
        s << vertex0.vertexId << "->" << pseudoVertexName << ";\n";
        s << pseudoVertexName << "->" << vertex1.vertexId << ";\n";
    }





    s << "}\n";
}



void AssemblyPathGraph::detangle()
{
    AssemblyPathGraph& graph = *this;

    // Detangle iteration.
    for(int iteration=0; ; ++iteration) {

        // Find the shortest edge that can be detangled.
        edge_descriptor eTangle;
        uint64_t tangleLength = std::numeric_limits<uint64_t>::max();

        BGL_FORALL_EDGES(e01, graph, AssemblyPathGraph) {
            const vertex_descriptor v0 = source(e01, graph);
            const vertex_descriptor v1 = target(e01, graph);
            if(in_degree(v0, graph) <2) {
                continue;   // Not a tagle.
            }
            if(out_degree(v1, graph) <2) {
                continue;   // Not a tangle.
            }

            const AssemblyPathGraphEdge& edge01 = graph[e01];
            const auto inDegree = in_degree(v0, graph);
            const auto outDegree = out_degree(v1, graph);

            cout << "Tangle edge " << edge01 <<
                " path length " << graph[e01].pathLength <<
                " in-degree " << inDegree <<
                " out-degree " << outDegree << endl;

            // Gather the in-edges and out-edges.
            vector<edge_descriptor> inEdges;
            vector<edge_descriptor> outEdges;
            BGL_FORALL_INEDGES(v0, e, graph, AssemblyPathGraph) {
                inEdges.push_back(e);
            }
            BGL_FORALL_OUTEDGES(v1, e, graph, AssemblyPathGraph) {
                outEdges.push_back(e);
            }


            // Count common oriented reads for each pair of in-edges and
            // out-edges.
            vector<OrientedReadId> commonOrientedReadIds;
            vector< vector<uint64_t> > tangleMatrix(inEdges.size(), vector<uint64_t>(outEdges.size()));
            cout << "Tangle matrix:" << endl;
            for(uint64_t inEdgeIndex=0; inEdgeIndex<inEdges.size(); inEdgeIndex++) {
                const AssemblyPathGraphEdge& inEdge = graph[inEdges[inEdgeIndex]];
                for(uint64_t outEdgeIndex=0; outEdgeIndex<outEdges.size(); outEdgeIndex++) {
                    const AssemblyPathGraphEdge& outEdge = graph[outEdges[outEdgeIndex]];
                    commonOrientedReadIds.clear();
                    std::set_intersection(
                        inEdge.orientedReadIds.begin(), inEdge.orientedReadIds.end(),
                        outEdge.orientedReadIds.begin(), outEdge.orientedReadIds.end(),
                        back_inserter(commonOrientedReadIds));
                    tangleMatrix[inEdgeIndex][outEdgeIndex] = commonOrientedReadIds.size();
                    cout << inEdge << " " << outEdge << " " << commonOrientedReadIds.size() << endl;
                }
            }



            // Count the non-zero elements in each row/column of the tangle.
            vector<uint64_t> inCounts(inDegree, 0);
            vector<uint64_t> outCounts(outDegree, 0);
            for(uint64_t inEdgeIndex=0; inEdgeIndex<inEdges.size(); inEdgeIndex++) {
                for(uint64_t outEdgeIndex=0; outEdgeIndex<outEdges.size(); outEdgeIndex++) {
                    if(tangleMatrix[inEdgeIndex][outEdgeIndex]) {
                        ++inCounts[inEdgeIndex];
                        ++outCounts[outEdgeIndex];
                    }
                }
            }
            cout << "inCounts ";
            copy(inCounts.begin(), inCounts.end(), ostream_iterator<uint64_t>(cout," "));
            cout << endl;
            cout << "outCounts ";
            copy(outCounts.begin(), outCounts.end(), ostream_iterator<uint64_t>(cout," "));
            cout << endl;

            const bool canDetangle =
                std::count(inCounts.begin(), inCounts.end(), 1) == int64_t(inCounts.size()) and
                std::count(outCounts.begin(), outCounts.end(), 1) == int64_t(outCounts.size());
            if(not canDetangle) {
                cout << "This edge cannot be detangled due to ambiguity." << endl;
            }

            if(not canDetangle) {
                continue;
            }

            if(graph[e01].pathLength >= tangleLength) {
                continue;   // We already have a shorter one.
            }

            eTangle = e01;
            tangleLength = graph[e01].pathLength;


        }

        cout << "This iteration will detangle edge " << graph[eTangle] << endl;

        break;  // For now.


    }
}
