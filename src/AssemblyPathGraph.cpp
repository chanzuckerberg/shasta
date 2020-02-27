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


        // Write is as a pseudo vertex.
        const string pseudoVertexName =
            "\"" +
            to_string(vertex0.vertexId) +
            "to" +
            to_string(vertex1.vertexId) +
            "\"";
        s << pseudoVertexName << " [";
        s << "shape=rectangle label=\"" << edge << "\\n" << edge.pathLength << "";
        if(edge.tangle != invalidTangleId) {
            s << "\\n" << edge.tangle;
        }
        s << "\"";



        // Color.
        if(edge.tangle != invalidTangleId) {
            // Tangle edge.
            SHASTA_ASSERT(edge.inTangle == invalidTangleId);
            SHASTA_ASSERT(edge.outTangle == invalidTangleId);
            s << " style=filled fillcolor=pink";
        } else if(edge.inTangle != invalidTangleId and edge.outTangle != invalidTangleId) {
            // The edge is an in-edge of a tangle and an out-edge of another tangle.
            s << " style=filled fillcolor=purple";
        } else if(edge.inTangle != invalidTangleId) {
            // The edge has an in-tangle, so it is an out-edge of a tangle.
            s << " style=filled fillcolor=red";
        } else if(edge.outTangle != invalidTangleId) {
            // The edge has an out-tangle, so it is an in-edge of a tangle.
            s << " style=filled fillcolor=green";
        }
        s << "];\n";



        // Write the arrows to/from the pseudovertex.
        s << vertex0.vertexId << "->" << pseudoVertexName << ";\n";
        s << pseudoVertexName << "->" << vertex1.vertexId << ";\n";
    }





    s << "}\n";
}



// Initial creation of the tangles.
void AssemblyPathGraph::createTangles()
{
    AssemblyPathGraph& graph = *this;

    // Just in case, clean up.
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph) {
        graph[e].clearTangles();
    }
    tangles.clear();
    nextTangleId = 0;


    // Consider all edges.
    BGL_FORALL_EDGES(e01, graph, AssemblyPathGraph) {
        const vertex_descriptor v0 = source(e01, graph);
        const vertex_descriptor v1 = target(e01, graph);

        // If the in-degree and out-degree are not at least 2, this edge
        // does not generate a tangle.
        if(in_degree(v0, graph) <2) {
            continue;
        }
        if(out_degree(v1, graph) <2) {
            continue;
        }

        const auto inDegree = in_degree(v0, graph);
        const auto outDegree = out_degree(v1, graph);

        Tangle tangle;
        tangle.edge = e01;
        SHASTA_ASSERT(graph[e01].tangle == invalidTangleId);
        graph[e01].tangle = nextTangleId;

        // Gather the in-edges and out-edges.
        BGL_FORALL_INEDGES(v0, e, graph, AssemblyPathGraph) {
            tangle.inEdges.push_back(e);
            SHASTA_ASSERT(graph[e].outTangle == invalidTangleId);
            graph[e].outTangle = nextTangleId;
        }
        BGL_FORALL_OUTEDGES(v1, e, graph, AssemblyPathGraph) {
            tangle.outEdges.push_back(e);
            SHASTA_ASSERT(graph[e].inTangle == invalidTangleId);
            graph[e].inTangle = nextTangleId;
        }



        // Compute the tangle matrix, which contains the number of common oriented reads
        // for each pair of in-edges and out-edges.
        vector<OrientedReadId> commonOrientedReadIds;
        tangle.matrix.resize(inDegree, vector<uint64_t>(outDegree));
        for(uint64_t inEdgeIndex=0; inEdgeIndex<inDegree; inEdgeIndex++) {
            const AssemblyPathGraphEdge& inEdge = graph[tangle.inEdges[inEdgeIndex]];
            for(uint64_t outEdgeIndex=0; outEdgeIndex<outDegree; outEdgeIndex++) {
                const AssemblyPathGraphEdge& outEdge = graph[tangle.outEdges[outEdgeIndex]];
                commonOrientedReadIds.clear();
                std::set_intersection(
                    inEdge.orientedReadIds.begin(), inEdge.orientedReadIds.end(),
                    outEdge.orientedReadIds.begin(), outEdge.orientedReadIds.end(),
                    back_inserter(commonOrientedReadIds));
                tangle.matrix[inEdgeIndex][outEdgeIndex] = commonOrientedReadIds.size();
            }
        }

#if 0


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
#endif

        tangle.tangleId = nextTangleId;
        tangles.insert(make_pair(nextTangleId++, tangle));

    }
    cout << "Found " << tangles.size() << " tangles." << endl;
}




void AssemblyPathGraph::writeTangles(const string& fileName) const
{
    ofstream file(fileName);
    writeTangles(file);
}
void AssemblyPathGraph::writeTangles(ostream& file) const
{
    const AssemblyPathGraph& graph = *this;

    for(const auto& p: tangles) {
        const Tangle& tangle = p.second;
        file << "Tangle " << tangle.tangleId << endl;
        file << "Tangle edge path " << graph[tangle.edge] << endl;
        file << "Tangle edge path length " << graph[tangle.edge].pathLength << " markers." << endl;
        file << "In-degree " << tangle.inDegree() << " out-degree " << tangle.outDegree() << endl;

        file << "In-edges: ";
        for(const edge_descriptor e: tangle.inEdges) {
            file << graph[e] << " ";
        }
        file << endl;

        file << "Out-edges: ";
        for(const edge_descriptor e: tangle.outEdges) {
            file << graph[e] << " ";
        }
        file << endl;

        file << "Tangle matrix:" << endl;
        for(uint64_t i=0; i<tangle.inDegree(); i++) {
            for(uint64_t j=0; j<tangle.outDegree(); j++) {
                file << graph[tangle.inEdges[i]] << " ";
                file << graph[tangle.outEdges[j]] << " ";
                file << tangle.matrix[i][j] << endl;
            }
        }

    }
}



void AssemblyPathGraph::detangle()
{
    // AssemblyPathGraph& graph = *this;

    // Detangle iteration.
    for(int iteration=0; ; ++iteration) {



    }
}
