// Shasta.
#include "AssemblyPathGraph.hpp"
#include "html.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <set>



AssemblyPathGraph::AssemblyPathGraph(const AssemblyGraph& assemblyGraph)
{
    AssemblyPathGraph& graph = *this;

    // Create a vertex for each assembly graph vertex.
    vector<vertex_descriptor> vertexDescriptors;
    for(AssemblyGraph::VertexId vertexId=0; vertexId<assemblyGraph.vertices.size(); vertexId++) {
        const vertex_descriptor v = add_vertex(AssemblyPathGraphVertex(vertexId), graph);
        vertexDescriptors.push_back(v);
    }

    // Fill in the reverse complemented vertex.
    for(AssemblyGraph::VertexId vertexId=0; vertexId<assemblyGraph.vertices.size(); vertexId++) {
        vertex_descriptor v = vertexDescriptors[vertexId];
        graph[v].reverseComplementVertex = vertexDescriptors[assemblyGraph.reverseComplementVertex[vertexId]];
    }

    // Sanity check.
    BGL_FORALL_VERTICES(v, graph, AssemblyPathGraph) {
        SHASTA_ASSERT(graph[graph[v].reverseComplementVertex].reverseComplementVertex == v);
    }



    // Create an edge for each assembly graph edge.
    vector<edge_descriptor> edgeDescriptors;
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[edgeId];
        const AssemblyGraph::VertexId vertexId0 = edge.source;
        const AssemblyGraph::VertexId vertexId1 = edge.target;
        const vertex_descriptor v0 = vertexDescriptors[vertexId0];
        const vertex_descriptor v1 = vertexDescriptors[vertexId1];
        edge_descriptor e;
        tie(e, ignore) = add_edge(v0, v1, AssemblyPathGraphEdge(edgeId), graph);
        edgeDescriptors.push_back(e);
    }

    // Fill in the reverse complemented edge.
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
        edge_descriptor e = edgeDescriptors[edgeId];
        graph[e].reverseComplementEdge = edgeDescriptors[assemblyGraph.reverseComplementEdge[edgeId]];
    }

    // Sanity check.
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph) {
        SHASTA_ASSERT(graph[graph[e].reverseComplementEdge].reverseComplementEdge == e);
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



Tangle& AssemblyPathGraph::getTangle(TangleId tangleId)
{
    auto it = tangles.find(tangleId);
    SHASTA_ASSERT(it != tangles.end());
    Tangle& tangle = it->second;
    SHASTA_ASSERT(tangle.tangleId == tangleId);
    return tangle;
}

// Const version.
const Tangle& AssemblyPathGraph::getTangle(TangleId tangleId) const
{
    auto it = tangles.find(tangleId);
    SHASTA_ASSERT(it != tangles.end());
    const Tangle& tangle = it->second;
    SHASTA_ASSERT(tangle.tangleId == tangleId);
    return tangle;
}



TangleId AssemblyPathGraph::getReverseComplementTangle(
    TangleId tangleId) const
{
    const AssemblyPathGraph& graph = *this;

    // Get the edge of this tangle.
    const edge_descriptor e = getTangle(tangleId).edge;
    const AssemblyPathGraphEdge& edge = graph[e];

    // Get the reverse complement edge.
    const edge_descriptor eReverseComplement = edge.reverseComplementEdge;
    const AssemblyPathGraphEdge& reverseComplementEdge = graph[eReverseComplement];

    // Return its tangle.
    const TangleId reverseComplementTangleId = reverseComplementEdge.tangle;
    SHASTA_ASSERT(reverseComplementTangleId != invalidTangleId);
    return reverseComplementTangleId;
}



void AssemblyPathGraph::detangle()
{
    AssemblyPathGraph& graph = *this;

    // Detangle iteration.
    for(int iteration=0; ; ++iteration) {

        const TangleId tangleId = findNextTangle();
        if(tangleId == invalidTangleId) {
            break;
        }
        Tangle& tangle = getTangle(tangleId);
        const TangleId reverseComplementTangleId = getReverseComplementTangle(tangleId);
        cout << "Detangle iteration " << iteration <<
            " begins, working on tangle " << tangleId <<
            " and its reverse complement tangle " <<
            reverseComplementTangleId << endl;


        // If the tangle collides with its reverse complement,
        // mark it as unsolvable, for now. This is restriction is
        // not fundamental - it just requires more coding.
        if(collidesWithReverseComplement(tangleId)) {
            tangle.unsolvable = true;
            getTangle(reverseComplementTangleId).unsolvable = true;
            cout << "Reverse complement tangles " << tangleId <<
                " and " << reverseComplementTangleId <<
                " marked unsolvable because they collide." << endl;
            continue;
        }

        // Write the graph at the beginning of this iteration.
        graph.writeGraphviz("AssemblyPathGraph-" + to_string(iteration) + ".dot");
        graph.writeHtml("AssemblyPathGraph-" + to_string(iteration) + ".html");

        // Detangle this tangle and its reverse complement.
        detangle(tangleId);
        detangle(reverseComplementTangleId);
    }
}


// Detangle a single tangle.
void AssemblyPathGraph::detangle(TangleId tangleId)
{
    cout << "Detangling tangle " << tangleId << endl;
    SHASTA_ASSERT(0);
}



// Return the next tangle to work on.
// This does a linear search, which coudl be eliminated
// with appropriated data structures if it becomes a
// performance problem.
// It currently returns the tangle with the shortest path
// on the tangle edge.
TangleId AssemblyPathGraph::findNextTangle() const
{
    const AssemblyPathGraph& graph = *this;

    TangleId bestTangleId = invalidTangleId;
    uint64_t bestTanglePathLength = std::numeric_limits<uint64_t>::max();
    for(const auto& p: tangles) {
        const Tangle& tangle = p.second;
        if(tangle.unsolvable) {
            continue;   // We gave up on this one.
        }
        const uint64_t pathLength = graph[tangle.edge].pathLength;
        if(pathLength < bestTanglePathLength) {
            bestTanglePathLength = pathLength;
            bestTangleId = tangle.tangleId;
        }
    }
    return bestTangleId;
}



// Return true if a tangle collides with its reverse complement.
bool AssemblyPathGraph::collidesWithReverseComplement(TangleId tangleId) const
{
    const AssemblyPathGraph& graph = *this;
    const Tangle& tangle = getTangle(tangleId);
    const TangleId reverseComplementTangleId = getReverseComplementTangle(tangleId);

    // If the tangle is the same as its reverse complement, we have a collision.
    // This is unusual but possible.
    if(reverseComplementTangleId == tangleId) {
        return true;
    }

    // Check the in-edges.
    for(const edge_descriptor e: tangle.inEdges) {
        if(graph[e].inTangle == reverseComplementTangleId) {
            return true;
        }
    }

    // Check the out-edges.
    for(const edge_descriptor e: tangle.outEdges) {
        if(graph[e].outTangle == reverseComplementTangleId) {
            return true;
        }
    }

    // If getting here, we did not find a collision between this tangle
    // and its reverse complement.
    return false;
}



void AssemblyPathGraph::writeHtml(const string& fileName) const
{
    ofstream html(fileName);
    writeHtml(html);
}



void AssemblyPathGraph::writeHtml(ostream& html) const
{
    writeHtmlBegin(html, "Assembly path graph");
    html << "<body>"
        "<h1>Assembly path graph</h1>";
    writeVerticesHtml(html);
    writeEdgesHtml(html);
    writeTanglesHtml(html);
    html << "</body>";
}



void AssemblyPathGraph::writeVerticesHtml(ostream& html) const
{
    const AssemblyPathGraph& graph = *this;

    html << "<h2>Vertices</h2>"
        "<p>Each vertex corresponds to a vertex of the assembly graph."
        "<p><table><tr>"
        "<th>Id"
        "<th>Id of<br>reverse<br>complement<br>vertex";

    BGL_FORALL_VERTICES(v, graph, AssemblyPathGraph) {
        const AssemblyPathGraphVertex& vertex = graph[v];
        const vertex_descriptor vReverseComplement = vertex.reverseComplementVertex;
        html <<
            "<tr id=v" << vertex.vertexId << ">" <<
            "<td class=centered>" << vertex.vertexId <<
            "<td class=centered>" << graph[vReverseComplement].vertexId;
    }

    html << "</table>";
}



void AssemblyPathGraph::writeEdgesHtml(ostream& html) const
{
    const AssemblyPathGraph& graph = *this;

    html << "<h2>Edges</h2>"
        "<p>Each edge corresponds to a path in the assembly graph vertex."
        "<p><table><tr>"
        "<th>Path"
        "<th>Path of<br>reverse<br>complement<br>edge"
        "<th>Source<br>vertex"
        "<th>Target<br>vertex"
        "<th>Path<br>length<br>(markers)"
        "<th>In-tangle"
        "<th>Tangle"
        "<th>Out-tangle";



    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const AssemblyGraph::VertexId vertexId0 = graph[v0].vertexId;
        const AssemblyGraph::VertexId vertexId1 = graph[v1].vertexId;
        const AssemblyPathGraphEdge& edge = graph[e];
        const edge_descriptor eReverseComplement = edge.reverseComplementEdge;
        html <<
            "<tr id='e" << edge << "'>"
            "<td class=centered>" << edge <<
            "<td class=centered>" << graph[eReverseComplement] <<
            "<td class=centered><a href='#v" << vertexId0 << "'>" << vertexId0 << "</a>" <<
            "<td class=centered><a href='#v" << vertexId1 << "'>" << vertexId1 << "</a>" <<
            "<td class=centered>" << edge.pathLength;

        html <<  "<td class=centered>";
        if(edge.inTangle != invalidTangleId) {
            html << "<a href='#t" << edge.inTangle << "'>" << edge.inTangle << "</a>";
        }

        html <<  "<td class=centered>";
        if(edge.tangle != invalidTangleId) {
            html << "<a href='#t" << edge.tangle << "'>" << edge.tangle << "</a>";
        }

        html <<  "<td class=centered>";
        if(edge.outTangle != invalidTangleId) {
            html << "<a href='#t" << edge.outTangle << "'>" << edge.outTangle << "</a>";
        }
    }

    html << "</table>";

}



void AssemblyPathGraph::writeTanglesHtml(ostream& html) const
{
    const AssemblyPathGraph& graph = *this;

    html << "<h2>Tangle</h2>"
        "A tangle is generated by each edge v<sub>0</sub>&rarr;v<sub>1</sub> "
        "for which the source vertex v<sub>0</sub> has in-degree greater than 1 "
        "and the target vertex v<sub>1</sub> has out-degree "
        "greater than 1."
        "<p><table><tr>"
        "<th>Id"
        "<th>In-edges"
        "<th>Tangle<br>edge"
        "<th>Out-edges"
        "<th>Tangle<br>matrix";

    for(const auto& p: tangles) {
        const Tangle& tangle = p.second;
        html <<
            "<tr id=t" << tangle.tangleId << ">"
            "<td class=centered>" << tangle.tangleId;

        // In-edges.
        html << "<td class=centered>";
        for(const edge_descriptor e: tangle.inEdges) {
            html << "<a href='#e" << graph[e] << "'>" <<
                graph[e] << "</a>" << " ";
        }

        // Tangle edge.
        html << "<td class=centered><a href='#e" << graph[tangle.edge] << "'>" <<
            graph[tangle.edge] << "</a>";

        // Out-edges.
        html << "<td class=centered>";
        for(const edge_descriptor e: tangle.outEdges) {
            html << "<a href='#e" << graph[e] << "'>" <<
                graph[e] << "</a>" << " ";
        }



        // Tangle matrix.
        html << "<td class=centered>";
        html << "<table style='margin-left:auto;margin-right:auto;'>";
        html << "<tr><td class=centered>";
        for(uint64_t j=0; j<tangle.outEdges.size(); j++) {
            const edge_descriptor e = tangle.outEdges[j];
            html << "<td class=centered>" << graph[e];
        }
        for(uint64_t i=0; i<tangle.inEdges.size(); i++) {
            const edge_descriptor e = tangle.inEdges[i];
            html << "<tr><td class=centered>" << graph[e];
            for(uint64_t j=0; j<tangle.outEdges.size(); j++) {
                const uint64_t value = tangle.matrix[i][j];
                html << "<td class=centered>";
                if(value) {
                    html << value;
                }
            }
        }
        html << "</table>";




        if(tangle.unsolvable) {
            html << "<td class=centered>Unsolvable";
        }

    }


    html << "</table>";

}
