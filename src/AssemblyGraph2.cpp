#include "AssemblyGraph2.hpp"
#include "orderPairs.hpp"
using namespace shasta;

#include <boost/graph/iteration_macros.hpp>

#include "fstream.hpp"


// The constructor creates an edge for each linear path
// in the marker graph. Therefore, immediately after construction,
// each edge has a single MarkerGraphPath (no bubbles).
AssemblyGraph2::AssemblyGraph2(const MarkerGraph& markerGraph) :
    markerGraph(markerGraph)
{
    const bool debug = false;
    ofstream debugOut;
    if(debug) {
        debugOut.open("AssemblyGraph2-Constructor.txt");
    }

    const MarkerGraph::EdgeId edgeCount = markerGraph.edges.size();
    vector<bool> wasFound(edgeCount, false);

    MarkerGraphPath nextEdges;
    MarkerGraphPath previousEdges;
    MarkerGraphPath path;
    MarkerGraphPath reverseComplementedPath;



    // Main loop over all edges of the marker graph.
    // At each iteration we find a new linear path of edges.
    for(MarkerGraph::EdgeId startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {
        if(debug) {
            const MarkerGraph::Edge& startEdge = markerGraph.edges[startEdgeId];
            debugOut << "Starting a new path at edge " << startEdgeId << " " <<
                startEdge.source << "->" << startEdge.target << "\n";
        }

        // If we already found this edge, skip it.
        // It is part of a path we already found.
        if(wasFound[startEdgeId]) {
            continue;
        }

        // Follow the path forward.
        nextEdges.clear();
        MarkerGraph::EdgeId edgeId = startEdgeId;
        bool isCircular = false;
        while(true) {
            const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
            const MarkerGraph::VertexId v1 = edge.target;
            const auto outEdges = markerGraph.edgesBySource[v1];
            if(outEdges.size() != 1) {
                break;
            }
            const auto inEdges = markerGraph.edgesByTarget[v1];
            if(inEdges.size() != 1) {
                break;
            }
            edgeId = outEdges[0];
            if(edgeId == startEdgeId) {
                isCircular = true;
                if(debug) {
                    cout << "Is circular." << endl;
                }
                break;
            }
            nextEdges.push_back(edgeId);
            SHASTA_ASSERT(not wasFound[edgeId]);
            if(debug) {
                debugOut << "Forward " << edgeId << " " <<
                    edge.source << "->" << edge.target << "\n";

            }
        }

        // Follow the path backward.
        if(!isCircular) {
            previousEdges.clear();
            edgeId = startEdgeId;
            while(true) {
                const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
                const MarkerGraph::VertexId v0 = edge.source;
                const auto outEdges = markerGraph.edgesBySource[v0];
                if(outEdges.size() != 1) {
                    break;
                }
                const auto inEdges = markerGraph.edgesByTarget[v0];
                if(inEdges.size() != 1) {
                    break;
                }
                edgeId = inEdges[0];
                previousEdges.push_back(edgeId);
                SHASTA_ASSERT(not wasFound[edgeId]);
                if(debug) {
                    debugOut << "Backward " << edgeId << " " <<
                        edge.source << "->" << edge.target << "\n";
                }
            }
        }

        // Gather the path.
        path.clear();
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter(path));
        path.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter(path));

        // Mark all the edges in the path as found.
        for(const MarkerGraph::EdgeId edgeId: path) {
            wasFound[edgeId] = true;
        }

        // Store this path as a new edge of the assembly graph.
        addEdge(path);
        if(debug) {
            for(const MarkerGraph::EdgeId edgeId: path) {
                const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
                debugOut << "Path " << edgeId << " " <<
                    edge.source << "->" << edge.target << "\n";
            }
        }

        // Also construct the reverse complemented path.
        reverseComplementedPath.clear();
        for(const MarkerGraph::EdgeId edgeId: path) {
            const MarkerGraph::EdgeId edgeIdRc = markerGraph.reverseComplementEdge[edgeId];
            reverseComplementedPath.push_back(edgeIdRc);
        }
        std::reverse(reverseComplementedPath.begin(), reverseComplementedPath.end());



        // Figure out if the reverse complemented chain is the same
        // as the original chain. This can happen in exceptional cases.
        bool isSelfComplementary = false;
        if(!isCircular) {
            isSelfComplementary = (path == reverseComplementedPath);
        } else {

            // For a circular path the test is more complex.
            // We check if the reverse complement of the first edge
            // is in the path.
            isSelfComplementary =
                find(path.begin(), path.end(), reverseComplementedPath.front()) != path.end();
        }


        // Store the reverse complemented path, if different from the original one.
        if(not isSelfComplementary) {
            addEdge(reverseComplementedPath);
            for(const MarkerGraph::EdgeId edgeIdRc: reverseComplementedPath) {
                SHASTA_ASSERT(not wasFound[edgeIdRc]);
                wasFound[edgeIdRc] = true;
            }

            if(debug) {
                for(const MarkerGraph::EdgeId edgeId: path) {
                    const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
                    debugOut << "Reverse complemented path " << edgeId << " " <<
                        edge.source << "->" << edge.target << "\n";
                }
            }
        }

    }



    // Check that all edges of the marker graph were found.;
    SHASTA_ASSERT(find(wasFound.begin(), wasFound.end(), false) == wasFound.end());

    cout << "The initial AssemblyGraph2 has " << num_vertices(*this) << " vertices and " <<
        num_edges(*this) << " edges." << endl;
    writeGfaNoSequence("AssemblyGraph2-0");

    // Gather bubble edges.
    gatherBubbles();
    cout << "After gathering bubbles, the AssemblyGraph2 has " << num_vertices(*this) << " vertices and " <<
        num_edges(*this) << " edges." << endl;
    writeGfaNoSequence("AssemblyGraph2-1");
}



// Get the vertex descriptor for the vertex corresponding to
// a given MarkerGraph::VertexId, creating the vertex if necessary.
AssemblyGraph2::vertex_descriptor AssemblyGraph2::getVertex(MarkerGraph::VertexId vertexId)
{
    const auto it = vertexMap.find(vertexId);
    if(it == vertexMap.end()) {
        const vertex_descriptor v = add_vertex(AssemblyGraph2Vertex(vertexId), *this);
        vertexMap.insert(make_pair(vertexId, v));
        return v;
    } else {
        return it->second;
    }
}



// Create a new edges corresponding to the given path.
// Also create the vertices if necessary.
void AssemblyGraph2::addEdge(const MarkerGraphPath& path)
{
    // Get the first and last edge of the path.
    const MarkerGraph::EdgeId edgeId0 = path.front();
    const MarkerGraph::EdgeId edgeId1 = path.back();
    const MarkerGraph::Edge edge0 = markerGraph.edges[edgeId0];
    const MarkerGraph::Edge edge1 = markerGraph.edges[edgeId1];

    // Get the first and last vertex of the path.
    // This creates the vertices if necessary.
    const MarkerGraph::VertexId vertexId0 = edge0.source;
    const MarkerGraph::VertexId vertexId1 = edge1.target;
    const vertex_descriptor v0 = getVertex(vertexId0);
    const vertex_descriptor v1 = getVertex(vertexId1);

    // Create the edge.
    add_edge(v0, v1, AssemblyGraph2Edge(nextEdgeId++, path), *this);
}



void AssemblyGraph2::writeCsv(const string& baseName) const
{
    writeVerticesCsv(baseName + "-Vertices.csv");
    writeEdgesCsv(baseName + "-Edges.csv");
    writeEdgeDetailsCsv(baseName + "-EdgeDetails.csv");
}



void AssemblyGraph2::writeVerticesCsv(const string& fileName) const
{
    const AssemblyGraph2& graph = *this;

    ofstream csv(fileName);
    csv << "VertexId0\n";
    BGL_FORALL_VERTICES(v, graph, AssemblyGraph2) {
        csv << graph[v].markerGraphVertexId << "\n";
    }
}



void AssemblyGraph2::writeEdgesCsv(const string& fileName) const
{
    const AssemblyGraph2& graph = *this;

    ofstream csv(fileName);
    csv << "VertexId0,VertexId1\n";
    BGL_FORALL_EDGES(e, graph, AssemblyGraph2) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        csv << graph[v0].markerGraphVertexId << ",";
        csv << graph[v1].markerGraphVertexId << "\n";
    }

}



void AssemblyGraph2::writeEdgeDetailsCsv(const string& fileName) const
{
    const AssemblyGraph2& graph = *this;

    ofstream csv(fileName);
    csv << "FirstVertexId,LastVertexId,Branch,Position,EdgeId,VertexId0,VertexId1\n";
    BGL_FORALL_EDGES(e, graph, AssemblyGraph2) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const MarkerGraph::VertexId vertexId0 = graph[v0].markerGraphVertexId;
        const MarkerGraph::VertexId vertexId1 = graph[v1].markerGraphVertexId;
        for(uint64_t i=0; i<graph[e].markerGraphPaths.size(); i++) {
            const MarkerGraphPath& path = graph[e].markerGraphPaths[i];
            for(uint64_t j=0; j<path.size(); j++) {
                const MarkerGraph::EdgeId edgeId = path[j];
                const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
                csv << vertexId0 << ",";
                csv << vertexId1 << ",";
                csv << i << ",";
                csv << j << ",";
                csv << edgeId << ",";
                csv << edge.source << ",";
                csv << edge.target << "\n";
            }
        }
    }

}



// Finds edges that form bubbles, then combine
// each of them into a single edge with multiple paths.
void AssemblyGraph2::gatherBubbles()
{
    using G = AssemblyGraph2;
    G& g = *this;

    // Vector to contain pairs (v1, e01) for a single v0.
    // Here, e01 is an edge v0->v1.
    vector< pair<vertex_descriptor, edge_descriptor> > next;

    // Ploidy histogram.
    vector<uint64_t> ploidyHistogram;


    // Look for sets of parallel edges v0->v1.
    BGL_FORALL_VERTICES(v0, g, G) {

        next.clear();
        BGL_FORALL_OUTEDGES(v0, e01, g, G) {
            const vertex_descriptor v1 = target(e01, g);
            next.push_back(make_pair(v1, e01));
        }

        // Sort them by v1.
        sort(next.begin(), next.end(), OrderPairsByFirstOnly<vertex_descriptor, edge_descriptor>());


        // Find streaks with the same v1.
        for(auto it=next.begin(); it!=next.end(); /*Incremented in the loop*/) {
            const vertex_descriptor v1 = it->first;

            auto jt = it;
            for(; jt!=next.end(); ++jt) {
                if(jt->first != v1) {
                    break;
                }
            }

            // Here, it and jt define a streak with the same v1.
            const uint64_t ploidy = jt - it;
            if(ploidy >= ploidyHistogram.size()) {
                ploidyHistogram.resize(ploidy + 1);
            }
            ++ploidyHistogram[ploidy];


            // Combine the edges in the streak into a new edge.
            if(ploidy > 1) {
                edge_descriptor eNew;
                bool edgeWasAdded = false;
                tie(eNew, edgeWasAdded) = add_edge(v0, v1, AssemblyGraph2Edge(nextEdgeId++), g);
                SHASTA_ASSERT(edgeWasAdded);
                AssemblyGraph2Edge& edgeNew = g[eNew];

                for(auto kt=it; kt!=jt; kt++) {
                    const edge_descriptor eOld = kt->second;
                    const AssemblyGraph2Edge& edgeOld = g[eOld];
                    copy(edgeOld.markerGraphPaths.begin(), edgeOld.markerGraphPaths.end(),
                       back_inserter(edgeNew.markerGraphPaths));
                    boost::remove_edge(eOld, g);
                }
            }

            // Prepare to process the next streak.
            it = jt;
        }
    }

    cout << "Ploidy histogram (counting both strands):" << endl;
    for(uint64_t ploidy=2; ploidy<ploidyHistogram.size(); ploidy++) {
        cout << "Ploidy " << ploidy << ": " << ploidyHistogram[ploidy] << " bubbles." << endl;
    }

}



void AssemblyGraph2::writeGfaNoSequence(const string& baseName) const
{
    using G = AssemblyGraph2;
    const G& g = *this;

    // Open the gfa and write the header.
    ofstream gfa(baseName + ".gfa");
    gfa << "H\tVN:Z:1.0\n";

    // Open the csv and write the header.
    ofstream csv(baseName + ".csv");
    csv << ",Color,First,Last\n";



    // Each edge of the AssemblyGraph2 generates a gfa Segment
    // for each of its marker graph paths.
    BGL_FORALL_EDGES(e, g, G) {
        const AssemblyGraph2Edge& edge = g[e];

        for(uint64_t branchId=0; branchId<edge.markerGraphPaths.size(); branchId++) {
            const MarkerGraphPath& path = edge.markerGraphPaths[branchId];

            gfa << "S\t" << edge.pathId(branchId) << "\t*\tLN:i:" << path.size() << "\n";
            csv << edge.pathId(branchId) << ",Green," <<
                path.front() << "," << path.back() << "\n";
        }
    }



    // Generate link record.
    // For each vertex, we generate a Link for each pair of
    // incoming/outgoing marker graph paths.
    BGL_FORALL_VERTICES(v, g, G) {

        // Loop over marker graph paths of incoming edges.
        BGL_FORALL_INEDGES(v, e0, g, G) {
            const AssemblyGraph2Edge& edge0 = g[e0];
            for(uint64_t i0=0; i0<edge0.markerGraphPaths.size(); i0++) {

                // Loop over marker graph paths of outgoing edges.
                BGL_FORALL_OUTEDGES(v, e1, g, G) {
                    const AssemblyGraph2Edge& edge1 = g[e1];
                    for(uint64_t i1=0; i1<edge1.markerGraphPaths.size(); i1++) {

                        gfa << "L\t" <<
                            edge0.pathId(i0) << "\t+\t" <<
                            edge1.pathId(i1) << "\t+\t\n";

                    }
                }
            }
        }
    }
}

