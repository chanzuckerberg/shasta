// Shasta.
#include "Assembler.hpp"
#include "LocalAssemblyGraph.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "chrono.hpp"
#include <queue>
#include <unordered_map>
#include "iterator.hpp"



// In the assembly graph, each vertex corresponds to a linear chain
// of edges in the pruned spanning subgraph of the marker graph.
// This code finds all linear chains of edges in
// the pruned spanning subgraph of the marker graph,
// and generates a vertex of the assembly graph
// for each chain it finds.
void Assembler::createAssemblyGraphVertices()
{
    // Some shorthands.
    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::EdgeId;

    // Check that we have what we need.
    checkMarkerGraphVerticesAreAvailable();
    checkMarkerGraphConnectivityIsOpen();
    const auto& edges = markerGraphConnectivity.edges;

    // Flag to control debug output.
    const bool debug = false;

    // Vector used to keep track of edges that were already found.
    const EdgeId edgeCount = markerGraphConnectivity.edges.size();
    MemoryMapped::Vector<bool> wasFound;
    wasFound.createNew(
        largeDataName("tmp-createAssemblyGraphVertices-wasFound"),
        largeDataPageSize);
    wasFound.resize(edgeCount);
    fill(wasFound.begin(), wasFound.end(), false);

    // Initialize the vertices of the assembly graph.
    assemblyGraph.vertices.createNew(
        largeDataName("AssemblyGraphVertices"),
        largeDataPageSize);



    // Work vectors reused for each chain.
    vector<EdgeId> nextEdges;
    vector<EdgeId> previousEdges;
    vector<EdgeId> chain;



    // Main loop over all edges of the pruned spanning subgraph
    // of the marker graph.
    // At each iteration we find a new linear chain of edges.
    for(EdgeId startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {
        const auto& startEdge = edges[startEdgeId];

        if(debug) {
            cout << "Working on start edge " << startEdgeId;
            cout << " " << startEdge.source << "->" << startEdge.target << endl;
        }

        // If this edge is not part of the pruned spanning subgraph, skip it.
        if(!startEdge.isInSpanningSubgraph) {
            if(debug) {
                cout << "Edge is not in the spanning subgraph." << endl;
            }
            continue;
        }
        if(startEdge.wasPruned) {
            if(debug) {
                cout << "Edge was pruned." << endl;
            }
            continue;
        }

        // If we already found this edge, skip it.
        // It is part of a chain we aready found.
        if(wasFound[startEdgeId]) {
            if(debug) {
                cout << "This edge is part of a chain we already found." << endl;
            }
            continue;
        }

        // Follow the chain forward.
        EdgeId edgeId = startEdgeId;
        while(true) {
            edgeId = nextEdgeInMarkerGraphPrunedSpanningSubgraphChain(edgeId);
            if(edgeId == invalidGlobalMarkerGraphEdgeId) {
                break;
            }
            nextEdges.push_back(edgeId);
        }

        // Follow the chain backward.
        edgeId = startEdgeId;
        while(true) {
            edgeId = previousEdgeInMarkerGraphPrunedSpanningSubgraphChain(edgeId);
            if(edgeId == invalidGlobalMarkerGraphEdgeId) {
                break;
            }
            previousEdges.push_back(edgeId);
        }

        // Gather the chain.
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter< vector<EdgeId> >(chain));
        chain.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter< vector<EdgeId> >(chain));

        // Mark all the edges in the chain as found.
        for(const EdgeId edgeId: chain) {
            wasFound[edgeId] = true;
        }

        // Store this chain as a new vertex of the assembly graph.
        assemblyGraph.vertices.appendVector(chain);

        // Write out the chain.
        if(debug) {
            cout << "Chain: ";
            cout << edges[chain.front()].source << " ";
            for(const EdgeId edgeId: chain) {
                cout << edges[edgeId].target << " ";
            }
            cout << endl;
        }

        // Cleanup.
        nextEdges.clear();
        previousEdges.clear();
        chain.clear();
    }



    // Check that all edges of the pruned spanning subgraph of the marker graph
    // were found.
    for(EdgeId edgeId=0; edgeId<edgeCount; edgeId++) {
        const auto& edge = markerGraphConnectivity.edges[edgeId];
        if(!edge.isInSpanningSubgraph) {
            CZI_ASSERT(!wasFound[edgeId]);
            continue;
        }
        if(edge.wasPruned) {
            CZI_ASSERT(!wasFound[edgeId]);
            continue;
        }
        CZI_ASSERT(wasFound[edgeId]);
    }

    wasFound.remove();

    cout << "The assembly graph has " << assemblyGraph.vertices.size() << " vertices." << endl;



    // Create a histogram of size (chain length) of assembly graph vertices.
    vector<size_t> histogram;
    for(VertexId vertexId=0; vertexId<assemblyGraph.vertices.size(); vertexId++) {
        const size_t size = assemblyGraph.vertices.size(vertexId);
        if(histogram.size() <= size) {
            histogram.resize(size+1);
        }
        ++(histogram[size]);
    }
    ofstream csv("AssemblyGraphChainLengthHistogram.csv");
    csv << "ChainLength, Frequency\n";
    for(size_t size=0; size<histogram.size(); size++) {
        const size_t frequency = histogram[size];
        if(frequency) {
            csv << size << "," << frequency << "\n";
        }
    }
}



void Assembler::accessAssemblyGraphVertices()
{
    assemblyGraph.vertices.accessExistingReadOnly(
        largeDataName("AssemblyGraphVertices"));
}



void Assembler::createAssemblyGraphEdges()
{
    cout << timestamp << "Creating assembly graph edges." << endl;

    // Check that we have what we need.
    auto& vertices = assemblyGraph.vertices;
    CZI_ASSERT(vertices.isOpen());
    checkMarkerGraphConnectivityIsOpen();

    // Shorthands for vertex and edge ids.
    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::EdgeId;

    // In the code below:
    // v indicates an assembly graph vertex id.
    // u indicates a marker graph vertex id.

    // Create a hash map with:
    // Key: Vertex id of the initial marker graph vertex
    // of the chain corresponding to an assembly graph vertex.
    // Value: the ids of the assembly graph vertices that begin there..
    // In the loop:
    // v: vertex in the assembly graph.
    // u: vertex in the marker graph.
    std::unordered_map<VertexId, vector<VertexId> > vertexMap;
    for(VertexId v=0; v<vertices.size(); v++) {
        const auto chain = vertices[v];
        CZI_ASSERT(chain.size() > 0);
        const EdgeId firstChainEdgeId = *(chain.begin());
        const MarkerGraphConnectivity::Edge& firstChainEdge =
            markerGraphConnectivity.edges[firstChainEdgeId];
        const VertexId u = firstChainEdge.source;
        vertexMap[u].push_back(v);
    }

    // Initialize the assembly graph edges.
    assemblyGraph.edges.createNew(
        largeDataName("AssemblyGraphEdges"),
        largeDataPageSize);

    // Now for each vertex in the assembly graph look up in the map
    // the last vertex of its chain.
    for(VertexId v0=0; v0<vertices.size(); v0++) {
        const auto chain = vertices[v0];
        CZI_ASSERT(chain.size() > 0);
        const EdgeId lastChainEdgeId = chain[chain.size()-1];
        const MarkerGraphConnectivity::Edge& lastChainEdge =
            markerGraphConnectivity.edges[lastChainEdgeId];
        const VertexId u = lastChainEdge.target;

        // Looking up the map gives the children of this assembly graph vertex.
        const vector<VertexId>& children = vertexMap[u];
        for(const VertexId v1: children) {
            assemblyGraph.edges.push_back(AssemblyGraph::Edge(v0, v1));
        }
    }
    cout << timestamp << "The assembly graph has " << vertices.size() << " vertices and ";
    cout << assemblyGraph.edges.size() << " edges." << endl;


    cout << timestamp << "Creating assembly graph edge by source and by target." << endl;

    assemblyGraph.edgesBySource.createNew(
        largeDataName("AssemblyGraphEdgesBySource"),
        largeDataPageSize);
    assemblyGraph.edgesByTarget.createNew(
        largeDataName("AssemblyGraphEdgesByTarget"),
        largeDataPageSize);
    assemblyGraph.edgesBySource.beginPass1(vertices.size());
    assemblyGraph.edgesByTarget.beginPass1(vertices.size());
    for(const AssemblyGraph::Edge& edge: assemblyGraph.edges) {
        assemblyGraph.edgesBySource.incrementCount(edge.source);
        assemblyGraph.edgesByTarget.incrementCount(edge.target);
    }
    assemblyGraph.edgesBySource.beginPass2();
    assemblyGraph.edgesByTarget.beginPass2();
    for(EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[edgeId];
        assemblyGraph.edgesBySource.store(edge.source, edgeId);
        assemblyGraph.edgesByTarget.store(edge.target, edgeId);
    }
    assemblyGraph.edgesBySource.endPass2();
    assemblyGraph.edgesByTarget.endPass2();

    cout << timestamp << "Done creating assembly graph edges." << endl;
}



void Assembler::accessAssemblyGraphEdges()
{
    assemblyGraph.edges.accessExistingReadOnly(
        largeDataName("AssemblyGraphEdges"));
    assemblyGraph.edgesBySource.accessExistingReadOnly(
        largeDataName("AssemblyGraphEdgesBySource"));
    assemblyGraph.edgesByTarget.accessExistingReadOnly(
        largeDataName("AssemblyGraphEdgesByTarget"));
}



// Extract a local assembly graph from the global assembly graph.
// This returns false if the timeout was exceeded.
bool Assembler::extractLocalAssemblyGraph(
    AssemblyGraph::VertexId startVertexId,
    int distance,
    double timeout,
    LocalAssemblyGraph& graph) const
{
    using vertex_descriptor = LocalAssemblyGraph::vertex_descriptor;
    using edge_descriptor = LocalAssemblyGraph::edge_descriptor;
    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::EdgeId;

    const bool debug = false;
    if(debug) {
        cout << "Begin extractLocalAssemblyGraph for vertex "
            << startVertexId << " distance " << distance << endl;
    }

    const auto startTime = steady_clock::now();

    // Add the start vertex.
    const vertex_descriptor vStart = graph.addVertex(startVertexId, 0);

    // Do the BFS.
    std::queue<vertex_descriptor> q;
    if(distance > 0) {
        q.push(vStart);
    }
    while(!q.empty()) {

        // See if we exceeded the timeout.
        if(timeout>0. && seconds(steady_clock::now() - startTime) > timeout) {
            graph.clear();
            return false;
        }

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const LocalAssemblyGraphVertex& vertex0 = graph[v0];
        const VertexId vertexId0 = vertex0.vertexId;
        const int distance0 = vertex0.distance;
        const int distance1 = distance0 + 1;

        if(debug) {
            cout << "Dequeued " << vertexId0 << " at distance " << distance0 << endl;
        }



        // Loop over children.
        const auto childEdges = assemblyGraph.edgesBySource[vertexId0];
        for(const EdgeId edgeId: childEdges) {
            const AssemblyGraph::Edge& globalEdge = assemblyGraph.edges[edgeId];
            const VertexId vertexId1 = globalEdge.target;

            if(debug) {
                cout << "Found child " << vertexId1 << endl;
            }

            // Find the vertex corresponding to this child, creating it if necessary.
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);
            if(!vertexExists) {
                v1 = graph.addVertex(vertexId1, distance1);
                if(distance1 < distance) {
                    q.push(v1);
                }
                if(debug) {
                    cout << "Vertex added " << vertexId1 << endl;
                }
            }

            // Create the edge v0->v1, if it does not already exist.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            if(!edgeExists) {
                tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
                CZI_ASSERT(edgeExists);
                if(debug) {
                    cout << "Edge added " << vertexId0 << "->" << vertexId1 << endl;
                }
            }
        }



        // Loop over parents.
        const auto parentEdges = assemblyGraph.edgesByTarget[vertexId0];
        for(const EdgeId edgeId: parentEdges) {
            const AssemblyGraph::Edge& globalEdge = assemblyGraph.edges[edgeId];
            const VertexId vertexId1 = globalEdge.source;

            if(debug) {
                cout << "Found parent " << vertexId1 << endl;
            }

            // Find the vertex corresponding to this child, creating it if necessary.
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);
            if(!vertexExists) {
                v1 = graph.addVertex(vertexId1, distance1);
                if(distance1 < distance) {
                    q.push(v1);
                }
                if(debug) {
                    cout << "Vertex added " << vertexId1 << endl;
                }
            }

            // Create the edge v1->v0, if it does not already exist.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v1, v0, graph);
            if(!edgeExists) {
                tie(e, edgeExists) = boost::add_edge(v1, v0, graph);
                CZI_ASSERT(edgeExists);
                if(debug) {
                    cout << "Edge added " << vertexId1 << "->" << vertexId0 << endl;
                }
            }
        }
    }



    // The BFS process did not create edges between vertices at maximum distance.
    // Do it now.
    // Loop over all vertices at maximum distance.
    BGL_FORALL_VERTICES(v0, graph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex0 = graph[v0];
        if(vertex0.distance != distance) {
            continue;
        }
        const AssemblyGraph::VertexId vertexId0 = vertex0.vertexId;

        // Loop over the children that exist in the local assembly graph
        // and are also at maximum distance.
        const auto childEdges = assemblyGraph.edgesBySource[vertexId0];
        for(uint64_t edgeId: childEdges) {
            const auto& edge = assemblyGraph.edges[edgeId];

            const AssemblyGraph::VertexId vertexId1 = edge.target;
            CZI_ASSERT(edge.source == vertexId0);
            CZI_ASSERT(vertexId1 < assemblyGraph.vertices.size());

            // See if we have a vertex for this global vertex id.
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);

            // If it does not exist in the local marker graph, skip.
            if(!vertexExists) {
                continue;
            }

            // If it is not at maximum distance, skip.
            const LocalAssemblyGraphVertex& vertex1 = graph[v1];
            if(vertex1.distance != distance) {
                continue;
            }

            // There is no way we already created this edge.
            // Check that this is the case.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            CZI_ASSERT(!edgeExists);

            // Add the edge.
            tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
            CZI_ASSERT(edgeExists);

        }
    }

    if(debug) {
        cout << "Vertices:" << endl;
        BGL_FORALL_VERTICES(v, graph, LocalAssemblyGraph) {
            cout << graph[v].vertexId << endl;
        }
        cout << "Edges:" << endl;
        BGL_FORALL_EDGES(e, graph, LocalAssemblyGraph) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            cout << graph[v0].vertexId << "->";
            cout << graph[v1].vertexId << endl;
        }

    }


    return true;
}



// Assemble sequence for a vertex of the assembly graph.
// Optionally outputs detailed assembly information
// in html (skipped if the html pointer is 0).
void Assembler::assembleAssemblyGraphVertex(
    AssemblyGraph::VertexId vertexId,
    vector<Base>&,
    vector<uint32_t>& repeatCounts,
    ostream* htmlPointer)
{
    const auto k = assemblerInfo->k;

    // The edges of this chain in the marker graph.
    const MemoryAsContainer<GlobalMarkerGraphEdgeId> edgeIds = assemblyGraph.vertices[vertexId];
    const size_t edgeCount = edgeIds.size();
    const size_t vertexCount = edgeCount + 1;

    // The vertices of this chain in the marker graph.
    vector<GlobalMarkerGraphVertexId> vertexIds;
    vertexIds.reserve(vertexCount);
    for(const GlobalMarkerGraphEdgeId edgeId: edgeIds) {
        const MarkerGraphConnectivity::Edge& edge =
            markerGraphConnectivity.edges[edgeId];
        vertexIds.push_back(edge.source);
    }
    const MarkerGraphConnectivity::Edge& lastEdge =
        markerGraphConnectivity.edges[edgeIds[edgeIds.size()-1]];
    vertexIds.push_back(lastEdge.target);

    // Vertex coverage.
    vector<uint32_t> vertexCoverage(vertexCount);
    for(size_t i=0; i<vertexCount; i++) {
        vertexCoverage[i] = uint32_t(globalMarkerGraphVertices.size(vertexIds[i]));
    }

    // Edge coverage.
    vector<uint32_t> edgeCoverage(edgeCount);
    for(size_t i=0; i<edgeCount; i++) {
        edgeCoverage[i] = uint32_t(markerGraphConnectivity.edgeMarkerIntervals.size(edgeIds[i]));
    }




    // Compute consensus sequence for the vertices of the chain.
    vector< vector<Base> > vertexSequences(vertexCount);
    vector< vector<uint32_t> > vertexRepeatCounts(vertexCount);
    for(size_t i=0; i<vertexCount; i++) {
        computeMarkerGraphVertexConsensusSequence(
            vertexIds[i], vertexSequences[i], vertexRepeatCounts[i]);
    }

    // Compute consensus sequence for the edges of the chain.
    vector< vector<Base> > edgeSequences(edgeCount);
    vector< vector<uint32_t> > edgeRepeatCounts(edgeCount);
    for(size_t i=0; i<edgeCount; i++) {
        computeMarkerGraphEdgeConsensusSequence(
            edgeIds[i], edgeSequences[i], edgeRepeatCounts[i]);
    }


    // Compute offsets for edges and vertices.
    // A vertex/edge offset is the position of the first base
    // of the vertex/edge consensus sequence (run-length)
    // relative to the first base of assembled sequence (run-length).
    vector<uint32_t> vertexOffsets(vertexCount);
    vector<uint32_t> edgeOffsets(edgeCount);
    vertexOffsets[0] = 0;
    for(size_t i=0; i<edgeCount; i++) {
        edgeOffsets[i] = vertexOffsets[i];
        vertexOffsets[i+1] = edgeOffsets[i] + uint32_t(edgeSequences[i].size()) - uint32_t(k);
    }



    // Compute, for each vertex, the portion of vertex sequence that contributes
    // to the assembly. This is the portion that does not overlap a vertex with greater coverage.
    // (Break ties using vertex ids).
    // An edge always contributes to the assembly all but its first and last k bases
    // (this can mean no bases, for edges between overlapping vertices).
    vector< pair<uint32_t, uint32_t> > vertexAssembledPortion(vertexCount);
    for(int i=0; i<int(vertexCount); i++) {

        // Check previous vertices.
        vertexAssembledPortion[i].first = 0;
        for(int j=i-1; j>=0; j--) {
            if(vertexOffsets[j]+k < vertexOffsets[i]) {
                break;
            }
            if(vertexCoverage[j]>vertexCoverage[i] ||
                (vertexCoverage[j]==vertexCoverage[i] && vertexIds[j]<vertexIds[i])) {
                vertexAssembledPortion[i].first = vertexOffsets[j] + uint32_t(k) - vertexOffsets[i];
                break;
            }
        }

        // Check following vertices.
        vertexAssembledPortion[i].second = uint32_t(k);
        for(int j=i+1; j<int(vertexCount); j++) {
            if(vertexOffsets[i]+k < vertexOffsets[j]) {
                break;
            }
            if(vertexCoverage[j]>vertexCoverage[i] ||
                (vertexCoverage[j]==vertexCoverage[i] && vertexIds[j]<vertexIds[i])) {
                vertexAssembledPortion[i].second = vertexOffsets[i] + uint32_t(k) - vertexOffsets[j];
                break;
            }
        }
    }



    // If requested, write out details in html format.
    if(htmlPointer) {
        ostream& html = *htmlPointer;

        // Write a title.
        html <<
            "<h1>Assembly graph vertex <a href="
            "'exploreAssemblyGraph?vertexId=" << vertexId <<
            "&maxDistance=6&detailed=on&sizePixels=1600&timeout=30'>" <<
            vertexId << "</a></h1>";

        // Write parent and child vertices in the assembly graph.
        html << "<p>Parent vertices in the assembly graph:";
        for(const auto parentEdge: assemblyGraph.edgesByTarget[vertexId]) {
            const AssemblyGraph::VertexId parent = assemblyGraph.edges[parentEdge].source;
            html <<
                " <a href='exploreAssemblyGraphVertex?vertexId=" << parent << "'>"
                << parent << "</a>";
        }

        html << "<p>Child vertices in the assembly graph:";
        for(const auto childEdge: assemblyGraph.edgesBySource[vertexId]) {
            const AssemblyGraph::VertexId child = assemblyGraph.edges[childEdge].target;
            html <<
                " <a href='exploreAssemblyGraphVertex?vertexId=" << child << "'>"
                << child << "</a>";
        }



        // Write a table with a row for each marker graph vertex or edge
        // in the marker graph chain.
        html <<
            "<p>This vertex of the assembly graph corresponds to a chain of " <<
            edgeIds.size() << " edges in the marker graph. "
            "The table below shows consensus sequences "
            "for the vertices and edges of this chain of the marker graph. "
            "All vertex and edge ids in the table refer to the marker graph."
            "<p><table><tr><th>Vertex<br>or<br>edge<th>Id<th>Coverage<th>Offset<th>Run-length<br>sequence<th>Sequence";
        const string urlPrefix = "exploreMarkerGraph?vertexId=";
        const string urlSuffix =
            "&maxDistance=5"
            "&detailed=on"
            "&minCoverage=3"
            "&minConsensus=3"
            "&sizePixels=600&timeout=10"
            "&useStoredConnectivity=on"
            "&onlyUseSpanningSubgraphEdges=on"
            "&dontUsePrunedEdges=on"
            "&showVertexId=on";
        for(size_t i=0; ; i++) {

            // Vertex.
            const GlobalMarkerGraphVertexId vertexId = vertexIds[i];
            const string url = urlPrefix + to_string(vertexId) + urlSuffix;
            const vector<Base>& vertexSequence = vertexSequences[i];
            const vector<uint32_t>& vertexRepeatCount = vertexRepeatCounts[i];
            const uint32_t maxVertexRepeatCount =
                *std::max_element(vertexRepeatCount.begin(), vertexRepeatCount.end());
            CZI_ASSERT(maxVertexRepeatCount < 10);  // For now. Add additional code when this fails.
            html <<
                 "<tr><td>Vertex" <<
                "<td class=centered><a href='" << url << "'>" << vertexId << "</a>"
                "<td class=centered>" << vertexCoverage[i] <<
                "<td class=centered>" << vertexOffsets[i] <<
                "<td style='font-family:courier'>";
            for(size_t j=0; j<vertexSequence.size(); j++) {
                const uint32_t repeatCount = vertexRepeatCount[j];
                if(j==vertexAssembledPortion[i].first) {
                    html << "<span style='background-color:LightGreen'>";
                }
                html << vertexSequence[j];
                if(j==vertexAssembledPortion[i].second-1) {
                    html << "</span>";
                }
            }
            html << "<br>";
            for(size_t j=0; j<vertexSequence.size(); j++) {
                const uint32_t repeatCount = vertexRepeatCount[j];
                if(j==vertexAssembledPortion[i].first) {
                    html << "<span style='background-color:LightGreen'>";
                }
                html << repeatCount % 10;
                if(j==vertexAssembledPortion[i].second-1) {
                    html << "</span>";
                }
            }
            html << "<td style='font-family:courier'>";
            for(size_t j=0; j<vertexSequence.size(); j++) {
                if(j==vertexAssembledPortion[i].first) {
                    html << "<span style='background-color:LightGreen'>";
                }
                const Base b = vertexSequence[j];
                const uint32_t repeatCount = vertexRepeatCount[j];
                for(uint32_t k=0; k<repeatCount; k++) {
                    html << b;
                }
                if(j==vertexAssembledPortion[i].second-1) {
                    html << "</span>";
                }
            }

            // This was the last vertex.
            if(i == edgeCount) {
                break;
            }

            // Edge.
            const GlobalMarkerGraphEdgeId edgeId = edgeIds[i];
            const MarkerGraphConnectivity::Edge& edge =
                markerGraphConnectivity.edges[edgeId];
            const string sourceUrl = urlPrefix + to_string(edge.source) + urlSuffix;
            const string targetUrl = urlPrefix + to_string(edge.target) + urlSuffix;
            const vector<Base>& edgeSequence = edgeSequences[i];
            const vector<uint32_t>& edgeRepeatCount = edgeRepeatCounts[i];
            const size_t edgeSequenceLength = edgeSequence.size();
            CZI_ASSERT(edgeRepeatCount.size() == edgeSequenceLength);
            const uint32_t maxEdgeRepeatCount =
                *std::max_element(edgeRepeatCount.begin(), edgeRepeatCount.end());
            CZI_ASSERT(maxEdgeRepeatCount < 10);  // For now. Add additional code when this fails.
            html <<
                "<tr><td>Edge<td class=centered>" << edgeId <<
                "<td class=centered>" << edgeCoverage[i] <<
                "<td class=centered>" << edgeOffsets[i] <<
                "<td style='font-family:courier'>";
            for(size_t j=0; j<edgeSequenceLength; j++) {
                if(edgeSequenceLength>2*k && j==k) {
                    html << "<span style='background-color:LightGreen'>";
                }
                html << edgeSequence[j];
                if(edgeSequenceLength>2*k && j==edgeSequenceLength-k-1) {
                    html << "</span>";
                }
            }
            html << "<br>";
            for(size_t j=0; j<edgeSequenceLength; j++) {
                const uint32_t repeatCount = edgeRepeatCount[j];
                if(edgeSequenceLength>2*k && j==k) {
                    html << "<span style='background-color:LightGreen'>";
                }
                html << repeatCount % 10;
                if(edgeSequenceLength>2*k && j==edgeSequenceLength-k-1) {
                    html << "</span>";
                }
            }
            html << "<td style='font-family:courier'>";
            for(size_t j=0; j<edgeSequenceLength; j++) {
                const Base b = edgeSequence[j];
                const uint32_t repeatCount = edgeRepeatCount[j];
                if(edgeSequenceLength>2*k && j==k) {
                    html << "<span style='background-color:LightGreen'>";
                }
                for(uint32_t k=0; k<repeatCount; k++) {
                    html << b;
                }
                if(edgeSequenceLength>2*k && j==edgeSequenceLength-k-1) {
                    html << "</span>";
                }
            }
         }
    }

}
