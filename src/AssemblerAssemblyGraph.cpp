// Shasta.
#include "Assembler.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Standard library.
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

        // Looking up the map gives the childre of this assembly graph vertex.
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
    assemblyGraph.edgesBySource.beginPass2();
    assemblyGraph.edgesByTarget.beginPass2();
    assemblyGraph.edgesBySource.endPass2();
    assemblyGraph.edgesByTarget.endPass2();

    cout << timestamp << "Done creating assembly graph edges." << endl;
}



void Assembler::accessAssemblyGraphEdges()
{
    assemblyGraph.edges.accessExistingReadOnly(
        largeDataName("AssemblyGraphEdges"));
}
