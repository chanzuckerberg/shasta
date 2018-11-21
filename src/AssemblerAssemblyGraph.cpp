// Shasta.
#include "Assembler.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Standard library.
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
    vector<EdgeId> nextVertices;
    vector<EdgeId> previousVertices;
    vector<EdgeId> chain;



    // Main loop over all edges of the pruned spanning subgraph
    // of the marker graph.
    // At each iteration we find a new linear chain of edges.
    for(EdgeId startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {
        const auto& startEdge = markerGraphConnectivity.edges[startEdgeId];

        // If this edge is not part of the pruned spanning subgraph, skip it.
        if(!startEdge.isInSpanningSubgraph) {
            continue;
        }
        if(startEdge.wasPruned) {
            continue;
        }

        // If we already found this edge, skip it.
        // It is part of a chain we aready found.
        if(wasFound[startEdgeId]) {
            continue;
        }

        // Follow the chain forward.
        EdgeId edgeId = startEdgeId;
        while(true) {
            edgeId = nextEdgeInMarkerGraphPrunedSpanningSubgraph(edgeId);
            if(edgeId == invalidGlobalMarkerGraphEdgeId) {
                break;
            }
            nextVertices.push_back(edgeId);
        }

        // Follow the chain backward.
        edgeId = startEdgeId;
        while(true) {
            edgeId = previousEdgeInMarkerGraphPrunedSpanningSubgraph(edgeId);
            if(edgeId == invalidGlobalMarkerGraphEdgeId) {
                break;
            }
            previousVertices.push_back(edgeId);
        }

        // Gather the chain.
        copy(previousVertices.rbegin(), previousVertices.rend(), back_inserter< vector<EdgeId> >(chain));
        chain.push_back(startEdgeId);
        copy(nextVertices.begin(), nextVertices.end(), back_inserter< vector<EdgeId> >(chain));

        // Mark all the edges in the chain as found.
        for(const EdgeId edgeId: chain) {
            wasFound[edgeId] = true;
        }

        // Store this chain as a new vertex of the assembly graph.
        assemblyGraph.vertices.appendVector(chain);

        // Cleanup.
        nextVertices.clear();
        previousVertices.clear();
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
