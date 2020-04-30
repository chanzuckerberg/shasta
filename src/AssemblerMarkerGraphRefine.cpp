#include "Assembler.hpp"
using namespace shasta;



// Refine the marker graph by removing vertices in tangle regions,
// then recreating edges. This must be called after
// transitive reduction. After this is called, the only
// two MarkerGraph field filled in are vertices and vertexTable.
// Everything else has to be recreated.
void Assembler::refineMarkerGraph(
    uint64_t refineThreshold,
    size_t threadCount)
{
    cout << timestamp << "Refine marker graph begins." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Check that we have what we need.
    checkMarkerGraphVerticesAreAvailable();
    SHASTA_ASSERT(markerGraph.vertexTable.isOpenWithWriteAccess);
    checkMarkerGraphEdgesIsOpen();
    SHASTA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    SHASTA_ASSERT(markerGraph.reverseComplementEdge.isOpen);
    cout << "The marker graph has " << markerGraph.vertices.size() <<
        " vertices and " << markerGraph.edges.size() << " edges." << endl;
    SHASTA_ASSERT(markerGraph.reverseComplementVertex.size() == markerGraph.vertices.size());
    SHASTA_ASSERT(markerGraph.reverseComplementEdge.size() == markerGraph.edges.size());

    // Create a temporary assembly graph.
    cout << timestamp << "Creating a temporary assembly graph." << endl;
    createAssemblyGraphEdges();
    createAssemblyGraphVertices();
    AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    cout << "The assembly graph has " << assemblyGraph.vertices.size() <<
        " vertices and " << assemblyGraph.edges.size() << " edges." << endl;

    // Vector to flag the vertices we want to remove.
    MemoryMapped::Vector<bool> isVertexToBeRemoved;
    isVertexToBeRemoved.createNew(
        largeDataName("tmp-VerticesToBeRemoved"), largeDataPageSize);
    isVertexToBeRemoved.resize(markerGraph.vertices.size());
    fill(isVertexToBeRemoved.begin(), isVertexToBeRemoved.end(), false);

    // Flag to be removed all marker graph vertices with out-degree
    // or in-degree greater than 1.
    uint64_t removedDueToDegreeCount = 0;
    for(MarkerGraph::VertexId vertexId=0;
        vertexId<markerGraph.vertices.size(); vertexId++) {
        if(
            (markerGraph.outDegree(vertexId) > 1) or
            (markerGraph.inDegree(vertexId) > 1)) {
            isVertexToBeRemoved[vertexId] = true;
            ++removedDueToDegreeCount;
        }
    }
    cout << removedDueToDegreeCount << " marker graph vertices "
        "were flagged for removal because they have "
        "out-degree or in-degree greater than 1." << endl;



    // Flag to be removed all marker graph vertices internal
    // to short assembly graph edges.
    for(AssemblyGraph::EdgeId aEdgeId=0;
        aEdgeId<assemblyGraph.edgeLists.size(); aEdgeId++) {

        // The marker graph edges corresponding to this assembly graph edge.
        const span<MarkerGraph::EdgeId> mEdgeIds =
            assemblyGraph.edgeLists[aEdgeId];

        // If long enough, skip.
        if(mEdgeIds.size() >= refineThreshold) {
            continue;
        }

        // This assembly graph edge is short.
        // Flag to be removed all marker graph vertices internal to it.
        for(uint64_t i=1; i<mEdgeIds.size(); i++) {
            const MarkerGraph::EdgeId mEdgeId = mEdgeIds[i];
            const MarkerGraph::VertexId mVertexId = markerGraph.edges[mEdgeId].source;
            SHASTA_ASSERT(mVertexId < isVertexToBeRemoved.size());
            isVertexToBeRemoved[mVertexId] = true;
        }
    }



    // Remove the marker graph vertices we flagged.
    MemoryMapped::Vector<MarkerGraph::VertexId> verticesToBeKept;
    verticesToBeKept.createNew(
        largeDataName("tmp-VerticesToBeKept"), largeDataPageSize);
    for(MarkerGraph::VertexId vertexId=0;
        vertexId<markerGraph.vertices.size(); vertexId++) {
        if(not isVertexToBeRemoved[vertexId]) {
            verticesToBeKept.push_back(vertexId);
        }
    }
    const uint64_t removedCount = markerGraph.vertices.size() - verticesToBeKept.size();
    cout << timestamp << "Out of " << markerGraph.vertices.size() <<
        " marker graph vertices, " << removedCount <<
        " will be removed." << endl;
    markerGraph.removeVertices(verticesToBeKept, largeDataPageSize, threadCount);




    // Clean up.
    assemblyGraph.remove();
    isVertexToBeRemoved.remove();
    verticesToBeKept.remove();

    cout << timestamp << "Refine marker graph ends." << endl;
}


