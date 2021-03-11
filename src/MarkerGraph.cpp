#include "MarkerGraph.hpp"
using namespace shasta;

const MarkerGraph::VertexId MarkerGraph::invalidVertexId = std::numeric_limits<VertexId>::max();
const MarkerGraph::EdgeId MarkerGraph::invalidEdgeId = std::numeric_limits<EdgeId>::max();
const MarkerGraph::CompressedVertexId
    MarkerGraph::invalidCompressedVertexId = std::numeric_limits<uint64_t>::max();



MarkerGraph::MarkerGraph() :
    MultithreadedObject<MarkerGraph>(*this)
    {}



void MarkerGraph::remove()
{
    destructVertices();
    if(vertexTable.isOpen) {
        vertexTable.remove();
    }
    if(reverseComplementVertex.isOpen) {
        reverseComplementVertex.remove();
    }
    if(edges.isOpen) {
        edges.remove();
    }
    if(reverseComplementEdge.isOpen) {
        reverseComplementEdge.remove();
    }
    if(edgeMarkerIntervals.isOpen()) {
        edgeMarkerIntervals.remove();
    }
    if(edgesBySource.isOpen()) {
        edgesBySource.remove();
    }
    if(edgesByTarget.isOpen()) {
        edgesByTarget.remove();
    }
    if(vertexRepeatCounts.isOpen) {
        vertexRepeatCounts.remove();
    }
    if(edgeConsensus.isOpen()) {
        edgeConsensus.remove();
    }
    if(edgeConsensusOverlappingBaseCount.isOpen) {
        edgeConsensusOverlappingBaseCount.remove();
    }
    if(vertexCoverageData.isOpen()) {
        vertexCoverageData.remove();
    }
    if(edgeCoverageData.isOpen()) {
        edgeCoverageData.remove();
    }
}



// Locate the edge given the vertices.
const MarkerGraph::Edge*
    MarkerGraph::findEdge(Uint40 source, Uint40 target) const
{
    const auto edgesWithThisSource = edgesBySource[source];
    for(const uint64_t i: edgesWithThisSource) {
        const Edge& edge = edges[i];
        if(edge.target == target) {
            return &edge;
        }
    }
    return 0;
}


MarkerGraph::EdgeId MarkerGraph::findEdgeId(Uint40 source, Uint40 target) const
{
    const Edge *edgePointer = findEdge(source, target);
    SHASTA_ASSERT(edgePointer);
    return edgePointer - edges.begin();
}

// Compute in-degree or out-degree of a vertex,
// counting only edges that were not removed.
uint64_t MarkerGraph::inDegree(VertexId vertexId) const
{
    uint64_t degree = 0;
    for(const EdgeId edgeId: edgesByTarget[vertexId]) {
        if(not edges[edgeId].wasRemoved()) {
            ++degree;
        }
    }
    return degree;
}
uint64_t MarkerGraph::outDegree(VertexId vertexId) const
{
    uint64_t degree = 0;
    for(const EdgeId edgeId: edgesBySource[vertexId]) {
        if(not edges[edgeId].wasRemoved()) {
            ++degree;
        }
    }
    return degree;
}



// Remove marker graph vertices and update vertices and vertexTable.
void MarkerGraph::removeVertices(
    const MemoryMapped::Vector<VertexId>& verticesToBeKept,
    uint64_t pageSize,
    uint64_t threadCount)
{
    // Get the names of the data structures we will work with.
    const string verticesName = vertices().getName();
    const string vertexTableName = vertexTable.fileName;
    const string newVerticesName =
        verticesName.empty() ? "" : (verticesName + "-tmp");
    const string newVertexTableName =
        vertexTableName.empty() ? "" : (vertexTableName + "-tmp");

    // Create the new vertices.
    removeVerticesData.verticesToBeKept = &verticesToBeKept;
    removeVerticesData.newVerticesPointer =
        make_shared<MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId> >();
    MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId>& newVertices
        = *removeVerticesData.newVerticesPointer;
    newVertices.createNew(newVerticesName, pageSize);
    newVertices.beginPass1(verticesToBeKept.size());
    const uint64_t batchCount = 10000;
    setupLoadBalancing(verticesToBeKept.size(), batchCount);
    runThreads(&MarkerGraph::removeVerticesThreadFunction1, threadCount);
    newVertices.beginPass2();
    setupLoadBalancing(verticesToBeKept.size(), batchCount);
    runThreads(&MarkerGraph::removeVerticesThreadFunction2, threadCount);
    newVertices.endPass2(false, true);



    // Replace the old vertices with the new ones.
    vertices().remove();
    destructVertices();
    newVertices.rename(verticesName);
    verticesPointer = removeVerticesData.newVerticesPointer;
    removeVerticesData.newVerticesPointer = 0;


    // Update the vertexTable, in place.
    fill(vertexTable.begin(), vertexTable.end(), invalidCompressedVertexId);
    setupLoadBalancing(vertexCount(), batchCount);
    runThreads(&MarkerGraph::removeVerticesThreadFunction3, threadCount);



    // Remove everything else.
    if(reverseComplementVertex.isOpen) {
        reverseComplementVertex.remove();
    }
    if(edges.isOpen) {
        edges.remove();
    }
    if(edgeMarkerIntervals.isOpen()) {
        edgeMarkerIntervals.remove();
    }
    if(edgesBySource.isOpen()) {
        edgesBySource.remove();
    }
    if(edgesByTarget.isOpen()) {
        edgesByTarget.remove();
    }
    if(reverseComplementEdge.isOpen) {
        reverseComplementEdge.remove();
    }
    if(vertexRepeatCounts.isOpen) {
        vertexRepeatCounts.remove();
    }
    if(edgeConsensus.isOpen()) {
        edgeConsensus.remove();
    }
    if(edgeConsensusOverlappingBaseCount.isOpen) {
        edgeConsensusOverlappingBaseCount.remove();
    }
    if(vertexCoverageData.isOpen()) {
        vertexCoverageData.remove();
    }
    if(edgeCoverageData.isOpen()) {
        edgeCoverageData.remove();
    }
}



void MarkerGraph::removeVerticesThreadFunction1(size_t threadId)
{
    const MemoryMapped::Vector<VertexId>& verticesToBeKept =
        *removeVerticesData.verticesToBeKept;
    MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId>& newVertices
        = *removeVerticesData.newVerticesPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices assigned to this thread.
        for(VertexId newVertexId=begin; newVertexId!=end; newVertexId++) {
            const VertexId oldVertexId = verticesToBeKept[newVertexId];
            newVertices.incrementCount(newVertexId, vertexCoverage(oldVertexId));
        }
    }
}



void MarkerGraph::removeVerticesThreadFunction2(size_t threadId)
{
    const MemoryMapped::Vector<VertexId>& verticesToBeKept =
        *removeVerticesData.verticesToBeKept;
    MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId>& newVertices
        = *removeVerticesData.newVerticesPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices assigned to this thread.
        for(VertexId newVertexId=begin; newVertexId!=end; newVertexId++) {
            const VertexId oldVertexId = verticesToBeKept[newVertexId];
            copy(vertices().begin(oldVertexId), vertices().end(oldVertexId),
                newVertices.begin(newVertexId));
        }
    }
}



void MarkerGraph::removeVerticesThreadFunction3(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices assigned to this thread.
        for(VertexId vertexId=begin; vertexId!=end; vertexId++) {
            const CompressedVertexId compressedVertexId = vertexId;
            const span<MarkerId> vertexMarkerIds = vertices()[vertexId];
            for(const MarkerId markerId: vertexMarkerIds) {
                vertexTable[markerId] = compressedVertexId;
            }
        }
    }
}



// This renumbers the vertex table to make sure that
// vertices are numbered contiguously starting at 0.
// This must be called after the vertexTable is changed,
// as in Assembler::cleanupDuplicateMarkers.
// After this is called, all other data structures
// are inconsistent and need to be recreated.
void MarkerGraph::renumberVertexTable(size_t threadCount)
{
    // Sanity check.
    SHASTA_ASSERT(threadCount > 0);
    SHASTA_ASSERT(vertexTable.isOpen);
    SHASTA_ASSERT(vertexTable.size() > 0);

    // Find the maximum vertex id.
    const VertexId maxVertexId = findMaxVertexTableEntry(threadCount);

    // Create a vector of bools that tells us which VertexId's are present.
    const string vertexTableName = vertexTable.fileName;
    renumberVertexTableData.isPresent.createNew(
        vertexTableName.empty() ? "" : (vertexTableName + "-tmp-isPresent"),
        vertexTable.getPageSize());
    renumberVertexTableData.isPresent.resize(maxVertexId - 1);
    fill(
        renumberVertexTableData.isPresent.begin(),
        renumberVertexTableData.isPresent.end(),
        false);


    // Clean up.
    renumberVertexTableData.isPresent.remove();
}



MarkerGraph::VertexId MarkerGraph::findMaxVertexTableEntry(size_t threadCount)
{
    // Sanity checks.
    SHASTA_ASSERT(threadCount > 0);
    SHASTA_ASSERT(vertexTable.isOpen);

    // Initialize the maximum VertexId found by each thread.
    findMaxVertexTableEntryData.threadMaxVertexId.resize(threadCount);
    fill(
        findMaxVertexTableEntryData.threadMaxVertexId.begin(),
        findMaxVertexTableEntryData.threadMaxVertexId.end(),
        0);

    // Each thread finds the maximum of a subset of the vertex table.
    const uint64_t batchSize = 100000;
    setupLoadBalancing(vertexTable.size(), batchSize);
    runThreads(&MarkerGraph::findMaxVertexTableEntryThreadFunction, threadCount);

    // Return the maximum value found by all threads.
    return *std::max_element(
        findMaxVertexTableEntryData.threadMaxVertexId.begin(),
        findMaxVertexTableEntryData.threadMaxVertexId.end());
}



void MarkerGraph::findMaxVertexTableEntryThreadFunction(size_t threadId)
{
    VertexId maxVertexId = 0;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertex table entries in this batch.
        for(uint64_t markerId=begin; markerId!=end; markerId++) {
            const CompressedVertexId compressedVertexId = vertexTable[markerId];
            if(compressedVertexId != invalidCompressedVertexId) {
                maxVertexId = max(maxVertexId, VertexId(compressedVertexId));
            }
        }
    }

    findMaxVertexTableEntryData.threadMaxVertexId[threadId] = maxVertexId;
}


