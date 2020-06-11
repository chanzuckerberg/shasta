#include "MarkerGraph.hpp"
using namespace shasta;

const MarkerGraph::VertexId MarkerGraph::invalidVertexId = std::numeric_limits<VertexId>::max();
const MarkerGraph::EdgeId MarkerGraph::invalidEdgeId = std::numeric_limits<EdgeId>::max();
const MarkerGraph::CompressedVertexId
    MarkerGraph::invalidCompressedVertexId = std::numeric_limits<uint64_t>::max();



MarkerGraph::MarkerGraph() :
    MultithreadedObject<MarkerGraph>(*this)
    {}



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
