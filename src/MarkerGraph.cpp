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
MarkerGraph::VertexId MarkerGraph::renumberVertexTable(size_t threadCount)
{
    // Sanity check.
    SHASTA_ASSERT(threadCount > 0);
    SHASTA_ASSERT(vertexTable.isOpen);
    SHASTA_ASSERT(vertexTable.size() > 0);

    // Find the maximum vertex id.
    const VertexId maxVertexId = findMaxVertexTableEntry(threadCount);

    // Call the lower level version.
    return renumberVertexTable(threadCount, maxVertexId);

}



// This second version can be called if the maximum vertex id
// present in the vertex table is already known, and is faster.
MarkerGraph::VertexId MarkerGraph::renumberVertexTable(size_t threadCount, VertexId maxVertexId)
{
    const bool debug = false;

    // Sanity check.
    SHASTA_ASSERT(threadCount > 0);
    SHASTA_ASSERT(vertexTable.isOpen);
    SHASTA_ASSERT(vertexTable.size() > 0);

    if(debug) {
        ofstream csv("VertexTable-Before.csv");
        csv << "MarkerId,VertexId\n";
        for(MarkerId markerId=0; markerId<vertexTable.size(); markerId++) {
            csv << markerId << "," << vertexTable[markerId] << "\n";
        }
    }

    cout << timestamp << "Renumbering the marker graph vertex table." << endl;

    // Create a vector of bools that tells us which VertexId's are present.
    const string vertexTableName = vertexTable.fileName;
    renumberVertexTableData.isPresent.createNew(
        vertexTableName.empty() ? "" : (vertexTableName + "-tmp-isPresent"),
        vertexTable.getPageSize());
    renumberVertexTableData.isPresent.resize(maxVertexId + 1);
    fill(
        renumberVertexTableData.isPresent.begin(),
        renumberVertexTableData.isPresent.end(),
        false);
    const uint64_t batchSize = 100000;
    setupLoadBalancing(vertexTable.size(), batchSize);
    runThreads(&MarkerGraph::renumberVertexTableThreadFunction1, threadCount);

    // Now we know what VertexId's are present, so we can compute the new VertexId
    // corresponding to each old VertexId.
    renumberVertexTableData.newVertexId.createNew(
        vertexTableName.empty() ? "" : (vertexTableName + "-tmp-newVertexId"),
        vertexTable.getPageSize());
    renumberVertexTableData.newVertexId.resize(maxVertexId + 1);
    VertexId newVertexId = 0;
    for(VertexId oldVertexId=0; oldVertexId<=maxVertexId; oldVertexId++) {
        if(renumberVertexTableData.isPresent[oldVertexId]) {
            renumberVertexTableData.newVertexId[oldVertexId] = newVertexId;
            ++newVertexId;
        } else {
            renumberVertexTableData.newVertexId[oldVertexId] = invalidVertexId;
        }
    }

    // Now we can renumber the vertex table.
    setupLoadBalancing(vertexTable.size(), batchSize);
    runThreads(&MarkerGraph::renumberVertexTableThreadFunction2, threadCount);

    // Clean up.
    renumberVertexTableData.newVertexId.remove();
    renumberVertexTableData.isPresent.remove();

    if(debug) {
        ofstream csv("VertexTable-After.csv");
        csv << "MarkerId,VertexId\n";
        for(MarkerId markerId=0; markerId<vertexTable.size(); markerId++) {
            csv << markerId << "," << vertexTable[markerId] << "\n";
        }
    }

    cout << timestamp << "Done renumbering the marker graph vertex table." << endl;
    return newVertexId - 1;
}



void MarkerGraph::renumberVertexTableThreadFunction1(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertex table entries in this batch.
        for(uint64_t markerId=begin; markerId!=end; markerId++) {
            const CompressedVertexId compressedVertexId = vertexTable[markerId];
            if(compressedVertexId != invalidCompressedVertexId) {
                renumberVertexTableData.isPresent[VertexId(compressedVertexId)] = true;
            }
        }
    }
}



void MarkerGraph::renumberVertexTableThreadFunction2(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertex table entries in this batch.
        for(uint64_t markerId=begin; markerId!=end; markerId++) {
            const CompressedVertexId compressedVertexId = vertexTable[markerId];
            if(compressedVertexId != invalidCompressedVertexId) {
                const VertexId oldVertexId = VertexId(compressedVertexId);
                const VertexId newVertexId = renumberVertexTableData.newVertexId[oldVertexId];
                vertexTable[markerId] = CompressedVertexId(newVertexId);
            }
        }
    }
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




// Recreate the vertices from the vertexTable.
// This assumes that valid VertexId's in the vertex table
// are numbered contiguously starting at 0 (call renumberVertexTable to ensure that).
void MarkerGraph::createVerticesFromVertexTable(size_t threadCount, VertexId maxVertexId)
{
    SHASTA_ASSERT(vertexTable.isOpen);
    const string vertexTableName = vertexTable.fileName;

    // Create our copy of the vertices.
    createVerticesFromVertexTableData.vertices.createNew(
        vertexTableName.empty() ? "" : (vertexTableName + "-tmp-isPresent"),
        vertexTable.getPageSize());

    // Pass 1: count the number of markers in each vertex.
    createVerticesFromVertexTableData.vertices.beginPass1(maxVertexId+1);
    uint64_t batchSize = 100000;
    setupLoadBalancing(vertexTable.size(), batchSize);
    runThreads(&MarkerGraph::createVerticesFromVertexTableThreadFunction1, threadCount);

    // Pass 2: fill in the markers of each vertex.
    createVerticesFromVertexTableData.vertices.beginPass2();
    setupLoadBalancing(vertexTable.size(), batchSize);
    runThreads(&MarkerGraph::createVerticesFromVertexTableThreadFunction2, threadCount);
    createVerticesFromVertexTableData.vertices.endPass2(false);

    // Pass 3: sort the markers of each vertex.
    batchSize = 1000;
    setupLoadBalancing(maxVertexId+1, batchSize);
    runThreads(&MarkerGraph::createVerticesFromVertexTableThreadFunction3, threadCount);

    // Finally, we copy to the main copy of the vertices.
    vertices().clear();
    const uint64_t vertexCount = createVerticesFromVertexTableData.vertices.size();
    for(VertexId vertexId=0; vertexId<vertexCount; vertexId++) {
        vertices().appendVector(createVerticesFromVertexTableData.vertices.size(vertexId));
    }
    setupLoadBalancing(vertexCount, batchSize);
    runThreads(&MarkerGraph::createVerticesFromVertexTableThreadFunction4, threadCount);

    // Cleanup.
    createVerticesFromVertexTableData.vertices.remove();
}



void MarkerGraph::createVerticesFromVertexTableThreadFunction1(size_t threadId)
{
    auto& vertices = createVerticesFromVertexTableData.vertices;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertex table entries in this batch.
        for(uint64_t markerId=begin; markerId!=end; markerId++) {
            const CompressedVertexId compressedVertexId = vertexTable[markerId];

            if(compressedVertexId != invalidCompressedVertexId) {
                vertices.incrementCountMultithreaded(VertexId(compressedVertexId));
            }
        }
    }

}



void MarkerGraph::createVerticesFromVertexTableThreadFunction2(size_t threadId)
{
    auto& vertices = createVerticesFromVertexTableData.vertices;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertex table entries in this batch.
        for(uint64_t markerId=begin; markerId!=end; markerId++) {
            const CompressedVertexId compressedVertexId = vertexTable[markerId];

            if(compressedVertexId != invalidCompressedVertexId) {
                vertices.storeMultithreaded(VertexId(compressedVertexId), markerId);
            }
        }
    }

}



void MarkerGraph::createVerticesFromVertexTableThreadFunction3(size_t threadId)
{
    auto& vertices = createVerticesFromVertexTableData.vertices;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over verteices in this batch.
        for(VertexId vertexId=begin; vertexId!=end; vertexId++) {
            auto vertexMarkers = vertices[vertexId];
            sort(vertexMarkers.begin(), vertexMarkers.end());
        }
    }

}



void MarkerGraph::createVerticesFromVertexTableThreadFunction4(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices in this batch.
        for(VertexId vertexId=begin; vertexId!=end; vertexId++) {

            // Make the copy.
            auto vertexMarkers = createVerticesFromVertexTableData.vertices[vertexId];
            copy(vertexMarkers.begin(), vertexMarkers.end(), verticesPointer->begin(vertexId));
        }
    }

}


