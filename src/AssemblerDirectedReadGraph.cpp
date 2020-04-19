#include "Assembler.hpp"
using namespace shasta;


void Assembler::createDirectedReadGraph(
    uint64_t maxTrim,
    uint64_t containedNeighborCount,
    uint64_t uncontainedNeighborCountPerDirection)
{
    // Initialize the directed read graph.
    directedReadGraph.createNew(largeDataName("DirectedReadGraph"), largeDataPageSize);
    directedReadGraph.createVertices(readCount());

    // Store the number of bases in each oriented read.
    for(ReadId readId=0; readId<readCount(); readId++) {
        const uint32_t baseCount = uint32_t(getReadRawSequenceLength(readId));
        const uint32_t markerCount = uint32_t(markers.size(OrientedReadId(readId, 0).getValue()));
        DirectedReadGraphVertex& vertex0 = directedReadGraph.getVertex(OrientedReadId(readId, 0).getValue());
        DirectedReadGraphVertex& vertex1 = directedReadGraph.getVertex(OrientedReadId(readId, 1).getValue());
        vertex0.baseCount = baseCount;
        vertex1.baseCount = baseCount;
        vertex0.markerCount = markerCount;
        vertex1.markerCount = markerCount;
    }

    // Add a pair of edges for each stored alignment.
    directedReadGraph.edges.reserve(2 * alignmentData.size());
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
        const AlignmentData& alignment = alignmentData[alignmentId];
        directedReadGraph.addEdgePair(alignment, alignmentId);
    }

    // Compute graph connectivity.
    directedReadGraph.computeConnectivity();

    // Flag contained vertices and set edge flags accordingly.
    directedReadGraph.flagContainedVertices(uint32_t(maxTrim));


    // Count the number of isolated reads and their bases.
    uint64_t isolatedReadCount = 0;
    uint64_t isolatedReadBaseCount = 0;
    for(ReadId readId=0; readId<readCount(); readId++) {
        const OrientedReadId orientedReadId(readId, 0);
        const DirectedReadGraph::VertexId vertexId = orientedReadId.getValue();
        if(directedReadGraph.totalDegree(vertexId) > 0) {
            continue;
        }
        ++isolatedReadCount;
        isolatedReadBaseCount += getReadRawSequenceLength(readId);
    }
    assemblerInfo->isolatedReadCount = isolatedReadCount;
    assemblerInfo->isolatedReadBaseCount = isolatedReadBaseCount;

    // Flag edges to be kept.
    // These are the edges that will be used to create the marker graph.
    directedReadGraph. flagEdgesToBeKept(
        containedNeighborCount,
        uncontainedNeighborCountPerDirection);

    // Write a csv file with information on the edges.
    directedReadGraph.writeEdges();

    // Make sure the read graph is invariant under reverse complementing.
    directedReadGraph.check();
}



void Assembler::accessDirectedReadGraphReadOnly()
{
    directedReadGraph.accessExistingReadOnly(largeDataName("DirectedReadGraph"));
}
void Assembler::accessDirectedReadGraphReadWrite()
{
    directedReadGraph.accessExistingReadWrite(largeDataName("DirectedReadGraph"));
}

