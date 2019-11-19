#include "Assembler.hpp"
using namespace shasta;


void Assembler::createDirectedReadGraph()
{
    // Initialize the directed read graph.
    directedReadGraph.createNew(largeDataName("DirectedReadGraph"), largeDataPageSize);
    directedReadGraph.createVertices(readCount());

    // Store the number of bases in each oriented read.
    for(ReadId readId=0; readId<readCount(); readId++) {
        const uint64_t baseCount = getReadRawSequenceLength(readId);
        const uint64_t markerCount = markers.size(OrientedReadId(readId, 0).getValue());
        DirectedReadGraphVertex& vertex0 = directedReadGraph.getVertex(OrientedReadId(readId, 0).getValue());
        DirectedReadGraphVertex& vertex1 = directedReadGraph.getVertex(OrientedReadId(readId, 1).getValue());
        vertex0.baseCount = baseCount;
        vertex1.baseCount = baseCount;
        vertex0.markerCount = markerCount;
        vertex1.markerCount = markerCount;
    }

    // Add a pair of edges for each stored alignment.
    directedReadGraph.edges.reserve(2 * alignmentData.size());
    for(const AlignmentData& alignment: alignmentData) {
        directedReadGraph.addEdgePair(alignment);
    }

    // Compute graph connectivity.
    directedReadGraph.computeConnectivity();

    // Make sure the read graph is invariant under reverse complementing.
    directedReadGraph.check();


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
}



void Assembler::accessDirectedReadGraphReadOnly()
{
    directedReadGraph.accessExistingReadOnly(largeDataName("DirectedReadGraph"));
}
void Assembler::accessDirectedReadGraphReadWrite()
{
    directedReadGraph.accessExistingReadWrite(largeDataName("DirectedReadGraph"));
}



void Assembler::directedReadGraphTransitiveReduction(
    double offsetTolerance0,
    double offsetTolerance1)
{
    SHASTA_ASSERT(directedReadGraph.isOpenWithWriteAccess());
    directedReadGraph.transitiveReduction(offsetTolerance0, offsetTolerance1);
}
