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
        directedReadGraph.getVertex(OrientedReadId(readId, 0).getValue()).baseCount = baseCount;
        directedReadGraph.getVertex(OrientedReadId(readId, 1).getValue()).baseCount = baseCount;
    }

    // Add a pair of edges for each stored alignment.
    directedReadGraph.edges.reserve(2 * alignmentData.size());
    for(const AlignmentData& alignment: alignmentData) {
        directedReadGraph.addEdgePair(alignment);
    }

    // Make sure the read graph is invariant under reverse complementing.
    directedReadGraph.check();


    SHASTA_ASSERT(0);
}



void Assembler::accessDirectedReadGraphReadOnly()
{
    directedReadGraph.accessExistingReadOnly(largeDataName("DirectedReadGraph"));
}
void Assembler::accessDirectedReadGraphReadWrite()
{
    directedReadGraph.accessExistingReadWrite(largeDataName("DirectedReadGraph"));
}
