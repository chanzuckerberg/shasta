#include "Assembler.hpp"
using namespace shasta;


void Assembler::createDirectedReadGraph()
{
    // Initialize the directed read graph.
    directedReadGraph.createNew(largeDataName("DirectedReadGraph"), largeDataPageSize);
    directedReadGraph.createVertices(readCount());
    directedReadGraph.edges.reserve(2 * alignmentData.size());

    // Add a pair of edges for each stored alignment.
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
