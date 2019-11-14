#include "Assembler.hpp"
using namespace shasta;


void Assembler::createDirectedReadGraph()
{
    // Initialize the directed read graph.
    directedReadGraph.createNew(largeDataName("DirectedReadGraph"), largeDataPageSize);
    directedReadGraph.vertices.resize(2 * readCount());
    directedReadGraph.edges.reserve(2 * alignmentData.size());

    // Add a pair of edges for each stored alignment.
    for(const AlignmentData& alignment: alignmentData) {
        directedReadGraph.addEdgePair(alignment);
    }


    SHASTA_ASSERT(0);
}



