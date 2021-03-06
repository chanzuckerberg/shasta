#include "Assembler.hpp"
#include "AssemblyGraph2.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;



void Assembler::createAssemblyGraph2()
{
    // Check that we have what we need.
    checkMarkerGraphVerticesAreAvailable();
    checkMarkerGraphEdgesIsOpen();
    SHASTA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    SHASTA_ASSERT(markerGraph.reverseComplementEdge.isOpen);
    SHASTA_ASSERT(markerGraph.edgesBySource.isOpen());
    SHASTA_ASSERT(markerGraph.edgesByTarget.isOpen());

    assemblyGraph2Pointer = make_shared<AssemblyGraph2>(markerGraph);
    assemblyGraph2Pointer->writeCsv("AssemblyGraph2");
}
