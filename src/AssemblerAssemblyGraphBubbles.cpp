#include "Assembler.hpp"
#include "Bubbles.hpp"
using namespace shasta;


void Assembler::analyzeAssemblyGraphBubbles(bool debug)
{
    // Check that we have what we need.
    SHASTA_ASSERT(assemblyGraphPointer);
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    SHASTA_ASSERT(assemblyGraph.vertices.isOpen);
    SHASTA_ASSERT(assemblyGraph.edges.isOpen);
    SHASTA_ASSERT(assemblyGraph.edgesBySource.isOpen());
    SHASTA_ASSERT(assemblyGraph.edgeLists.isOpen());

    // Create the Bubbles.
    bubbles = make_shared<Bubbles>(*this, debug);
}
