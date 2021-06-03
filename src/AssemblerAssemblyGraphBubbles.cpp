#include "Assembler.hpp"
#include "Bubbles.hpp"
using namespace shasta;


void Assembler::analyzeAssemblyGraphBubbles()
{
    // Check that we have what we need.
    SHASTA_ASSERT(assemblyGraphPointer);
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    SHASTA_ASSERT(assemblyGraph.vertices.isOpen);
    SHASTA_ASSERT(assemblyGraph.edges.isOpen);
    SHASTA_ASSERT(assemblyGraph.edgesBySource.isOpen());

    // Create the Bubbles.
    Bubbles bubbles(assemblyGraph);
}
