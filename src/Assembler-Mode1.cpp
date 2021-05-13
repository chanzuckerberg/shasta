#include "Assembler.hpp"
#include "Mode1-AssemblyGraph.hpp"
using namespace shasta;
using namespace Mode1;

void Assembler::createMode1AssemblyGraph(
    uint64_t minEdgeCoverage,
    uint64_t minEdgeCoveragePerStrand)
{
    mode1Data.assemblyGraphPointer = make_unique<Mode1::AssemblyGraph>(
        minEdgeCoverage,
        minEdgeCoveragePerStrand,
        markers,
        markerGraph);
}
