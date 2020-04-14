#include "Assembler.hpp"
#include "CompressedAssemblyGraph.hpp"
using namespace shasta;


void Assembler::createCompressedAssemblyGraph()
{
    CompressedAssemblyGraph graph(*assemblyGraphPointer);

    // GFA output (without sequence).
    const double basesPerMarker =
        double(assemblerInfo->baseCount) /
        double(markers.totalSize()/2);
    graph.writeGfa("CompressedAssemblyGraph.gfa", basesPerMarker);
    graph.writeHtml("CompressedAssemblyGraph.html");
}
