#include "Assembler.hpp"
#include "CompressedAssemblyGraph.hpp"
using namespace shasta;


void Assembler::createCompressedAssemblyGraph()
{
    compressedAssemblyGraph = make_shared<CompressedAssemblyGraph>(*this);
    CompressedAssemblyGraph& graph = *compressedAssemblyGraph;

    // GFA output (without sequence).
    const double basesPerMarker =
        double(assemblerInfo->baseCount) /
        double(markers.totalSize()/2);
    graph.writeGfa("CompressedAssemblyGraph.gfa", basesPerMarker);

    // Write everything in csv format.
    graph.writeCsv();
}


void Assembler::colorCompressedAssemblyGraph(const string& s)
{
    compressedAssemblyGraph->color(s,
        "CompressedAssemblyGraphColor.csv");
}
