#include "Assembler.hpp"
#include "subgraph.hpp"
using namespace shasta;



void Assembler::exploreCompressedAssemblyGraph(
    const vector<string>& request,
    ostream& html)
{
    html << "<h1>Compressed assembly graph</h1>"
        "<p>Each edge of the compressed assembly graph corresponds to "
        " a linear sequence of bubbles in the assembly graph, "
        "which can degenerate to a single edge.";

    // Create the CompressedAssemblyGraph, if necessary.
    if(not compressedAssemblyGraph) {
        html << "<pre>" << timestamp << "Creating the compressed assembly graph.\n";
        createCompressedAssemblyGraph();
        html << timestamp << "Done creating the compressed assembly graph.";
    }
}
