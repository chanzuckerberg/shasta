#include "Assembler.hpp"
#include "mode3.hpp"
#include "mode3-LocalAssemblyGraph.hpp"
using namespace shasta;
using namespace mode3;


void Assembler::exploreMode3AssemblyGraph(
    const vector<string>& request,
    ostream& html)
{
    html << "Not implemented.";
    SHASTA_ASSERT(assemblyGraph3Pointer);

    const mode3::LocalAssemblyGraph localAssemblyGraph(
        *assemblyGraph3Pointer, 200, 10);
    localAssemblyGraph.writeGraphviz("LocalAssemblyGraph.dot");

}
