#include "Assembler.hpp"
#include "mode3.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3;



void Assembler::mode3Assembly(
    size_t threadCount)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const DynamicAssemblyGraph dynamicAssemblyGraph(
        reads->getFlags(),
        markers,
        markerGraph,
        largeDataFileNamePrefix,
        largeDataPageSize,
        threadCount);

    mode3AssemblyGraph = std::make_shared<mode3::AssemblyGraph>(
        dynamicAssemblyGraph,
        largeDataFileNamePrefix,
        largeDataPageSize);

    const mode3::LocalAssemblyGraph localAssemblyGraph(
        *mode3AssemblyGraph, 200, 10);
    localAssemblyGraph.writeGraphviz("LocalAssemblyGraph.dot");

}
