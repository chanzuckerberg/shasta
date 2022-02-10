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

    assemblyGraph3Pointer = std::make_shared<mode3::AssemblyGraph>(
        dynamicAssemblyGraph,
        largeDataFileNamePrefix,
        largeDataPageSize);
    assemblyGraph3Pointer->writeGfa("AssemblyGraph.gfa");

}



void Assembler::accessMode3AssemblyGraph()
{
    assemblyGraph3Pointer = std::make_shared<mode3::AssemblyGraph>(largeDataFileNamePrefix);
}

