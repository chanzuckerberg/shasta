#include "Assembler.hpp"
#include "mode3.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3;



void Assembler::mode3Assembly(
    size_t threadCount)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minClusterSize = 3;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    assemblyGraph3Pointer = std::make_shared<mode3::AssemblyGraph>(
        largeDataFileNamePrefix,
        largeDataPageSize,
        threadCount,
        markers,
        markerGraph);
    auto& assemblyGraph3 = *assemblyGraph3Pointer;
    assemblyGraph3.writeGfa("AssemblyGraph.gfa");
    assemblyGraph3.clusterSegments(threadCount, minClusterSize);

}



void Assembler::accessMode3AssemblyGraph()
{
    assemblyGraph3Pointer = std::make_shared<mode3::AssemblyGraph>(largeDataFileNamePrefix, markers, markerGraph);
}



void Assembler::analyzeMode3Subgraph(const vector<uint64_t>& segmentIds)
{
    SHASTA_ASSERT(assemblyGraph3Pointer);
    assemblyGraph3Pointer->analyzeSubgraph(segmentIds);
}



void Assembler::createMode3AssemblyPath(
    uint64_t segmentId,
    uint64_t direction,    // 0 = forward, 1 = backward
    vector<uint64_t>& path // The segmentId's of the path.
    ) const
{
    SHASTA_ASSERT(assemblyGraph3Pointer);
    assemblyGraph3Pointer->createAssemblyPath(segmentId, direction, path);

}

