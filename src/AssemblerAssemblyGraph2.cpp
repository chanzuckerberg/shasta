#include "Assembler.hpp"
#include "AssemblyGraph2.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"
using namespace shasta;



void Assembler::createAssemblyGraph2(
    uint64_t pruneLength,
    const Mode2AssemblyOptions& mode2Options,
    size_t threadCount,
    bool debug
    )
{
    // Check that we have what we need.
    checkMarkerGraphVerticesAreAvailable();
    checkMarkerGraphEdgesIsOpen();
    SHASTA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    SHASTA_ASSERT(markerGraph.reverseComplementEdge.isOpen);
    SHASTA_ASSERT(markerGraph.edgesBySource.isOpen());
    SHASTA_ASSERT(markerGraph.edgesByTarget.isOpen());

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    performanceLog << timestamp << "Assembler::createAssemblyGraph2 begins." << endl;

    assemblyGraph2Pointer = make_shared<AssemblyGraph2>(
        assemblerInfo->readRepresentation,
        assemblerInfo->k,
        getReads().getFlags(),
        markers,
        markerGraph,
        pruneLength,
        mode2Options,
        assemblerInfo->assemblyGraph2Statistics,
        threadCount,
        debug
        );

    performanceLog << timestamp << "Assembler::createAssemblyGraph2 ends." << endl;
}
