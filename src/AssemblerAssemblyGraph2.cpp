#include "Assembler.hpp"
#include "AssemblyGraph2.hpp"
#include "Reads.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;



void Assembler::createAssemblyGraph2(
    double bubbleRemovalDiscordantRatioThreshold,
    double bubbleRemovalAmbiguityThreshold,
    uint64_t bubbleRemovalMaxPeriod,
    uint64_t superbubbleRemovalEdgeLengthThreshold,
    uint64_t pruneLength,
    uint64_t phasingMinReadCount,
    bool suppressGfaOutput,
    bool suppressFastaOutput,
    bool suppressDetailedOutput,
    bool suppressPhasedOutput,
    bool suppressHaploidOutput,
    size_t threadCount
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

    cout << timestamp << "Assembler::createAssemblyGraph2 begins." << endl;

    assemblyGraph2Pointer = make_shared<AssemblyGraph2>(
        assemblerInfo->readRepresentation,
        assemblerInfo->k,
        getReads().getFlags(),
        markers,
        markerGraph,
        bubbleRemovalDiscordantRatioThreshold,
        bubbleRemovalAmbiguityThreshold,
        bubbleRemovalMaxPeriod,
        superbubbleRemovalEdgeLengthThreshold,
        pruneLength,
        phasingMinReadCount,
        suppressGfaOutput,
        suppressFastaOutput,
        suppressDetailedOutput,
        suppressPhasedOutput,
        suppressHaploidOutput,
        threadCount
        );

    cout << timestamp << "Assembler::createAssemblyGraph2 ends." << endl;
}
