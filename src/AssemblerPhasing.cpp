#include "Assembler.hpp"
using namespace shasta;


void Assembler::createPhasingGraph(
    size_t threadCount,
    double phasingSimilarityThreshold,
    int maxNeighborCount)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store information used by the phasing graph to create
    // binary data.
    // phasingGraph.dataFileNamePrefix = largeDataFileNamePrefix;
    // phasingGraph.dataPageSize = largeDataPageSize;

    // Find the oriented reads internal to each assembly graph edge.
    phasingGatherOrientedReads(threadCount);

#if 0
    // Find the assembly graph edges that each oriented read is internal to.
    phasingGatherAssemblyGraphEdges(threadCount);
    phasingSortAssemblyGraphEdges(threadCount);

    // Find oriented read pairs with phasing similarity greater than the threshold.
    phasingGraph.findSimilarPairs(threadCount, phasingSimilarityThreshold);

    // Only keep up to maxNeighborCount neighbors.
    phasingGraph.keepBestSimilarPairs(maxNeighborCount);

    // Write out the global phasing graph in graphviz format.
    phasingGraph.writeGraphviz();
#endif
}



void Assembler::accessPhasingGraph()
{
    phasingGraph.orientedReads.accessExistingReadOnly(
        largeDataName("PhasingGraphOrientedReads"));
#if 0
    phasingGraph.assemblyGraphEdges.accessExistingReadOnly(
        largeDataName("PhasingGraphAssemblyGraphEdges"));
#endif
}



// Find the oriented reads internal to each assembly graph edge.
void Assembler::phasingGatherOrientedReads(size_t threadCount)
{
    phasingGraph.orientedReads.createNew(
        largeDataName("PhasingGraphOrientedReads"), largeDataPageSize);
    phasingGraph.orientedReads.beginPass1(assemblyGraph.edges.size());
    setupLoadBalancing(assemblyGraph.edges.size(), 1000);
    runThreads(&Assembler::phasingGatherOrientedReadsPass1, threadCount);
    phasingGraph.orientedReads.beginPass2();
    setupLoadBalancing(assemblyGraph.edges.size(), 1000);
    runThreads(&Assembler::phasingGatherOrientedReadsPass2, threadCount);
    phasingGraph.orientedReads.endPass2();
}



void Assembler::phasingGatherOrientedReadsPass1(size_t threadId)
{
    phasingGatherOrientedReadsPass(1);
}
void Assembler::phasingGatherOrientedReadsPass2(size_t threadId)
{
    phasingGatherOrientedReadsPass(2);
}
void Assembler::phasingGatherOrientedReadsPass(int pass)
{

    // Define this here to reduce memory allocation activity.
    vector<OrientedReadId> orientedReadIds;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while (getNextBatch(begin, end)) {

        // Loop over assembly graph edges assigned to this batch.
        for (AssemblyGraph::EdgeId assemblyGraphEdgeId = begin;
            assemblyGraphEdgeId != end; ++assemblyGraphEdgeId) {

            // Access the marker graph edges corresponding
            // to this assembly graph edge.
            const MemoryAsContainer<MarkerGraph::EdgeId> markerGraphEdges =
                assemblyGraph.edgeLists[assemblyGraphEdgeId];
            const uint64_t n = markerGraphEdges.size();
            SHASTA_ASSERT(n > 0);

            // Gather the oriented read ids.
            orientedReadIds.clear();
            if (n == 1) {

                // There is only one marker graph edge.
                // The OrientedReadId's are the ones in that one edge.
                const MarkerGraph::EdgeId markerGraphEdgeId =
                    markerGraphEdges[0];
                const MemoryAsContainer<MarkerInterval> markerIntervals =
                    markerGraph.edgeMarkerIntervals[markerGraphEdgeId];
                for (const MarkerInterval markerInterval : markerIntervals) {
                    orientedReadIds.push_back(markerInterval.orientedReadId);
                }

            } else {

                // The OrientedReadId's are the ones from all of the
                // marker graph vertices internal to this chain.
                // These are the target vertices of every edge except the last in the chain.
                for (uint64_t i = 0; i < n - 1; i++) {
                    const MarkerGraph::EdgeId markerGraphEdgeId =
                        markerGraphEdges[i];
                    const MarkerGraph::Edge &markedGraphEdge =
                        markerGraph.edges[markerGraphEdgeId];
                    const MarkerGraph::VertexId markerGraphVertexId =
                        markedGraphEdge.target;

                    // Loop over the markers in this marker graph vertex.
                    const MemoryAsContainer<MarkerId> markerIds =
                        markerGraph.vertices[markerGraphVertexId];
                    for (const MarkerId markerId : markerIds) {
                        OrientedReadId orientedReadId;
                        tie(orientedReadId, ignore) = findMarkerId(markerId);
                        orientedReadIds.push_back(orientedReadId);
                    }

                }
            }

            // Deduplicate.
            sort(orientedReadIds.begin(), orientedReadIds.end(),
                std::greater<OrientedReadId>());
            orientedReadIds.resize(
                unique(orientedReadIds.begin(), orientedReadIds.end())
                    - orientedReadIds.begin());

            // Store.
            if (pass == 1) {
                phasingGraph.orientedReads.incrementCount(assemblyGraphEdgeId,
                    orientedReadIds.size());
            } else {
                for (const OrientedReadId orientedReadId : orientedReadIds) {
                    phasingGraph.orientedReads.store(assemblyGraphEdgeId,
                        orientedReadId);
                }
            }

        }

    }

}


#if 0
// Find the assembly graph edges that each oriented read is internal to..
void Assembler::phasingGatherAssemblyGraphEdges(size_t threadCount)
{
    const uint64_t orientedReadCount = 2 * reads.size();

    phasingGraph.assemblyGraphEdges.createNew(
        largeDataName("PhasingGraphAssemblyGraphEdges"), largeDataPageSize);
    phasingGraph.assemblyGraphEdges.beginPass1(orientedReadCount);
    setupLoadBalancing(assemblyGraph.edges.size(), 1000);
    runThreads(&Assembler::phasingGatherAssemblyGraphEdgesPass1, threadCount);
    phasingGraph.assemblyGraphEdges.beginPass2();
    setupLoadBalancing(assemblyGraph.edges.size(), 1000);
    runThreads(&Assembler::phasingGatherAssemblyGraphEdgesPass2, threadCount);
    phasingGraph.assemblyGraphEdges.endPass2();
}

void Assembler::phasingGatherAssemblyGraphEdgesPass1(size_t threadId)
{
    phasingGatherAssemblyGraphEdgesPass(1);
}
void Assembler::phasingGatherAssemblyGraphEdgesPass2(size_t threadId)
{
    phasingGatherAssemblyGraphEdgesPass(2);
}



void Assembler::phasingGatherAssemblyGraphEdgesPass(int pass)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while (getNextBatch(begin, end)) {

        // Loop over assembly graph edges assigned to this batch.
        for (AssemblyGraph::EdgeId assemblyGraphEdgeId = begin;
            assemblyGraphEdgeId != end; ++assemblyGraphEdgeId) {

            // Access the oriented reads internal to this assembly graph edge.
            const MemoryAsContainer<OrientedReadId> orientedReadIds =
                phasingGraph.orientedReads[assemblyGraphEdgeId];

            // Loop over these oriented reads.
            for (const OrientedReadId orientedReadId : orientedReadIds) {
                if (pass == 1) {
                    phasingGraph.assemblyGraphEdges.incrementCountMultithreaded(
                        orientedReadId.getValue());
                } else {
                    phasingGraph.assemblyGraphEdges.storeMultithreaded(
                        orientedReadId.getValue(), assemblyGraphEdgeId);
                }
            }

        }

    }

}


void Assembler::phasingSortAssemblyGraphEdges(size_t threadCount)
{
    const uint64_t orientedReadCount = 2 * reads.size();
    setupLoadBalancing(orientedReadCount, 1000);
    runThreads(&Assembler::phasingSortAssemblyGraphEdgesThreadFunction,
        threadCount);

}
void Assembler::phasingSortAssemblyGraphEdgesThreadFunction(size_t threadId)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while (getNextBatch(begin, end)) {

        // Loop over oriented reads assigned to this batch.
        for (uint64_t orientedReadId = begin; orientedReadId != end;
            ++orientedReadId) {

            //  Sort the assembly graph edges that this oriented read is internal to.
            MemoryAsContainer<AssemblyGraph::EdgeId> edges =
                phasingGraph.assemblyGraphEdges[orientedReadId];
            sort(edges.begin(), edges.end());
        }
    }

}


double Assembler::computePhasingSimilarity(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1)
{
    return computePhasingSimilarity(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1)
        );
}
double Assembler::computePhasingSimilarity(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1)
{
    return phasingGraph.computePhasingSimilarity(orientedReadId0, orientedReadId1);
}
#endif



double Assembler::computePhasingSimilarity(
    AssemblyGraph::EdgeId edgeId0,
    AssemblyGraph::EdgeId edgeId1)
{
    return phasingGraph.computePhasingSimilarity(edgeId0, edgeId1);
}



#if 0
uint64_t Assembler::countCommonInternalOrientedReads(
    AssemblyGraph::EdgeId edgeId0,
    AssemblyGraph::EdgeId edgeId1)
{
    return phasingGraph.countCommonInternalOrientedReads(edgeId0, edgeId1);
}
#endif

