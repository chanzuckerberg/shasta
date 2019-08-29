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
    phasingGraph.dataFileNamePrefix = largeDataFileNamePrefix;
    phasingGraph.dataPageSize = largeDataPageSize;

    // Find the oriented reads internal to each assembly graph edge.
    phasingGatherOrientedReads(threadCount);

    // Find the assembly graph edges that each oriented read is internal to.
    phasingGatherAssemblyGraphEdges(threadCount);
    phasingSortAssemblyGraphEdges(threadCount);

    // Find oriented read pairs with phasing similarity greater than the threshold.
    phasingGraph.findSimilarPairs(threadCount, phasingSimilarityThreshold);

    // Only keep up to maxNeighborCount neighbors.
    phasingGraph.keepBestSimilarPairs(maxNeighborCount);

    // Write out the global phasing graph in graphviz format.
    phasingGraph.writeGraphviz();
}



void Assembler::accessPhasingGraph()
{
    phasingGraph.orientedReads.accessExistingReadOnly(
        largeDataName("PhasingGraphOrientedReads"));
    phasingGraph.assemblyGraphEdges.accessExistingReadOnly(
        largeDataName("PhasingGraphAssemblyGraphEdges"));
#if 0
    phasingGraph.turns.accessExistingReadOnly(
        largeDataName("PhasingGraphTurns"));
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
double Assembler::computePhasingSimilarity(
    AssemblyGraph::EdgeId edgeId0,
    AssemblyGraph::EdgeId edgeId1)
{
    return phasingGraph.computePhasingSimilarity(edgeId0, edgeId1);
}
uint64_t Assembler::countCommonInternalOrientedReads(
    AssemblyGraph::EdgeId edgeId0,
    AssemblyGraph::EdgeId edgeId1)
{
    return phasingGraph.countCommonInternalOrientedReads(edgeId0, edgeId1);
}


#if 0
pair<uint64_t, uint64_t> Assembler::countCommonTurns(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1)
{
    return countCommonTurns(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1)
        );
}
pair<uint64_t, uint64_t> Assembler::countCommonTurns(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1)
{
    return phasingGraph.countCommonTurns(orientedReadId0, orientedReadId1);
}



// Find turns in the phasing graph.
// See the definition of class PhasingGraph::Turn for more information.
void Assembler::phasingFindTurns()
{
    // Initialize the data structure to contain turns for each oriented read.
    phasingGraph.turns.createNew(
        largeDataName("PhasingGraphTurns"), largeDataPageSize);

    // Work areas defined here to reduce memnory allocation activity.
    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::EdgeId;
    vector<EdgeId> sources;
    vector<EdgeId> targets;
    vector<VertexId> candidateHinges;

    // Loop over all oriented reads.
    const ReadId readCount = ReadId(reads.size());
    for(ReadId readId=0; readId!=readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);

            // Add a new vector of turns for this oriented reads.
            phasingGraph.turns.appendVector();

            // Access the assembly graph edges that this oriented read
            // is internal to.
            MemoryAsContainer<EdgeId> edgeIds =
                phasingGraph.assemblyGraphEdges[orientedReadId.getValue()];

            // Each turn vertex must be the target of one of these
            // edges and also the source of one of these edges.
            // Therefore, to find possible turn vertices
            // we compute the intersection of the
            // sources and targets of these edges.
            sources.clear();
            targets.clear();
            for(const EdgeId edgeId: edgeIds) {
                const AssemblyGraph::Edge edge = assemblyGraph.edges[edgeId];
                sources.push_back(edge.source);
                targets.push_back(edge.target);
            }
            // Sort and deduplicate.
            sort(sources.begin(), sources.end());
            sort(targets.begin(), targets.end());
            sources.resize(unique(sources.begin(), sources.end()) - sources.begin());
            targets.resize(unique(targets.begin(), targets.end()) - targets.begin());
            // Compute the intersection.
            candidateHinges.clear();
            std::set_intersection(
                sources.begin(), sources.end(),
                targets.begin(), targets.end(),
                back_inserter(candidateHinges));

            // Loop over all candidate hinges.
            // The candidate hinges are sorted, so the turns are generated
            // in order of increasing v.
            for(VertexId v: candidateHinges) {

                // Find the one and only incoming edge
                // that contains our oriented read. If more than one found,
                // this is not a valid hinge because of condition 3.
                // in the definition of a Turn.
                EdgeId e0 = AssemblyGraph::invalidEdgeId;
                bool isHinge = true;
                const MemoryAsContainer<EdgeId> incomingEdges =
                    assemblyGraph.edgesByTarget[v];
                for(const EdgeId edgeId: incomingEdges) {
                    if(std::binary_search(edgeIds.begin(), edgeIds.end(), edgeId)) {
                        // This edge contains our oriented read.
                        if(e0 == AssemblyGraph::invalidEdgeId) {
                            e0 = edgeId;
                        } else {
                            // This is the second one we found.
                            isHinge = false;
                            break;
                        }
                    }
                }
                if(!isHinge) {
                    continue;
                }

                // Find the one and only outgoing edge
                // that contains our oriented read. If more than one found,
                // this is not a valid hinge because of condition 4.
                // in the definition of a Turn.
                EdgeId e1 = AssemblyGraph::invalidEdgeId;
                const MemoryAsContainer<EdgeId> outgoingEdges =
                    assemblyGraph.edgesBySource[v];
                for(const EdgeId edgeId: outgoingEdges) {
                    if(std::binary_search(edgeIds.begin(), edgeIds.end(), edgeId)) {
                        // This edge contains our oriented read.
                        if(e1 == AssemblyGraph::invalidEdgeId) {
                            e1 = edgeId;
                        } else {
                            // This is the second one we found.
                            isHinge = false;
                            break;
                        }
                    }
                }
                if(!isHinge) {
                    continue;
                }

                phasingGraph.turns.append(PhasingGraph::Turn(e0, e1, v));
                /*
                cout << orientedReadId << " ";
                cout << e0 << " ";
                cout << e1 << " ";
                cout << v << "\n";
                */
            }
        }
    }
}
#endif

