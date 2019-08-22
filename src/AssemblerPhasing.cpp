#include "Assembler.hpp"
using namespace shasta;


void Assembler::createPhasingGraph(size_t threadCount)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Find the oriented reads internal to each assembly graph edge.
    phasingGraph.orientedReads.createNew(
		largeDataName("PagingGraphOrientedReads"), largeDataPageSize);
	phasingGraph.orientedReads.beginPass1(assemblyGraph.edges.size());
    setupLoadBalancing(assemblyGraph.edges.size(), 1000);
    runThreads(&Assembler::createPhasingGraphGatherOrientedReadsPass1, threadCount);
	phasingGraph.orientedReads.beginPass2();
    setupLoadBalancing(assemblyGraph.edges.size(), 1000);
    runThreads(&Assembler::createPhasingGraphGatherOrientedReadsPass2, threadCount);
	phasingGraph.orientedReads.endPass2();
}



void Assembler::createPhasingGraphGatherOrientedReadsPass1(size_t threadId)
{
	createPhasingGraphGatherOrientedReads(1);
}
void Assembler::createPhasingGraphGatherOrientedReadsPass2(size_t threadId)
{
	createPhasingGraphGatherOrientedReads(2);
}
void Assembler::createPhasingGraphGatherOrientedReads(int pass)
{

	// Define this here to reduce memory allocation activity.
	vector<OrientedReadId> orientedReadIds;

	// Loop over all batches assigned to this thread.
	uint64_t begin, end;
	while(getNextBatch(begin, end)) {

		// Loop over assembly graph edges assigned to this batch.
		for(AssemblyGraph::EdgeId assemblyGraphEdgeId=begin;
				assemblyGraphEdgeId!=end; ++assemblyGraphEdgeId) {

			// Access the marker graph edges corresponding
			// to this assembly graph edge.
			const MemoryAsContainer<MarkerGraph::EdgeId> markerGraphEdges =
				assemblyGraph.edgeLists[assemblyGraphEdgeId];
			const uint64_t n = markerGraphEdges.size();
			SHASTA_ASSERT(n > 0);


			// Gather the oriented read ids.
			orientedReadIds.clear();
			if(n == 1) {

				// There is only one marker graph edge.
				// The OrientedReadId's are the ones in that one edge.
				const MarkerGraph::EdgeId markerGraphEdgeId = markerGraphEdges[0];
				const MemoryAsContainer<MarkerInterval> markerIntervals =
					markerGraph.edgeMarkerIntervals[markerGraphEdgeId];
				for(const MarkerInterval markerInterval: markerIntervals) {
					orientedReadIds.push_back(markerInterval.orientedReadId);
				}

			} else {

				// The OrientedReadId's are the ones from all of the
				// marker graph vertices internal to this chain.
				// These are the target vertices of every edge except the last in the chain.
				for(uint64_t i=0; i<n-1; i++) {
					const MarkerGraph::EdgeId markerGraphEdgeId = markerGraphEdges[i];
					const MarkerGraph::Edge& markedGraphEdge = markerGraph.edges[markerGraphEdgeId];
					const MarkerGraph::VertexId markerGraphVertexId = markedGraphEdge.target;

					// Loop over the markers in this marker graph vertex.
					const MemoryAsContainer<MarkerId> markerIds = markerGraph.vertices[markerGraphVertexId];
					for(const MarkerId markerId: markerIds) {
						OrientedReadId orientedReadId;
						tie(orientedReadId, ignore)  = findMarkerId(markerId);
						orientedReadIds.push_back(orientedReadId);
					}

				}
			}

			// Deduplicate.
			sort(orientedReadIds.begin(), orientedReadIds.end(), std::greater<OrientedReadId>());
			orientedReadIds.resize(
				unique(orientedReadIds.begin(), orientedReadIds.end()) - orientedReadIds.begin());

			// Store.
			if(pass == 1) {
				phasingGraph.orientedReads.incrementCount(assemblyGraphEdgeId, orientedReadIds.size());
			} else {
				for(const OrientedReadId orientedReadId: orientedReadIds) {
					phasingGraph.orientedReads.store(assemblyGraphEdgeId, orientedReadId);
				}
			}

		}

	}

}
