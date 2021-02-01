#include "Assembler.hpp"
#include "LocalAlignmentCandidateGraph.hpp"
using namespace shasta;

#include "chrono.hpp"
#include <queue>




bool Assembler::createLocalAlignmentCandidateGraph(
        vector<OrientedReadId>& starts,
        uint32_t maxDistance,           // How far to go from starting oriented read.
        bool allowChimericReads,
        double timeout,                 // Or 0 for no timeout.
        LocalAlignmentCandidateGraph& graph)
{
    const auto startTime = steady_clock::now();

    // Initialize a BFS starting at the start vertex.
    std::queue<OrientedReadId> q;

    for (auto& start: starts) {
        // If the starting read is chimeric and we don't allow chimeric reads, do nothing.
        if (!allowChimericReads && reads->getFlags(start.getReadId()).isChimeric) {
            continue;
        }

        // Add the starting vertex.
        graph.addVertex(start, uint32_t(reads->getRead(start.getReadId()).baseCount), 0);

        // Add each starting vertex to the BFS queue
        q.push(start);
    }

    // Do the BFS.
    while(!q.empty()) {

        // See if we exceeded the timeout.
        if(seconds(steady_clock::now() - startTime) > timeout) {
            graph.clear();
            return false;
        }

        // Dequeue a vertex.
        const OrientedReadId orientedReadId0 = q.front();
        // cout << "Dequeued " << orientedReadId0;
        // cout << " with " << candidateTable.size(orientedReadId0.getValue()) << " overlaps." << endl;
        q.pop();
        const uint32_t distance0 = graph.getDistance(orientedReadId0);
        const uint32_t distance1 = distance0 + 1;

        // Loop over overlaps involving this vertex.
        for(const uint64_t i: alignmentCandidates.candidateTable[orientedReadId0.getValue()]) {
            const OrientedReadPair& pair = alignmentCandidates.candidates[i];

            // Get the other oriented read involved in this overlap.
            const OrientedReadId orientedReadId1 = pair.getOther(orientedReadId0);

            bool inAlignments = false;
            bool inReadgraph = false;
            bool inReferenceAlignments = false;

            // Search the AlignmentTable to see if this pair exists
            for(const ReadId alignmentIndex: alignmentTable[orientedReadId0.getValue()]) {
                const AlignmentData& ad = alignmentData[alignmentIndex];

                // Check if the pair matches the current candidate pair of interest
                 if (ad.getOther(orientedReadId0) == orientedReadId1){
                     inAlignments = true;
                 }
            }

            // Search the ReadGraph to see if this pair exists
            for(const ReadId readGraphIndex: readGraph.connectivity[orientedReadId0.getValue()]) {
                const ReadGraphEdge& e = readGraph.edges[readGraphIndex];

                // Check if the pair matches the current candidate pair of interest
                if (e.getOther(orientedReadId0) == orientedReadId1){
                    inReadgraph = true;
                }
            }

            // TODO add reference alignment table and check for this edge's existence


            // Update our BFS.
            // Note that we are pushing to the queue vertices at maxDistance,
            // so we can find all of their edges to other vertices at maxDistance.
            if(distance0 < maxDistance) {
                if(!graph.vertexExists(orientedReadId1)) {
                    graph.addVertex(orientedReadId1,
                                    uint32_t(reads->getRead(orientedReadId1.getReadId()).baseCount), distance1);
                    q.push(orientedReadId1);
                }

                graph.addEdge(orientedReadId0, orientedReadId1, inAlignments, inReadgraph, inReferenceAlignments);
            } else {
                SHASTA_ASSERT(distance0 == maxDistance);
                if(graph.vertexExists(orientedReadId1)) {
                    graph.addEdge(orientedReadId0, orientedReadId1, inAlignments, inReadgraph, inReferenceAlignments);
                }
            }
        }
    }

    return true;
}


// Compute candidateTable from alignmentCandidates.
// This could be made multithreaded if it becomes a bottleneck.
void AlignmentCandidates::computeCandidateTable(ReadId readCount, string largeDataName, size_t largeDataPageSize)
{
    candidateTable.createNew(largeDataName, largeDataPageSize);
    candidateTable.beginPass1(ReadId(2 * readCount));
    for(const OrientedReadPair& pair: candidates) {
        OrientedReadId orientedReadId0(pair.readIds[0], 0);
        OrientedReadId orientedReadId1(pair.readIds[1], pair.isSameStrand ? 0 : 1);
        candidateTable.incrementCount(orientedReadId0.getValue());
        candidateTable.incrementCount(orientedReadId1.getValue());
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        candidateTable.incrementCount(orientedReadId0.getValue());
        candidateTable.incrementCount(orientedReadId1.getValue());
    }
    candidateTable.beginPass2();
    for(uint32_t i=0; i<candidates.size(); i++) {
        const OrientedReadPair& pair = candidates[i];
        OrientedReadId orientedReadId0(pair.readIds[0], 0);
        OrientedReadId orientedReadId1(pair.readIds[1], pair.isSameStrand ? 0 : 1);
        candidateTable.store(orientedReadId0.getValue(), i);
        candidateTable.store(orientedReadId1.getValue(), i);
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        candidateTable.store(orientedReadId0.getValue(), i);
        candidateTable.store(orientedReadId1.getValue(), i);
    }
    candidateTable.endPass2();



    // Sort each section of the candidate table by OrientedReadId.
    vector< pair<OrientedReadId, uint32_t> > v;
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        for(Strand strand0=0; strand0<2; strand0++) {
            const OrientedReadId orientedReadId0(readId0, strand0);

            // Access the section of the candidate table for this oriented read.
            const span<uint32_t> candidateTableSection =
                    candidateTable[orientedReadId0.getValue()];

            // Store pairs(OrientedReadId, candidateIndex).
            v.clear();
            for(uint32_t candidateIndex: candidateTableSection) {
                const OrientedReadPair& pair = candidates[candidateIndex];
                const OrientedReadId orientedReadId1 = pair.getOther(orientedReadId0);
                v.push_back(make_pair(orientedReadId1, candidateIndex));
            }

            // Sort.
            sort(v.begin(), v.end());

            // Store the sorted candidateIndex.
            for(size_t i=0; i<v.size(); i++) {
                candidateTableSection[i] = v[i].second;
            }
        }
    }

    candidateTable.unreserve();

}
