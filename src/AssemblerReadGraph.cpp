// shasta.
#include "Assembler.hpp"
#include "LocalReadGraph.hpp"
#include "orderPairs.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard libraries.
#include "chrono.hpp"
#include <queue>
#include "tuple.hpp"



// Create a local read graph starting from a given oriented read
// and walking out a given distance on the global read graph.
bool Assembler::createLocalReadGraph(
    OrientedReadId orientedReadIdStart,
    size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
    size_t maxTrim,                 // Maximum left/right trim (expressed in bases) to generate an edge.
    uint32_t maxDistance,           // How far to go from starting oriented read.
    double timeout,                 // Or 0 for no timeout.
    LocalReadGraph& graph)
{
    const auto startTime = steady_clock::now();

    // Add the starting vertex.
    graph.addVertex(orientedReadIdStart,
        uint32_t(reads[orientedReadIdStart.getReadId()].baseCount), 0);

    // Initialize a BFS starting at the start vertex.
    std::queue<OrientedReadId> q;
    q.push(orientedReadIdStart);



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
        // cout << " with " << alignmentTable.size(orientedReadId0.getValue()) << " overlaps." << endl;
        q.pop();
        const uint32_t distance0 = graph.getDistance(orientedReadId0);
        const uint32_t distance1 = distance0 + 1;

        // Loop over overlaps/alignments involving this vertex.
        for(const uint64_t i: alignmentTable[orientedReadId0.getValue()]) {
            CZI_ASSERT(i < alignmentData.size());
            const AlignmentData& ad = alignmentData[i];

            // If the alignment involves too few markers, skip.
            if(ad.info.markerCount < minAlignedMarkerCount) {
                continue;
            }

            // To compute the trim, keep into account the fact
            // that the stored AlignmentInfo was computed for
            // the ReadId's stored in the Overlap, with the first one on strand 0.
            const OrientedReadId overlapOrientedReadId0(ad.readIds[0], 0);
            const OrientedReadId overlapOrientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
            uint32_t leftTrim;
            uint32_t rightTrim;
            tie(leftTrim, rightTrim) = computeTrim(
                overlapOrientedReadId0,
                overlapOrientedReadId1,
                ad.info);
            if(leftTrim>maxTrim || rightTrim>maxTrim) {
                continue;
            }

            // The overlap and the alignment satisfy our criteria.
            // Get the other oriented read involved in this overlap.
            const OrientedReadId orientedReadId1 = ad.getOther(orientedReadId0);


            // Update our BFS.
            // Note that we are pushing to the queue vertices at maxDistance,
            // so we can find all of their edges to other vertices at maxDistance.
            if(distance0 < maxDistance) {
                if(!graph.vertexExists(orientedReadId1)) {
                    graph.addVertex(orientedReadId1,
                        uint32_t(reads[orientedReadId1.getReadId()].baseCount), distance1);
                    q.push(orientedReadId1);
                }
                graph.addEdge(orientedReadId0, orientedReadId1,
                    ad.info);
            } else {
                CZI_ASSERT(distance0 == maxDistance);
                if(graph.vertexExists(orientedReadId1)) {
                    graph.addEdge(orientedReadId0, orientedReadId1,
                        ad.info);
                }
            }


#if 0
            // THE OLD CODE DOES NOT CREATE EDGES
            // BETWEEN VERTICES AT MAXDISTANCE.
            // Add the vertex for orientedReadId1, if necessary.
            // Also add orientedReadId1 to the queue, unless
            // we already reached the maximum distance.
            if(!graph.vertexExists(orientedReadId1)) {
                graph.addVertex(orientedReadId1,
                    uint32_t(reads[orientedReadId1.getReadId()].baseCount), distance1);
                if(distance1 < maxDistance) {
                    q.push(orientedReadId1);
                    // cout << "Enqueued " << orientedReadId1 << endl;
                }
            }

            // Add the edge.
            graph.addEdge(orientedReadId0, orientedReadId1,
                ad.info);
#endif
        }

    }

    return true;
}



// Create the global read graph.
void Assembler::createReadGraph(uint32_t maxTrim)
{
    // Find the number of reads and oriented reads.
    const ReadId orientedReadCount = uint32_t(markers.size());
    CZI_ASSERT((orientedReadCount % 2) == 0);
    const ReadId readCount = orientedReadCount / 2;



    // Find which oriented reads are contained.
    vector<bool> isContained(orientedReadCount);
    uint32_t containedCount = 0;
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        for(Strand strand0=0; strand0<2; strand0++) {

            // Get the information we need, including the alignments for this oriented read.
            const OrientedReadId orientedReadId0(readId0, strand0);
            const vector< pair<OrientedReadId, AlignmentInfo> > alignments =
                findOrientedAlignments(orientedReadId0);
            const uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());

            // Check all the alignments.
            isContained[orientedReadId0.getValue()] = false;
            for(const auto& p: alignments) {
                const AlignmentInfo& alignmentInfo = p.second;
                const uint32_t leftTrim0 = alignmentInfo.firstOrdinals.first;
                const uint32_t rightTrim0 = markerCount0 - 1 - alignmentInfo.lastOrdinals.first;
                if(leftTrim0<=maxTrim && rightTrim0<=maxTrim) {
                    isContained[orientedReadId0.getValue()] = true;
                    ++containedCount;
                    break;
                }
            }
        }
    }

    cout << "Found " << containedCount << " contained oriented reads out of ";
    cout << orientedReadCount << " total." << endl;
    cout << "Number of non-contained oriented reads is ";
    cout << orientedReadCount - containedCount << "." << endl;


    ofstream graphOut("ReadGraph.dot");
    graphOut << "graph G {\n";
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        for(Strand strand0=0; strand0<2; strand0++) {
            const OrientedReadId orientedReadId0(readId0, strand0);
            if(isContained[orientedReadId0.getValue()]) {
                continue;
            }

            // Get the information we need, including the alignments for this oriented read.
            const vector< pair<OrientedReadId, AlignmentInfo> > alignments =
                findOrientedAlignments(orientedReadId0);
            const uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());

            // Check all the alignments.
            for(const auto& p: alignments) {
                const OrientedReadId orientedReadId1 = p.first;
                if(isContained[orientedReadId1.getValue()]) {
                    continue;
                }
                if(orientedReadId1 < orientedReadId0) {
                    continue;
                }

                const uint32_t markerCount1 = uint32_t(markers[orientedReadId0.getValue()].size());
                const AlignmentInfo& alignmentInfo = p.second;

                // If there is too much trim on both sides of orientedReadId0,
                // this is is a suspicious alignment. Discard.
                const uint32_t leftTrim0 = alignmentInfo.firstOrdinals.first;
                const uint32_t rightTrim0 = markerCount0 - 1 - alignmentInfo.lastOrdinals.first;
                if(leftTrim0 > maxTrim && rightTrim0 > maxTrim) {
                    continue;
                }

                // If there is too much trim on both sides of orientedReadId1,
                // this is is a suspicious alignment. Discard.
                const uint32_t leftTrim1 = alignmentInfo.firstOrdinals.second;
                const uint32_t rightTrim1 = markerCount1 - 1 - alignmentInfo.lastOrdinals.second;
                if(leftTrim1 > maxTrim && rightTrim1 > maxTrim) {
                    continue;
                }

                graphOut << orientedReadId0.getValue() << "--";
                graphOut << orientedReadId1.getValue() << ";\n";
            }
        }
    }
    graphOut << "}\n";
}
