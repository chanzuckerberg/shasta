/*******************************************************************************

The read graph is a bidirected graph in which each vertex
represents a read.

Using a bidirected graph, we only need
one vertex per read, not two, like in the naive approach,
which uses one vertex for each of the two orientations
of a read.

For more information on the bidirected approach, see
Medvedev, Georgiou, Myers, and Brud,
"Computability of Models for Sequence Assembly",
International Workshop on Algorithms in Bioinformatics, pp 289-301 (2007).
http://medvedevgroup.com/papers/wabi07.pdf
In particular, see Fig. 1 and Section 2.2.

In summary:
- Each edge has two directions ("arrows"), one for each end.
- In a valid path, the two edges into and from a given vertex
  must have opposite orientations at that vertex.
- If the two edge orientations out of a vertex agree
  with the direction of the path,
  the read associated with the vertex is on strand 0 (unchanged).
  Otherwise, it is on strand 1 (reverse complemented).

We store edges with the lowest numbered read as the first read.

We use the standard approach to construct string graphs
(Myers, "The fragment assembly string graph" (2005),
doi:10.1093/bioinformatics/bti111,
http://www.cs.utoronto.ca/~brudno/csc2427/myers.pdf),
which prescribes that contained reads should not be included in the graph.
A contained read is a read that has an alignment
with a longer read covering the entire read.

However, in order not to lose coverage from contained reads
(which would result in a unacceptable reduction in coverage),
for each contained read we keep track of the containing read
that achieves the best alignment. This information is used later,
to incorporate contained reads in the marker graph.

*******************************************************************************/

// Shasta.
#include "Assembler.hpp"
#include "LocalReadGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Standard libraries.
#include "chrono.hpp"
#include <queue>



// Create the global read graph.
void Assembler::createReadGraph(uint32_t maxTrim)
{
    const bool debug = false;

    // Find the number of reads and oriented reads.
    const ReadId orientedReadCount = uint32_t(markers.size());
    CZI_ASSERT((orientedReadCount % 2) == 0);
    const ReadId readCount = orientedReadCount / 2;

    // Allocate containingOrientedReadId.
    containingOrientedReadId.createNew(largeDataName("ContainingOrientedReadId"), largeDataPageSize);
    containingOrientedReadId.resize(readCount);



    // Find which reads are contained.
    vector< pair<OrientedReadId, AlignmentInfo> > alignments;
    uint32_t containedCount = 0;
    for(ReadId readId0=0; readId0<readCount; readId0++) {

        // Get the information we need.
        const Strand strand0 = 0;
        const OrientedReadId orientedReadId0(readId0, strand0);
        const uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());

        // Get the alignments for this reads, properly oriented.
        alignments = findOrientedAlignments(orientedReadId0);

        // Check all the alignments, looking for the best containing oriented read.
        OrientedReadId bestContaining = OrientedReadId::invalid();
        uint32_t bestMarkerCount = 0;
        for(const auto& p: alignments) {
            const OrientedReadId orientedReadId1 = p.first;
            const AlignmentInfo& alignmentInfo = p.second;
            const uint32_t markerCount1 = uint32_t(markers[orientedReadId1.getValue()].size());

            // Compute the trim values.
            const uint32_t leftTrim0 = alignmentInfo.firstOrdinals.first;
            const uint32_t rightTrim0 = markerCount0 - 1 - alignmentInfo.lastOrdinals.first;
            const uint32_t leftTrim1 = alignmentInfo.firstOrdinals.second;
            const uint32_t rightTrim1 = markerCount1 - 1 - alignmentInfo.lastOrdinals.second;

            // Sanity check on the trim values.
            if( (leftTrim0 > maxTrim && rightTrim0 > maxTrim) &&
                (leftTrim1 > maxTrim && rightTrim1 > maxTrim)) {
                cout << "Found an invalid alignment:" << endl;
                cout << "orientedReadId0 " << orientedReadId0 << endl;
                cout << "orientedReadId1 " << orientedReadId1 << endl;
                cout << "markerCount0 " << markerCount0 << endl;
                cout << "markerCount1 " << markerCount1 << endl;
                cout << "leftTrim0 " << leftTrim0 << endl;
                cout << "rightTrim0 " << rightTrim0 << endl;
                cout << "leftTrim1 " << leftTrim1 << endl;
                cout << "rightTrim1 " << rightTrim1 << endl;
                throw runtime_error("Found a bad alignment.");
            }

            // If this is a containing alignment, see if it is better
            // that the one currently stored.
            if(leftTrim0<=maxTrim && rightTrim0<=maxTrim) {
                if(alignmentInfo.markerCount > bestMarkerCount) {
                    bestContaining = orientedReadId1;
                    bestMarkerCount = alignmentInfo.markerCount;
                }
            }
        }

        // Store the best containing oriented read.
        // If we did not find any, this will remain set to OrientedReadId::invalid().
        containingOrientedReadId[readId0] = bestContaining;
        if(bestContaining != OrientedReadId::invalid()) {
            ++containedCount;
        }
    }
    cout << "Found " << containedCount << " contained reads out of ";
    cout << readCount << " total." << endl;
    cout << "Number of non-contained reads is ";
    cout << readCount - containedCount << "." << endl;



    // Open the file for output of the global reead graph
    // (debugging only).
    ofstream graphOut;
    if(debug) {
        graphOut.open("ReadGraph.dot");
        graphOut << "graph G {" << endl;
    }



    // Create read graph edges.
    // Each alignment between non-contained reads generates an edge.
    readGraphEdges.createNew(largeDataName("ReadGraphEdges"), largeDataPageSize);
    for(size_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
        const AlignmentData& alignment = alignmentData[alignmentId];
        const AlignmentInfo& alignmentInfo = alignment.info;

        // If this alignment is not between contained reads, skip it.
        const ReadId readId0 = alignment.readIds[0];
        if(isContainedRead(readId0)) {
            continue;
        }
        const ReadId readId1 = alignment.readIds[1];
        if(isContainedRead(readId1)) {
            continue;
        }

        // Create the read graph edge.
        ReadGraphEdge edge;
        edge.alignmentId = alignmentId & 0x3fffffffffffffff;    // To suppress compiler warning

        // Compute the left and right trim for the two reads.
        const uint32_t markerCount0 = uint32_t(markers[OrientedReadId(readId0, 0).getValue()].size());
        const uint32_t markerCount1 = uint32_t(markers[OrientedReadId(readId1, 0).getValue()].size());
        CZI_ASSERT(alignmentInfo.lastOrdinals.first < markerCount0);
        CZI_ASSERT(alignmentInfo.lastOrdinals.second < markerCount1);
        const uint32_t leftTrim0 = alignmentInfo.firstOrdinals.first;
        const uint32_t rightTrim0 = markerCount0 - 1 - alignmentInfo.lastOrdinals.first;
        const uint32_t leftTrim1 = alignmentInfo.firstOrdinals.second;
        const uint32_t rightTrim1 = markerCount1 - 1 - alignmentInfo.lastOrdinals.second;

        // Sanity check on the left and right trim.
        if(
            (leftTrim0 > maxTrim && rightTrim0 > maxTrim) ||
            (leftTrim1 > maxTrim && rightTrim1 > maxTrim)) {

            cout << "Found a bad alignment:" << endl;
            cout << "Read ids: " << readId0 << " " << readId1 << endl;
            cout << "Is same strand: " << alignment.isSameStrand << endl;
            cout << "Trims for first read: " << leftTrim0 << " " << rightTrim0 << endl;
            cout << "Trims for second read: " << leftTrim1 << " " << rightTrim1 << endl;
            cout << "Marker counts: " << markerCount0 << " " << markerCount1 << endl;
            cout << "First ordinals: " << alignmentInfo.firstOrdinals.first << " " << alignmentInfo.firstOrdinals.second << endl;
            cout << "Last ordinals: " << alignmentInfo.lastOrdinals.first << " " << alignmentInfo.lastOrdinals.second << endl;
            throw runtime_error("Found a bad alignment.");
        }



        // Fill in the direction bits.
        // See comments at the beginning of this file for more information.
        // This code is untested.
        if(alignment.isSameStrand) {

            // The two reads are on the same strand.
            // Case 1 in Fig. 1 of the paper referenced above.
            // The "leftmost" read has direction away from the vertex
            // and the "rightmost" read has direction towards the vertex.
            if(leftTrim1 <= maxTrim) {
                edge.direction0 = ReadGraphEdge::awayFromVertex;
                edge.direction1 = ReadGraphEdge::towardsVertex;
            } else {
                edge.direction0 = ReadGraphEdge::towardsVertex;
                edge.direction1 = ReadGraphEdge::awayFromVertex;
            }

        } else {

            // The two reads are on opposite strands.
            if(leftTrim0 <= maxTrim) {
                // Case 2 in Fig. 1 of the paper referenced above.
                // Both directions are towards the vertex.
                edge.direction0 = ReadGraphEdge::towardsVertex;
                edge.direction1 = ReadGraphEdge::towardsVertex;

            } else {
                // Case 3 in Fig. 1 of the paper referenced above.
                // Both directions are away from the vertex.
                edge.direction0 = ReadGraphEdge::awayFromVertex;
                edge.direction1 = ReadGraphEdge::awayFromVertex;
            }

        }

        if(debug) {
            graphOut << readId0 << "--" << readId1 << ";\n";
        }

        // Store the new edge.
        readGraphEdges.push_back(edge);
    }

    if(debug) {
        graphOut << "}" << endl;
    }



    // Create read graph connectivity.
    readGraphConnectivity.createNew(largeDataName("ReadGraphConnectivity"), largeDataPageSize);
    readGraphConnectivity.beginPass1(readCount);
    for(const ReadGraphEdge& edge: readGraphEdges) {
        const AlignmentData& alignment = alignmentData[edge.alignmentId];
        readGraphConnectivity.incrementCount(alignment.readIds[0]);
        readGraphConnectivity.incrementCount(alignment.readIds[1]);
    }
    readGraphConnectivity.beginPass2();
    for(size_t i=0; i<readGraphEdges.size(); i++) {
        const ReadGraphEdge& edge = readGraphEdges[i];
        const AlignmentData& alignment = alignmentData[edge.alignmentId];
        readGraphConnectivity.store(alignment.readIds[0], uint32_t(i));
        readGraphConnectivity.store(alignment.readIds[1], uint32_t(i));
    }
    readGraphConnectivity.endPass2();
}



void Assembler::accessReadGraph()
{
    containingOrientedReadId.accessExistingReadOnly(largeDataName("ContainingOrientedReadId"));
    readGraphEdges.accessExistingReadOnly(largeDataName("ReadGraphEdges"));
    readGraphConnectivity.accessExistingReadOnly(largeDataName("ReadGraphConnectivity"));
}



// Follow the chain of containing reads until we reach a non-contained read.
OrientedReadId Assembler::findContainingReadRecursive(OrientedReadId orientedReadId) const
{
    while(true) {
        const ReadId readId = orientedReadId.getReadId();
        if(!isContainedRead(readId)) {
            return orientedReadId;
        }
        OrientedReadId containing = containingOrientedReadId[readId];
        if(orientedReadId.getStrand() == 1) {
            containing.flipStrand();
        }
        orientedReadId = containing;
    }
}



// Create a local subgraph of the global read graph,
// starting at a given vertex and extending out to a specified
// distance (number of edges).
// If the specified Readid corresponds to a contained read,
// which does not have a corresponding vertex in the read graph,
// the local subgraph starts instead from the containing read
// of the specified read.
bool Assembler::createLocalReadGraph(
    ReadId& readIdStart,            // If the specified read is contained, modified to the containing read.
    uint32_t maxDistance,           // How far to go from starting oriented read.
    double timeout,                 // Or 0 for no timeout.
    LocalReadGraph& graph)
{
    const auto startTime = steady_clock::now();

    // If the specified read is a contained read, find the containing read.
    if(isContainedRead(readIdStart)) {
        readIdStart = findContainingReadRecursive(OrientedReadId(readIdStart, 0)).getReadId();
    }

    // Add the starting vertex.
    graph.addVertex(readIdStart, uint32_t(markers[OrientedReadId(readIdStart, 0).getValue()].size()), 0);

    // Initialize a BFS starting at the start vertex.
    std::queue<ReadId> q;
    q.push(readIdStart);



    // Do the BFS.
    while(!q.empty()) {

        // See if we exceeded the timeout.
        if(seconds(steady_clock::now() - startTime) > timeout) {
            graph.clear();
            return false;
        }

        // Dequeue a vertex.
        const ReadId readId0 = q.front();
        q.pop();
        const uint32_t distance0 = graph.getDistance(readId0);
        const uint32_t distance1 = distance0 + 1;

        // Loop over edges of the global read graph involving this vertex.
        for(const uint64_t i: readGraphConnectivity[readId0]) {
            CZI_ASSERT(i < readGraphEdges.size());
            const ReadGraphEdge& globalEdge = readGraphEdges[i];
            const AlignmentData& alignment = alignmentData[globalEdge.alignmentId];

            // Get the other read involved in this edge of the read graph.
            const ReadId readId1 = alignment.getOther(readId0);


            // Update our BFS.
            // Note that we are pushing to the queue vertices at maxDistance,
            // so we can find all of their edges to other vertices at maxDistance.
            if(distance0 < maxDistance) {
                if(!graph.vertexExists(readId1)) {
                    graph.addVertex(readId1,
                        uint32_t(markers[OrientedReadId(readId1, 0).getValue()].size()), distance1);
                    q.push(readId1);
                }
                graph.addEdge(
                    alignment.readIds[0],
                    alignment.readIds[1],
                    globalEdge.alignmentId,
                    globalEdge.direction0,
                    globalEdge.direction1);
            } else {
                CZI_ASSERT(distance0 == maxDistance);
                if(graph.vertexExists(readId1)) {
                    graph.addEdge(
                        alignment.readIds[0],
                        alignment.readIds[1],
                        globalEdge.alignmentId,
                        globalEdge.direction0,
                        globalEdge.direction1);
                }
            }

        }

    }
    return true;
}
