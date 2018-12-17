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

The read graph uses the best maxAlignmentCount alignments for each read.
Additional alignments do not generate an edge of the read graph.

*******************************************************************************/

// Shasta.
#include "Assembler.hpp"
#include "LocalReadGraph.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard libraries.
#include "chrono.hpp"
#include <queue>


#if 0
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
#endif



// For each read, keep only the best maxAlignmentCount alignments.
// Note that the connectivity of the resulting read graph can
// be more than maxAlignmentCount.
void Assembler::createReadGraph(uint32_t maxTrim)
{
    // Number of alignments to be kept for each read.
    // This should be passed in as an argument instead.
    const size_t maxAlignmentCount = 6;

    // Find the number of reads and oriented reads.
    const ReadId orientedReadCount = uint32_t(markers.size());
    CZI_ASSERT((orientedReadCount % 2) == 0);
    const ReadId readCount = orientedReadCount / 2;

    // Mark all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentData.size(), false);

    // Vector to keep the alignments for each read,
    // with their number of markers.
    // Contains pairs(marker count, alignment id).
    vector< pair<uint32_t, uint32_t> > readAlignments;



    // Loop over reads.
    for(ReadId readId=0; readId<readCount; readId++) {

        // Gather the alignments for this read, each with its number of markers.
        readAlignments.clear();
        for(const uint32_t alignmentId: alignmentTable[OrientedReadId(readId, 0).getValue()]) {
            const AlignmentData& alignment = alignmentData[alignmentId];
            readAlignments.push_back(make_pair(alignment.info.markerCount, alignmentId));
        }

        // Keep the best maxAlignmentCount.
        if(readAlignments.size() > maxAlignmentCount) {
            std::nth_element(
                readAlignments.begin(),
                readAlignments.begin() + maxAlignmentCount,
                readAlignments.end(),
                std::greater< pair<uint32_t, uint32_t> >());
            readAlignments.resize(maxAlignmentCount);
        }

        // Mark the surviving alignments as to be kept.
        for(const auto& p: readAlignments) {
            const uint32_t alignmentId = p.second;
            keepAlignment[alignmentId] = true;
        }
    }
    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;



    // Now we can create the read graph.
    // Only the alignments we marker as "keep" generate an edge in the read graph.
    readGraphEdges.createNew(largeDataName("ReadGraphEdges"), largeDataPageSize);
    for(size_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
        if(!keepAlignment[alignmentId]) {
            continue;
        }
        const AlignmentData& alignment = alignmentData[alignmentId];
        const AlignmentInfo& alignmentInfo = alignment.info;

        // Create the read graph edge.
        ReadGraphEdge edge;
        edge.alignmentId = alignmentId & 0x3fffffffffffffff;    // To suppress compiler warning

        // Compute the left and right trim for the two reads.
        const ReadId readId0 = alignment.readIds[0];
        const ReadId readId1 = alignment.readIds[1];
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
            (leftTrim0 > maxTrim && rightTrim0 > maxTrim) &&
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

        // Store the new edge.
        readGraphEdges.push_back(edge);
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
    readGraphEdges.accessExistingReadOnly(largeDataName("ReadGraphEdges"));
    readGraphConnectivity.accessExistingReadOnly(largeDataName("ReadGraphConnectivity"));
}
void Assembler::checkReadGraphIsOpen()
{
    if(!readGraphEdges.isOpen) {
        throw runtime_error("Read graph edges are not accessible.");
    }
    if(!readGraphConnectivity.isOpen()) {
        throw runtime_error("Read graph connectivity is not accessible.");
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
    uint32_t maxTrim,               // To define alignment containment
    bool allowChimericReads,
    double timeout,                 // Or 0 for no timeout.
    LocalReadGraph& graph)
{
    const auto startTime = steady_clock::now();

    // If the starting read is chimeric and we don't allow chimeric reads, do nothing.
    if(!allowChimericReads && isChimericRead[readIdStart]) {
        return true;
    }

    // Add the starting vertex.
    graph.addVertex(readIdStart, uint32_t(markers[OrientedReadId(readIdStart, 0).getValue()].size()),
        isChimericRead[readIdStart], 0);

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

            // If this read is flagged chimeric and we don't allow chimeric reads, skip.
            if(!allowChimericReads && isChimericRead[readId1]) {
                continue;
            }


            // Update our BFS.
            // Note that we are pushing to the queue vertices at maxDistance,
            // so we can find all of their edges to other vertices at maxDistance.
            if(distance0 < maxDistance) {
                if(!graph.vertexExists(readId1)) {
                    graph.addVertex(readId1,
                        uint32_t(markers[OrientedReadId(readId1, 0).getValue()].size()),
                        isChimericRead[readId1], distance1);
                    q.push(readId1);
                }
                graph.addEdge(
                    alignment.readIds[0],
                    alignment.readIds[1],
                    globalEdge.alignmentId,
                    globalEdge.direction0,
                    globalEdge.direction1,
                    isContainmentAlignment(globalEdge.alignmentId, maxTrim));
            } else {
                CZI_ASSERT(distance0 == maxDistance);
                if(graph.vertexExists(readId1)) {
                    graph.addEdge(
                        alignment.readIds[0],
                        alignment.readIds[1],
                        globalEdge.alignmentId,
                        globalEdge.direction0,
                        globalEdge.direction1,
                        isContainmentAlignment(globalEdge.alignmentId, maxTrim));
                }
            }

        }

    }
    return true;
}



// Use the read graph to flag chimeric reads.
// For each read and corresponding vertex v0, we do
// a BFS in the read graph up to the specified maxDistance.
// We then compute connected components of the subgraph
// consisting of the vertices reached by the bfs, minus v0.
// If not all the vertices at maximum distance are
// in the same component, the read corresponding to v0
// is flagged as chimeric.
void Assembler::flagChimericReads(size_t maxDistance, size_t threadCount)
{
    cout << timestamp << "Begin flagging chimeric reads, max distance " << maxDistance << endl;

    // Check that we have what we need.
    checkReadGraphIsOpen();
    const size_t readCount = readGraphConnectivity.size();

    // Create the flags to indicate which reads are chimeric.
    isChimericRead.createNew(
        largeDataName("IsChimericRead"),
        largeDataPageSize);
    isChimericRead.resize(readCount);

    // Store the argument so it is accessible by all threads.
    CZI_ASSERT(maxDistance < 255);
    flagChimericReadsData.maxDistance = maxDistance;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;
    // Multithreaded loop over all reads.
    cout << timestamp << "Processing " << readCount << " reads." << endl;
    setupLoadBalancing(readCount, 10000);
    runThreads(&Assembler::flagChimericReadsThreadFunction, threadCount,
        "threadLogs/flagChimericReads");

    cout << timestamp << "Done flagging chimeric reads." << endl;

    const size_t chimericReadCount = std::count(isChimericRead.begin(), isChimericRead.end(), true);
    cout << timestamp << "Flagged " << chimericReadCount << " reads as chimeric out of ";
    cout << readCount << " total." << endl;
    cout << "Chimera rate is " << double(chimericReadCount) / double(readCount) << endl;
}



void Assembler::flagChimericReadsThreadFunction(size_t threadId)
{
    ostream& out = getLog(threadId);
    const size_t maxDistance = flagChimericReadsData.maxDistance;

    // Vector used for BFS searches by this thread.
    // It stores the local vertex id in the current BFS assigned to each vertex,
    // or notReached for vertices not yet reached by the current BFS.
    // This is of size equal to the number of reads, and each thread has its own copy.
    // This is not prohibitive. For example, for a large human size run with
    // 20 million reads and 100 threads, the total space is only 8 GB.
    const ReadId notReached = std::numeric_limits<ReadId>::max();
    MemoryMapped::Vector<ReadId> vertexTable;
    vertexTable.createNew(
        largeDataName("tmp-FlagChimericReads-VertexTable" + to_string(threadId)),
        largeDataPageSize);
    vertexTable.resize(readGraphConnectivity.size());
    fill(vertexTable.begin(), vertexTable.end(), notReached);

    // Vector to contain the vertices we found in the current BFS,
    // each with the distance from the start vertex.
    vector< pair<ReadId, uint32_t> > localVertices;

    // The queue used for the BFS.
    std::queue<ReadId> q;

    // Vectors used to compute connected components after each BFS>
    vector<ReadId> rank;
    vector<ReadId> parent;


    // Loop over all batches assigned to this thread.
    size_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << "Working on batch " << begin << " " << end << endl;

        // Loop over all reads assigned to this batch.
        for(ReadId vStart=ReadId(begin); vStart!=ReadId(end); vStart++) {
            // out<< timestamp << "Working on read " << vStart << endl;

            // Check that there is no garbage left by the previous BFS.
            CZI_ASSERT(localVertices.empty());
            CZI_ASSERT(q.empty());

            // Begin by flagging this read as not chimeric.
            isChimericRead[vStart] = false;



            // Do the BFS for this read.
            ReadId localVertexId = 0;
            q.push(vStart);
            localVertices.push_back(make_pair(vStart, 0));
            vertexTable[vStart] = localVertexId++;
            while(!q.empty()) {

                // Dequeue a vertex.
                const ReadId v0 = q.front();
                q.pop();
                const uint32_t distance0 = localVertices[vertexTable[v0]].second;
                const uint32_t distance1 = distance0 + 1;
                // out << "Dequeued " << v0 << endl;

                // Loop over edges involving this vertex.
                const auto edges = readGraphConnectivity[v0];
                for(const uint32_t edgeId: edges) {
                    const ReadGraphEdge& edge = readGraphEdges[edgeId];
                    const AlignmentData& alignment = alignmentData[edge.alignmentId];
                    const ReadId v1 = alignment.getOther(v0);
                    // out << "Found " << v1 << endl;

                    // If we already encountered this read in this BFS, do nothing.
                    if(vertexTable[v1] != notReached) {
                        // out << "Previously reached." << endl;
                        continue;
                    }

                    // Record this vertex.
                    // out << "Recording " << v1 << endl;
                    localVertices.push_back(make_pair(v1, distance1));
                    vertexTable[v1] = localVertexId++;

                    // If at distance less than maxDistance, also enqueue it.
                    if(distance1 < maxDistance) {
                        // out << "Enqueueing " << v1 << endl;
                        q.push(v1);
                    }
                }
            }
            // out << "BFS found " << localVertices.size() << " vertices." << endl;



            // Now that we have the list of vertices with maxDistance of vStart,
            // compute connected components, disregarding edges that involve v0.

            // Initialize the disjoint set data structures.
            const ReadId n = ReadId(localVertices.size());
            rank.resize(n);
            parent.resize(n);
            boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
            for(ReadId i=0; i<n; i++) {
                disjointSets.make_set(i);
            }

            // Loop over all edges involving the vertices we found during the BFS,
            // but disregarding vertices involving vStart.
            for(const auto& p: localVertices) {
                const ReadId v0 = p.first;
                if(v0 == vStart) {
                    continue;   // Skip edges involving vStart.
                }
                const ReadId u0 = vertexTable[v0];
                CZI_ASSERT(u0 != notReached);
                const auto edges = readGraphConnectivity[v0];
                for(const uint32_t edgeId: edges) {
                    const ReadGraphEdge& edge = readGraphEdges[edgeId];
                    const AlignmentData& alignment = alignmentData[edge.alignmentId];
                    const ReadId v1 = alignment.getOther(v0);
                    if(v1 == vStart) {
                        continue;   // Skip edges involving vStart.
                    }
                    const ReadId u1 = vertexTable[v1];
                    if(u1 != notReached) {
                        disjointSets.union_set(u0, u1);
                    }
                }
            }


            // Now check the vertices at maximum distance.
            // If they belong to more than one connected component,
            // removing vStart affects the large scale connectivity of the
            // read graph, and therefore we flag vStart as chimeric.
            ReadId component = std::numeric_limits<ReadId>::max();
            for(const auto& p: localVertices) {
                if(p.second != maxDistance) {
                    continue;
                }
                const ReadId v = p.first;
                const ReadId u = vertexTable[v];
                CZI_ASSERT(u != notReached);
                const ReadId uComponent = disjointSets.find_set(u);
                if(component == std::numeric_limits<ReadId>::max()) {
                    component = uComponent;
                } else {
                    if(uComponent != component) {
                        isChimericRead[vStart] = true;
                        out << "Flagged read " << vStart << " as chimeric." << endl;
                        break;
                    }
                }
            }


            // Before processing the next read, we need to reset
            // all entries of the distance vector to not found,
            // then clear the verticesFound vector.
            for(const auto& p: localVertices) {
                const ReadId readId = p.first;
                vertexTable[readId] = notReached;;
            }
            localVertices.clear();
        }
    }

    // Remove our work vector.
    vertexTable.remove();

}



void Assembler::accessChimericReadsFlags()
{
    isChimericRead.accessExistingReadOnly(
        largeDataName("IsChimericRead"));
}



// Compute connected components of the read graph.
// This treats chimeric reads as isolated.
void Assembler::computeReadGraphConnectedComponents()
{
    // Check that we have what we need.
    checkReadGraphIsOpen();
    const size_t readCount = reads.size();
    CZI_ASSERT(readGraphConnectivity.size() == readCount);
    checkAlignmentDataAreOpen();

    // Compute the raw sequence length of each vertex (read).
    cout << timestamp << "Computing the raw read length of each read." << endl;
    vector<size_t> vertexRawSequenceLength(readCount, 0);
    for(ReadId readId=0; readId<readCount; readId++) {
        const MemoryAsContainer<uint8_t> repeatCounts = readRepeatCounts[readId];
        size_t rawSequenceLength = 0;
        for(uint8_t repeatCount: repeatCounts) {
            rawSequenceLength += repeatCount;
        }
        vertexRawSequenceLength[readId] = rawSequenceLength;
    }
    const size_t totalRawSequenceLength = std::accumulate(
        vertexRawSequenceLength.begin(), vertexRawSequenceLength.end(), 0ULL);
    cout << "Found " << readCount;
    cout << " reads for a total " << totalRawSequenceLength << " bases." << endl;



    // Compute connected components of the read graph,
    // threating chimeric reads as isolated.
    vector<ReadId> rank(readCount);
    vector<ReadId> parent(readCount);
    boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    cout << timestamp << "Computing connected components of the read graph." << endl;
    for(ReadId readId=0; readId<readCount; readId++) {
        disjointSets.make_set(readId);
    }
    for(const ReadGraphEdge& edge: readGraphEdges) {
        const AlignmentData& alignment = alignmentData[edge.alignmentId];
        const ReadId readId0 = alignment.readIds[0];
        const ReadId readId1 = alignment.readIds[1];
        if(isChimericRead[readId0]) {
            continue;
        }
        if(isChimericRead[readId1]) {
            continue;
        }
        disjointSets.union_set(readId0, readId1);
    }



    // Gather the vertices of each component.
    std::map<ReadId, vector<ReadId> > componentMap;
    for(ReadId readId=0; readId<readCount; readId++) {
        const ReadId componentId = disjointSets.find_set(readId);
        componentMap[componentId].push_back(readId);
    }
    cout << "The read graph has " << componentMap.size() <<
        " connected components." << endl;



    // Sort the components by decreasing size,
    // where size = total raw sequence length.
    // componentTable contains pairs(size, componentId as key in componentMap).
    vector< pair<size_t, ReadId> > componentTable;
    for(const auto& p: componentMap) {
        const vector<ReadId>& component = p.second;

        size_t componentRawSequenceLength = 0;
        for(const ReadId readId: component) {
            componentRawSequenceLength += vertexRawSequenceLength[readId];
        }

        componentTable.push_back(make_pair(componentRawSequenceLength, p.first));
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<size_t, ReadId>>());



    // Store components in this order of decreasing size.
    vector< vector<ReadId> > components;
    for(const auto& p: componentTable) {
        components.push_back(componentMap[p.second]);
    }
    cout << timestamp << "Done computing connected components of the read graph." << endl;



    // Write information for each component.
    ofstream csv("ReadGraphComponents.csv");
    csv << "Component,RepresentingRead,ReadCount,RawSequenceLength,"
        "AccumulatedReadCount,AccumulatedRawSequenceLength,"
        "AccumulatedReadCountFraction,AccumulatedRawSequenceLengthFraction\n";
    size_t accumulatedReadCount = 0;
    size_t accumulatedRawSequenceLength = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<ReadId>& component = components[componentId];
        accumulatedReadCount += component.size();
        const size_t componentRawSequenceLength = componentTable[componentId].first;
        accumulatedRawSequenceLength += componentRawSequenceLength;
        const double accumulatedReadCountFraction =
            double(accumulatedReadCount)/double(readCount);
        const double accumulatedRawSequenceFraction =
            double(accumulatedRawSequenceLength)/double(totalRawSequenceLength);


        // Write out.
        csv << componentId << ",";
        csv << component.front() << ",";
        csv << component.size() << ",";
        csv << componentRawSequenceLength << ",";
        csv << accumulatedReadCount << ",";
        csv << accumulatedRawSequenceLength << ",";
        csv << accumulatedReadCountFraction << ",";
        csv << accumulatedRawSequenceFraction << "\n";
    }
}
