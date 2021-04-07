
// Shasta.
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "LocalReadGraph.hpp"
#include "orderPairs.hpp"
#include "shastaLapack.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/maximum_adjacency_search.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/property_map/property_map.hpp>



// Standard libraries.
#include "chrono.hpp"
#include "iterator.hpp"
#include <numeric>
#include <queue>
#include <random>
#include <stack>



// For each read, keep only the best maxAlignmentCount alignments.
// Note that the connectivity of the resulting read graph can
// be more than maxAlignmentCount.
void Assembler::createReadGraph(
    uint32_t maxAlignmentCount,
    uint32_t maxTrim)
{
    // Find the number of reads and oriented reads.
    const ReadId orientedReadCount = uint32_t(markers.size());
    SHASTA_ASSERT((orientedReadCount % 2) == 0);
    const ReadId readCount = orientedReadCount / 2;

    // Mark all alignments as not to be kept.
    vector<bool> keepAlignment(alignmentData.size(), false);

    // Vector to keep the alignments for each read,
    // with their number of markers.
    // Contains pairs(marker count, alignment id).
    vector< pair<uint32_t, uint32_t> > readAlignments;

    const bool debug = false;
    if(debug) {
        cout << "createReadGraph begins, maxAlignmentCount " << maxAlignmentCount << endl;
    }


    // Loop over reads.
    for(ReadId readId=0; readId<readCount; readId++) {
        if(debug) {
            cout << "Working on read " << readId << endl;
        }

        // Gather the alignments for this read, each with its number of markers.
        readAlignments.clear();
        for(const uint32_t alignmentId: alignmentTable[OrientedReadId(readId, 0).getValue()]) {
            const AlignmentData& alignment = alignmentData[alignmentId];
            readAlignments.push_back(make_pair(alignment.info.markerCount, alignmentId));
        }
        if(debug) {
            cout << "Found " << readAlignments.size() << " alignments." << endl;
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
        if(debug) {
            cout << "Kept " << readAlignments.size() << " alignments." << endl;
        }

        // Mark the surviving alignments as to be kept.
        for(const auto& p: readAlignments) {
            const uint32_t alignmentId = p.second;
            keepAlignment[alignmentId] = true;
            if(debug) {
                const AlignmentData& alignment = alignmentData[alignmentId];
                cout << "Marked alignment " << alignment.readIds[0] << " " <<
                    alignment.readIds[1] << (alignment.isSameStrand ? " same strand" : " opposite strand") << endl;
            }
        }
    }
    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;

    // Create the read graph using the alignments we selected.
    createReadGraphUsingSelectedAlignments(keepAlignment);
}



// This is called for ReadGraph.creationMethod 0 and 2.
void Assembler::createReadGraphUsingSelectedAlignments(vector<bool>& keepAlignment)
{


    // Now we can create the read graph.
    // Only the alignments we marked as "keep" generate edges in the read graph.
    readGraph.edges.createNew(largeDataName("ReadGraphEdges"), largeDataPageSize);
    for(size_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {

        // Record whether this alignment is used in the read graph.
        const bool keepThisAlignment = keepAlignment[alignmentId];
        AlignmentData& alignment = alignmentData[alignmentId];
        alignment.info.isInReadGraph = uint8_t(keepThisAlignment);

        // If this alignment is not used in the read graph, we are done.
        if(not keepThisAlignment) {
            continue;
        }

        // Create the edge corresponding to this alignment.
        ReadGraphEdge edge;
        edge.alignmentId = alignmentId & 0x3fff'ffff'ffff'ffff;
        edge.crossesStrands = 0;
        edge.hasInconsistentAlignment = 0;
        edge.orientedReadIds[0] = OrientedReadId(alignment.readIds[0], 0);
        edge.orientedReadIds[1] = OrientedReadId(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
        SHASTA_ASSERT(edge.orientedReadIds[0] < edge.orientedReadIds[1]);
        readGraph.edges.push_back(edge);

        // Also create the reverse complemented edge.
        edge.orientedReadIds[0].flipStrand();
        edge.orientedReadIds[1].flipStrand();
        SHASTA_ASSERT(edge.orientedReadIds[0] < edge.orientedReadIds[1]);
        readGraph.edges.push_back(edge);
    }

    // Release unused allocated memory
    readGraph.unreserve();

    // Create read graph connectivity.
    readGraph.connectivity.createNew(largeDataName("ReadGraphConnectivity"), largeDataPageSize);
    readGraph.connectivity.beginPass1(2 * reads->readCount());
    for(const ReadGraphEdge& edge: readGraph.edges) {
        readGraph.connectivity.incrementCount(edge.orientedReadIds[0].getValue());
        readGraph.connectivity.incrementCount(edge.orientedReadIds[1].getValue());
    }
    readGraph.connectivity.beginPass2();
    for(size_t i=0; i<readGraph.edges.size(); i++) {
        const ReadGraphEdge& edge = readGraph.edges[i];
        readGraph.connectivity.store(edge.orientedReadIds[0].getValue(), uint32_t(i));
        readGraph.connectivity.store(edge.orientedReadIds[1].getValue(), uint32_t(i));
    }
    readGraph.connectivity.endPass2();

    // Count the number of isolated reads and their bases.
    uint64_t isolatedReadCount = 0;
    uint64_t isolatedReadBaseCount = 0;
    for(ReadId readId=0; readId<reads->readCount(); readId++) {
        const OrientedReadId orientedReadId(readId, 0);
        const uint64_t neighborCount = readGraph.connectivity.size(orientedReadId.getValue());
        if(neighborCount > 0) {
            continue;
        }
        ++isolatedReadCount;
        isolatedReadBaseCount += reads->getReadRawSequenceLength(readId);
    }
    assemblerInfo->isolatedReadCount = isolatedReadCount;
    assemblerInfo->isolatedReadBaseCount = isolatedReadBaseCount;
}



void Assembler::accessReadGraph()
{
    readGraph.edges.accessExistingReadOnly(largeDataName("ReadGraphEdges"));
    readGraph.connectivity.accessExistingReadOnly(largeDataName("ReadGraphConnectivity"));
}
void Assembler::accessReadGraphReadWrite()
{
    readGraph.edges.accessExistingReadWrite(largeDataName("ReadGraphEdges"));
    readGraph.connectivity.accessExistingReadWrite(largeDataName("ReadGraphConnectivity"));
}
void Assembler::checkReadGraphIsOpen()
{
    if(!readGraph.edges.isOpen) {
        throw runtime_error("Read graph edges are not accessible.");
    }
    if(!readGraph.connectivity.isOpen()) {
        throw runtime_error("Read graph connectivity is not accessible.");
    }

}


// Create a local subgraph of the global read graph,
// starting at a given vertex and extending out to a specified
// distance (number of edges).
bool Assembler::createLocalReadGraph(
        OrientedReadId start,
        uint32_t maxDistance,           // How far to go from starting oriented read.
        bool allowChimericReads,
        bool allowCrossStrandEdges,
        double timeout,                 // Or 0 for no timeout.
        LocalReadGraph& graph)
{
    const vector<OrientedReadId> starts = {start};
    bool success = createLocalReadGraph(
            starts,
            maxDistance,           // How far to go from starting oriented read.
            allowChimericReads,
            allowCrossStrandEdges,
            timeout,                 // Or 0 for no timeout.
            graph
    );

    return success;
}


// Create a local subgraph of the global read graph,
// starting at a given vertex and extending out to a specified
// distance (number of edges).
bool Assembler::createLocalReadGraph(
    const vector<OrientedReadId>& starts,
    uint32_t maxDistance,           // How far to go from starting oriented read.
    bool allowChimericReads,
    bool allowCrossStrandEdges,
    double timeout,                 // Or 0 for no timeout.
    LocalReadGraph& graph)
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
        graph.addVertex(start, uint32_t(markers[start.getValue()].size()),
                        reads->getFlags(start.getReadId()).isChimeric, 0);

        // Add each starting vertex to the BFS queue
        q.push(start);
    }

    // Do the BFS.
    while (!q.empty()) {

        // See if we exceeded the timeout.
        if (timeout > 0. && (seconds(steady_clock::now() - startTime) > timeout)) {
            graph.clear();
            return false;
        }

        // Dequeue a vertex.
        const OrientedReadId orientedReadId0 = q.front();
        q.pop();
        const uint32_t distance0 = graph.getDistance(orientedReadId0);
        const uint32_t distance1 = distance0 + 1;

        // Loop over edges of the global read graph involving this vertex.
        for (const uint64_t i: readGraph.connectivity[orientedReadId0.getValue()]) {
            SHASTA_ASSERT(i < readGraph.edges.size());
            const ReadGraphEdge& globalEdge = readGraph.edges[i];

            if (!allowCrossStrandEdges && globalEdge.crossesStrands) {
                continue;
            }

            // Get the other oriented read involved in this edge of the read graph.
            const OrientedReadId orientedReadId1 = globalEdge.getOther(orientedReadId0);

            // If this read is flagged chimeric and we don't allow chimeric reads, skip.
            if (!allowChimericReads && reads->getFlags(orientedReadId1.getReadId()).isChimeric) {
                continue;
            }

            // Get alignment information.
            const AlignmentData& alignment = alignmentData[globalEdge.alignmentId];
            OrientedReadId alignmentOrientedReadId0(alignment.readIds[0], 0);
            OrientedReadId alignmentOrientedReadId1(alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
            AlignmentInfo alignmentInfo = alignment.info;
            if (alignmentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
                swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
                alignmentInfo.swap();
            }
            if (alignmentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
                alignmentOrientedReadId0.flipStrand();
                alignmentOrientedReadId1.flipStrand();
                alignmentInfo.reverseComplement();
            }
            SHASTA_ASSERT(alignmentOrientedReadId0 == orientedReadId0);
            const uint32_t markerCount = alignmentInfo.markerCount;

            // Update our BFS.
            // Note that we are pushing to the queue vertices at maxDistance,
            // so we can find all of their edges to other vertices at maxDistance.
            if (distance0 < maxDistance) {
                if (!graph.vertexExists(orientedReadId1)) {
                    graph.addVertex(orientedReadId1,
                                    uint32_t(markers[orientedReadId1.getValue()].size()),
                                    reads->getFlags(orientedReadId1.getReadId()).isChimeric, distance1);
                    q.push(orientedReadId1);
                }
                graph.addEdge(
                        orientedReadId0,
                        orientedReadId1,
                        markerCount,
                        i,
                        globalEdge.crossesStrands == 1);
            } else {
                SHASTA_ASSERT(distance0 == maxDistance);
                if (graph.vertexExists(orientedReadId1)) {
                    graph.addEdge(
                            orientedReadId0,
                            orientedReadId1,
                            markerCount,
                            i,
                            globalEdge.crossesStrands == 1);
                }
            }
        }
    }
    return true;
}



// Use the read graph to flag chimeric reads.
// For each oriented read and corresponding vertex v0, we do
// a BFS in the read graph up to the specified maxDistance.
// We then compute connected components of the subgraph
// consisting of the vertices reached by the bfs, minus v0
// and possibly its reverse complement.
// If not all the vertices at maximum distance are
// in the same component, the read corresponding to v0
// is flagged as chimeric.
void Assembler::flagChimericReads(size_t maxDistance, size_t threadCount)
{
    cout << timestamp << "Begin flagging chimeric reads, max distance " << maxDistance << endl;

    // Check that we have what we need.
    checkReadGraphIsOpen();
    const size_t orientedReadCount = readGraph.connectivity.size();
    SHASTA_ASSERT((orientedReadCount % 2) == 0);
    const size_t readCount = orientedReadCount / 2;

    // If maxDistance is zero, just flag all reads as not chimeric.
    if(maxDistance == 0) {
        for(ReadId readId=0; readId<readCount; readId++) {
            reads->setChimericFlag(readId, false);
        }
        return;
    }

    // Store the argument so it is accessible by all threads.
    SHASTA_ASSERT(maxDistance < 255);
    flagChimericReadsData.maxDistance = maxDistance;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Multithreaded loop over all reads.
    cout << timestamp << "Processing " << readCount << " reads." << endl;
    setupLoadBalancing(readCount, 10000);
    runThreads(&Assembler::flagChimericReadsThreadFunction, threadCount);

    cout << timestamp << "Done flagging chimeric reads." << endl;

    size_t chimericReadCount = 0;
    for(ReadId readId=0; readId!=readCount; readId++) {
        if(reads->getFlags(readId).isChimeric) {
            ++chimericReadCount;
        }
    }
    assemblerInfo->chimericReadCount = chimericReadCount;
    cout << timestamp << "Flagged " << chimericReadCount << " reads as chimeric out of ";
    cout << readCount << " total." << endl;
    cout << "Chimera rate is " << double(chimericReadCount) / double(readCount) << endl;
}



void Assembler::flagChimericReadsThreadFunction(size_t threadId)
{
    const size_t maxDistance = flagChimericReadsData.maxDistance;

    // Vector used for BFS searches by this thread.
    // It stores the local vertex id in the current BFS assigned to each vertex,
    // or notReached for vertices not yet reached by the current BFS.
    // Indexed by orientedRead.getValue().
    // This is of size equal to the number of oriented reads, and each thread has its own copy.
    // This is not prohibitive. For example, for a large human size run with
    // 20 million reads and 100 threads, the total space is only 16 GB.
    MemoryMapped::Vector<uint32_t> vertexTable;
    vertexTable.createNew(
        largeDataName("tmp-FlagChimericReads-VertexTable" + to_string(threadId)),
        largeDataPageSize);
    vertexTable.resize(readGraph.connectivity.size());
    const uint32_t notReached = std::numeric_limits<uint32_t>::max();
    fill(vertexTable.begin(), vertexTable.end(), notReached);

    // Vector to contain the vertices we found in the current BFS,
    // each with the distance from the start vertex.
    vector< pair<OrientedReadId, uint32_t> > localVertices;

    // The queue used for the BFS.
    std::queue<OrientedReadId> q;

    // Vectors used to compute connected components after each BFS.
    vector<uint32_t> rank;
    vector<uint32_t> parent;


    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.
        for(ReadId startReadId=ReadId(begin); startReadId!=ReadId(end); startReadId++) {

            // Check that there is no garbage left by the previous BFS.
            SHASTA_ASSERT(localVertices.empty());
            SHASTA_ASSERT(q.empty());

            // Begin by flagging this read as not chimeric.
            reads->setChimericFlag(startReadId, false);



            // Do the BFS for this read and strand 0.
            const OrientedReadId startOrientedReadId(startReadId, 0);
            uint32_t localVertexId = 0;
            q.push(startOrientedReadId);
            localVertices.push_back(make_pair(startOrientedReadId, 0));
            vertexTable[startOrientedReadId.getValue()] = localVertexId++;
            while(!q.empty()) {

                // Dequeue a vertex.
                const OrientedReadId v0 = q.front();
                q.pop();
                const uint32_t distance0 = localVertices[vertexTable[v0.getValue()]].second;
                const uint32_t distance1 = distance0 + 1;
                // out << "Dequeued " << v0 << endl;

                // Loop over edges involving this vertex.
                const auto edgeIds = readGraph.connectivity[v0.getValue()];
                for(const uint32_t edgeId: edgeIds) {
                    const ReadGraphEdge& edge = readGraph.edges[edgeId];
                    if(edge.crossesStrands) {
                        continue;
                    }
                    const OrientedReadId v1 = edge.getOther(v0);
                    // out << "Found " << v1 << endl;

                    // If we already encountered this read in this BFS, do nothing.
                    if(vertexTable[v1.getValue()] != notReached) {
                        // out << "Previously reached." << endl;
                        continue;
                    }

                    // Record this vertex.
                    // out << "Recording " << v1 << endl;
                    localVertices.push_back(make_pair(v1, distance1));
                    vertexTable[v1.getValue()] = localVertexId++;

                    // If at distance less than maxDistance, also enqueue it.
                    if(distance1 < maxDistance) {
                        // out << "Enqueueing " << v1 << endl;
                        q.push(v1);
                    }
                }
            }
            // out << "BFS found " << localVertices.size() << " vertices." << endl;



            // Now that we have the list of vertices with maxDistance of vStart,
            // compute connected components, disregarding edges that involve v0
            // and possibly its reverse complement.

            // Initialize the disjoint set data structures.
            const ReadId n = ReadId(localVertices.size());
            rank.resize(n);
            parent.resize(n);
            boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
            for(ReadId i=0; i<n; i++) {
                disjointSets.make_set(i);
            }

            // Loop over all edges involving the vertices we found during the BFS,
            // but disregarding vertices involving vStart or its reverse complement.
            for(const auto& p: localVertices) {
                const OrientedReadId v0 = p.first;
                if(v0.getReadId() == startOrientedReadId.getReadId()) {
                    continue;   // Skip edges involving vStart or its reverse complement.
                }
                const uint32_t u0 = vertexTable[v0.getValue()];
                SHASTA_ASSERT(u0 != notReached);
                const auto edges = readGraph.connectivity[v0.getValue()];
                for(const uint32_t edgeId: edges) {
                    const ReadGraphEdge& edge = readGraph.edges[edgeId];
                    if(edge.crossesStrands) {
                        continue;
                    }
                    const OrientedReadId v1 = edge.getOther(v0);
                    if(v1.getReadId() == startOrientedReadId.getReadId()) {
                        continue;   // Skip edges involving startOrientedReadId.
                    }
                    const uint32_t u1 = vertexTable[v1.getValue()];
                    if(u1 != notReached) {
                        disjointSets.union_set(u0, u1);
                    }
                }
            }


            // Now check the vertices at maximum distance.
            // If they belong to more than one connected component,
            // removing vStart affects the large scale connectivity of the
            // read graph, and therefore we flag vStart as chimeric.
            uint32_t component = std::numeric_limits<uint32_t>::max();
            for(const auto& p: localVertices) {
                if(p.second != maxDistance) {
                    continue;
                }
                const OrientedReadId v = p.first;
                if(v.getReadId() == startOrientedReadId.getReadId()) {
                    // Skip the reverse complement of the start vertex.
                    continue;
                }
                const uint32_t u = vertexTable[v.getValue()];
                SHASTA_ASSERT(u != notReached);
                const uint32_t uComponent = disjointSets.find_set(u);
                if(component == std::numeric_limits<ReadId>::max()) {
                    component = uComponent;
                } else {
                    if(uComponent != component) {
                        reads->setChimericFlag(startReadId, true);
                        // Also flag all alignments involving this read as not in the read graph.
                        const span<uint32_t> alignmentIds = alignmentTable[OrientedReadId(startReadId, 0).getValue()];
                        for(const uint32_t alignmentId: alignmentIds) {
                            alignmentData[alignmentId].info.isInReadGraph = 0;
                        }
                        break;
                    }
                }
            }


            // Before processing the next read, we need to reset
            // all entries of the distance vector to notReached,
            // then clear the verticesFound vector.
            for(const auto& p: localVertices) {
                const OrientedReadId orientedReadId = p.first;
                vertexTable[orientedReadId.getValue()] = notReached;
            }
            localVertices.clear();
        }
    }

    // Remove our work vector.
    vertexTable.remove();

}



// Compute connected components of the read graph.
// This treats chimeric reads as isolated.
// Components with fewer than minComponentSize are considered
// small and excluded from assembly by setting the
// isInSmallComponent for all the reads they contain.
void Assembler::computeReadGraphConnectedComponents(
    size_t minComponentSize
    )
{
    // Check that we have what we need.
    reads->checkReadFlagsAreOpenForWriting();
    checkReadGraphIsOpen();
    const size_t readCount = reads->readCount();
    const size_t orientedReadCount = 2*readCount;
    SHASTA_ASSERT(readGraph.connectivity.size() == orientedReadCount);
    checkAlignmentDataAreOpen();



    // Compute connected components of the read graph,
    // treating chimeric reads as isolated.
    vector<ReadId> rank(orientedReadCount);
    vector<ReadId> parent(orientedReadCount);
    boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    cout << timestamp << "Computing connected components of the read graph." << endl;
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            disjointSets.make_set(OrientedReadId(readId, strand).getValue());
        }
    }
    for(const ReadGraphEdge& edge: readGraph.edges) {
        if(edge.crossesStrands) {
            continue;
        }
        const OrientedReadId orientedReadId0 = edge.orientedReadIds[0];
        const OrientedReadId orientedReadId1 = edge.orientedReadIds[1];
        const ReadId readId0 = orientedReadId0.getReadId();
        const ReadId readId1 = orientedReadId1.getReadId();
        if(reads->getFlags(readId0).isChimeric) {
            continue;
        }
        if(reads->getFlags(readId1).isChimeric) {
            continue;
        }
        disjointSets.union_set(orientedReadId0.getValue(), orientedReadId1.getValue());
    }



    // Gather the vertices of each component.
    std::map<ReadId, vector<OrientedReadId> > componentMap;
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const ReadId componentId = disjointSets.find_set(orientedReadId.getValue());
            componentMap[componentId].push_back(orientedReadId);
        }
    }
    cout << "The read graph has " << componentMap.size() <<
        " connected components." << endl;



    // Sort the components by decreasing size (number of reads).
    // componentTable contains pairs(size, componentId as key in componentMap).
    vector< pair<size_t, uint32_t> > componentTable;
    for(const auto& p: componentMap) {
        const vector<OrientedReadId>& component = p.second;
        componentTable.push_back(make_pair(component.size(), p.first));
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<size_t, uint32_t>>());



    // Store components in this order of decreasing size.
    vector< vector<OrientedReadId> > components;
    for(const auto& p: componentTable) {
        components.push_back(componentMap[p.second]);
    }
    cout << timestamp << "Done computing connected components of the read graph." << endl;



    // Write information for each component.
    ofstream csv("ReadGraphComponents.csv");
    csv << "Component,RepresentingRead,OrientedReadCount,IsSmall,IsSelfComplementary,"
        "AccumulatedOrientedReadCount,"
        "AccumulatedOrientedReadCountFraction\n";
    size_t accumulatedOrientedReadCount = 0;
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];
        accumulatedOrientedReadCount += component.size();
        const double accumulatedOrientedReadCountFraction =
            double(accumulatedOrientedReadCount)/double(orientedReadCount);

        const bool isSelfComplementary =
            component.size() > 1 &&
            (component[0].getReadId() == component[1].getReadId());


        // Write out.
        csv << componentId << ",";
        csv << component.front() << ",";
        csv << component.size() << ",";
        csv << ((component.size() < minComponentSize) ? "Yes" : "No") << ",";
        csv << (isSelfComplementary ? "Yes" : "No") << ",";
        csv << accumulatedOrientedReadCount << ",";
        csv << accumulatedOrientedReadCountFraction << "\n";
    }



    // Clear the read flags that will be set below.
    // Note that we are not changing the isChimeric flags.
    reads->setIsInSmallComponentFlagForAll(false);
    reads->setStrandFlagForAll(Strand(0));


    // Strand separation. Process the connected components one at a time.
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // If this component is small, set the isInSmallComponent flag for all
        // the reads it contains.
        if(component.size() < minComponentSize) {
            for(const OrientedReadId orientedReadId: component) {
                const ReadId readId = orientedReadId.getReadId();
                reads->setIsInSmallComponentFlag(readId, true);
            }
            continue;
        }

        // Find out if this component is self-complementary.
        const bool isSelfComplementary =
            component.size() > 1 &&
            (component[0].getReadId() == component[1].getReadId());
        if(isSelfComplementary) {
            SHASTA_ASSERT((component.size() % 2) == 0);
        }

        // If this component is not self-complementary,
        // use it for assembly only if its first read is on strand 0.
        // Set the strand of all the reads as
        // the strand present in this component.
        if(!isSelfComplementary) {
            if(component[0].getStrand() == 0) {
                for(const OrientedReadId orientedReadId: component) {
                    const ReadId readId = orientedReadId.getReadId();
                    const Strand strand = orientedReadId.getStrand();
                    reads->setStrandFlag(readId, strand);
                }
            } else {
                // No need to set any strand flags here.
                // They will be set when processing the complementary component.
            }
            continue;
        }

        // If getting here, the component is self-complementary
        // and we need to do strand separation.
        SHASTA_ASSERT(isSelfComplementary);
        cout << "Processing self-complementary component " << componentId <<
            " with " << component.size() << " oriented reads." << endl;

    }



    // Check that any read flagged isChimeric is also flagged isInSmallComponent.
    reads->checkIfAChimericIsAlsoInSmallComponent();
}



// Write a FASTA file containing all reads that appear in
// the local read graph.
void Assembler::writeLocalReadGraphReads(
    ReadId readId,
    Strand strand,
    uint32_t maxDistance,
    bool allowChimericReads,
    bool allowCrossStrandEdges)
{
    // Create the requested local read graph.
    LocalReadGraph localReadGraph;
    SHASTA_ASSERT(createLocalReadGraph(
        OrientedReadId(readId, strand),
        maxDistance,
        allowChimericReads,
        allowCrossStrandEdges,
        0.,
        localReadGraph));

    // Gather the reads.
    std::set<ReadId> readsSet;
    BGL_FORALL_VERTICES(v, localReadGraph, LocalReadGraph) {
        readsSet.insert(localReadGraph[v].orientedReadId.getReadId());
    }



    // Write the fasta file.
    const string fileName = "LocalReadGraph.fasta";
    ofstream fasta(fileName);
    for(const ReadId readId: readsSet) {

        // Write the header line with the read name.
        const auto readName = reads->getReadName(readId);
        fasta << ">" << readId << " ";
        copy(readName.begin(), readName.end(), ostream_iterator<char>(fasta));
        const auto metaData = reads->getReadMetaData(readId);
        if(metaData.size() > 0) {
            fasta << " ";
            copy(metaData.begin(), metaData.end(), ostream_iterator<char>(fasta));
        }
        fasta << "\n";

        // Write the sequence.
        const auto& sequence = reads->getRead(readId);
        const auto& counts = reads->getReadRepeatCounts(readId);
        const size_t n = sequence.baseCount;
        SHASTA_ASSERT(counts.size() == n);
        for(size_t i=0; i<n; i++) {
            const Base base = sequence[i];
            const uint8_t count = counts[i];
            for(size_t k=0; k<count; k++) {
                fasta << base;
            }
        }
        fasta << "\n";
    }
    cout << "Wrote " << readsSet.size() << " reads to " << fileName << endl;


}



void Assembler::flagCrossStrandReadGraphEdges(int maxDistance, size_t threadCount)
{
    const bool debug = false;

    // Initial message.
    cout << timestamp << "Begin flagCrossStrandReadGraphEdges." << endl;

    // Check that we have what we need.
    checkReadGraphIsOpen();
    const size_t readCount = reads->readCount();
    const size_t orientedReadCount = 2*readCount;
    SHASTA_ASSERT(readGraph.connectivity.size() == orientedReadCount);
    checkAlignmentDataAreOpen();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Clear the crossesStrands flag for all read graph edges.
    const size_t edgeCount = readGraph.edges.size();
    for(size_t edgeId=0; edgeId!=edgeCount; edgeId++) {
        readGraph.edges[edgeId].crossesStrands = 0;
    }

    // If maxDistance is 0, don't flag any edges as cross strand edges.
    if(maxDistance == 0) {
        cout << "Skipped flagCrossStrandReadGraphEdges." << endl;
        return;
    }

    // Store the maximum distance so all threads can see it.
    flagCrossStrandReadGraphEdgesData.maxDistance = maxDistance;

    // Find which vertices are close to their reverse complement.
    // "Close" means that there is a path of distance up to maxDistance.
    flagCrossStrandReadGraphEdgesData.isNearStrandJump.clear();
    flagCrossStrandReadGraphEdgesData.isNearStrandJump.resize(orientedReadCount, false);
    const size_t batchSize = 10000;
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::flagCrossStrandReadGraphEdgesThreadFunction, threadCount);
    const auto& isNearStrandJump = flagCrossStrandReadGraphEdgesData.isNearStrandJump;


    size_t nearStrandJumpVertexCount = 0;
    for(ReadId readId=0; readId<readCount; readId++) {
        if(isNearStrandJump[readId]) {
            ++nearStrandJumpVertexCount;
        }
    }
    cout << "Of " << orientedReadCount << " vertices in the read graph, " <<
        nearStrandJumpVertexCount << " are within distance " <<
        maxDistance << " of their reverse complement." << endl;

    // Find connected components of the subgraph consisting of
    // vertices that are close to their reverse complement.
    vector<ReadId> rank(orientedReadCount);
    vector<ReadId> parent(orientedReadCount);
    boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            disjointSets.make_set(OrientedReadId(readId, strand).getValue());
        }
    }
    for(const ReadGraphEdge& edge: readGraph.edges) {
        const OrientedReadId orientedReadId0 = edge.orientedReadIds[0];
        const OrientedReadId orientedReadId1 = edge.orientedReadIds[1];
        const auto v0 = orientedReadId0.getValue();
        const auto v1 = orientedReadId1.getValue();
        if(isNearStrandJump[v0] && isNearStrandJump[v1]) {
            disjointSets.union_set(v0, v1);
        }
    }

    // Gather the vertices in each connected component.
    vector<vector<OrientedReadId> > componentVertices(orientedReadCount);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto v = orientedReadId.getValue();
            if(isNearStrandJump[v]) {
                const ReadId componentId = disjointSets.find_set(v);
                componentVertices[componentId].push_back(orientedReadId);
            }
        }
    }



    // Loop over connected components.
    // Each connected component corresponds to a region of the read graph
    // that has a strand jump.
    // For each such region we process edges in order of decreasing
    // number of markers. We mark an edge as crossing strands if adding it
    // would cause a vertex to become reachable from its reverse complement.
    size_t strandJumpCount = 0;
   for(ReadId componentId=0; componentId!=orientedReadCount; componentId++) {
        const vector<OrientedReadId>& vertices = componentVertices[componentId];
        const size_t vertexCount = vertices.size();
        if(vertexCount <2) {
            continue;
        }
        ++strandJumpCount;
        if(debug) {
            cout << "Found a strand jump region with " << vertexCount <<
                " vertices near read " << vertices.front().getReadId() << "." << endl;
        }

        // Verify that the vertices are a self-complementary set.
        SHASTA_ASSERT((vertexCount %2) == 0);
        for(size_t i=0; i<vertexCount; i+=2) {
            const OrientedReadId orientedReadId0 = vertices[i];
            const OrientedReadId orientedReadId1 = vertices[i+1];
            SHASTA_ASSERT(orientedReadId0.getReadId() == orientedReadId1.getReadId());
            SHASTA_ASSERT(orientedReadId0.getStrand() == 0);
            SHASTA_ASSERT(orientedReadId1.getStrand() == 1);
        }

        // Map the vertices to integers in (0, vertexCount-1).
        std::map<OrientedReadId, uint32_t> vertexMap;
        for(uint32_t i=0; i<vertexCount; i++) {
            vertexMap.insert(make_pair(vertices[i], i));
        }

        // Gather the edges within this region.
        // Store its edge with its alignmentId.
        // This allows us later to match edges into reverse complemented pairs.
        vector< pair<uint32_t, uint64_t> > edgeIds;  // pair(edgeId, alignmentId).
        for(const OrientedReadId orientedReadId0: vertices) {
            const OrientedReadId::Int v0 = orientedReadId0.getValue();
            for(const uint32_t edgeId: readGraph.connectivity[v0]) {
                const ReadGraphEdge& edge = readGraph.edges[edgeId];
                const OrientedReadId orientedReadId1 = edge.getOther(orientedReadId0);
                if(vertexMap.find(orientedReadId1) == vertexMap.end()) {
                    continue;
                }
                if(edge.orientedReadIds[0] == orientedReadId0) { // So we don't add it twice.
                    edgeIds.push_back(make_pair(edgeId, edge.alignmentId));
                }
            }
        }
        if(debug) {
            cout << "This strand jump region contains " << edgeIds.size() << " edges." << endl;
        }

        // Sort them by alignment  id, so pairs of reverse complemented edges come together.
        SHASTA_ASSERT((edgeIds.size() %2) == 0);
        sort(edgeIds.begin(), edgeIds.end(),
            OrderPairsBySecondOnly<uint32_t, uint64_t>());
        for(size_t i=0; i<edgeIds.size(); i+=2){
            SHASTA_ASSERT(edgeIds[i].second == edgeIds[i+1].second);
        }

        // Gather pairs of reverse complemented edges, each with their
        // number of markers.
        vector< pair< array<uint32_t, 2>, uint32_t> > edgePairs;
        for(size_t i=0; i<edgeIds.size(); i+=2){
            const uint64_t alignmentId = edgeIds[i].second;
            SHASTA_ASSERT(alignmentId == edgeIds[i+1].second);
            const uint32_t markerCount = alignmentData[alignmentId].info.markerCount;
            const array<uint32_t, 2> edgePair = {edgeIds[i].first, edgeIds[i+1].first};
            edgePairs.push_back(make_pair(edgePair, markerCount));
        }
        sort(edgePairs.begin(), edgePairs.end(),
            OrderPairsBySecondOnlyGreater<array<uint32_t, 2>, uint32_t>());

        // Initialize a disjoint set data structure for this region.
        vector<ReadId> rank(vertexCount);
        vector<ReadId> parent(vertexCount);
        boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
        for(size_t i=0; i<vertexCount; i++) {
            disjointSets.make_set(i);
        }


        // Process the edge pairs, in order of decreasing number of markers.
        for(const auto& p: edgePairs) {
            const array<uint32_t, 2>& edgeIds = p.first;
            for(const uint32_t edgeId: edgeIds) {
                const ReadGraphEdge& edge = readGraph.edges[edgeId];

                // Get the oriented reads of this edge.
                const OrientedReadId orientedReadId0 = edge.orientedReadIds[0];
                const OrientedReadId orientedReadId1 = edge.orientedReadIds[1];
                const uint32_t i0 = vertexMap[orientedReadId0];
                const uint32_t i1 = vertexMap[orientedReadId1];

                // Get their reverse complemented oriented reads.
                OrientedReadId orientedReadId0rc = orientedReadId0;
                orientedReadId0rc.flipStrand();
                OrientedReadId orientedReadId1rc = orientedReadId1;
                orientedReadId1rc.flipStrand();
                const uint32_t i0rc = vertexMap[orientedReadId0rc];
                const uint32_t i1rc = vertexMap[orientedReadId1rc];

                // Get everybody's component.
                const uint32_t component0 = disjointSets.find_set(i0);
                const uint32_t component1 = disjointSets.find_set(i1);
                const uint32_t component0rc = disjointSets.find_set(i0rc);
                const uint32_t component1rc = disjointSets.find_set(i1rc);

                // Check that we have not already screwed up earlier.
                SHASTA_ASSERT(component0 != component0rc);
                SHASTA_ASSERT(component1 != component1rc);

                // If adding this edge would bring (orientedReadId0, orientedReadId1rc)
                // or (orientedReadId1, orientedReadId0rc)
                // in the same component, mark it as a cross strand edge.
                if(component0==component1rc || component1==component0rc) {
                    ReadGraphEdge& edge = readGraph.edges[edgeId];
                    edge.crossesStrands = 1;
                    // Also mark the corresponding alignment as not in the read graph.
                    const uint64_t alignmentId = edge.alignmentId;
                    alignmentData[alignmentId].info.isInReadGraph = 0;
                } else {
                    disjointSets.union_set(i0, i1);
                    disjointSets.union_set(i0rc, i1rc);
                }
            }
        }
    }
    cout << "Found " << strandJumpCount << " strand jump regions." << endl;

    // Count the number of edges we flagged as cross-strand.
    size_t crossStrandEdgeCount = 0;
    for(size_t edgeId=0; edgeId!=edgeCount; edgeId++) {
        if(readGraph.edges[edgeId].crossesStrands) {
            ++crossStrandEdgeCount;
        }
    }
    cout << "Marked " << crossStrandEdgeCount << " read graph edges out of " <<
        edgeCount <<
        " total as cross-strand." << endl;

    // Done.
    cout << timestamp << "End flagCrossStrandReadGraphEdges." << endl;
}



void Assembler::flagCrossStrandReadGraphEdgesThreadFunction(size_t threadId)
{
    const size_t readCount = reads->readCount();
    const size_t maxDistance = flagCrossStrandReadGraphEdgesData.maxDistance;
    auto& isNearStrandJump = flagCrossStrandReadGraphEdgesData.isNearStrandJump;
    vector<uint32_t> distance(2*readCount, ReadGraph::infiniteDistance);
    vector<OrientedReadId> reachedVertices;
    vector<uint32_t> parentEdges(2*readCount);
    vector<uint32_t> shortestPath;
    uint64_t begin, end;

    while(getNextBatch(begin, end)) {

        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {
            if((readId %100000) == 0) {
                std::lock_guard<std::mutex> lock(mutex);
                cout << timestamp << threadId << " " << readId << "/" << readCount << endl;
            }
            const OrientedReadId orientedReadId0(readId, 0);
            const OrientedReadId orientedReadId1(readId, 1);
            readGraph.computeShortPath(orientedReadId0, orientedReadId1,
                maxDistance, shortestPath,
                distance, reachedVertices, parentEdges);
            if(!shortestPath.empty()) {
                isNearStrandJump[orientedReadId0.getValue()] = true;
                isNearStrandJump[orientedReadId1.getValue()] = true;
            }
        }
    }

}



void Assembler::removeReadGraphBridges(uint64_t maxDistance)
{
    // Check that we have what we need.
    SHASTA_ASSERT(alignmentData.isOpen);
    SHASTA_ASSERT(readGraph.edges.isOpen);
    SHASTA_ASSERT(readGraph.connectivity.isOpen());

    // Flag alignments that are currently in the read graph.
    vector<bool> keepAlignment(alignmentData.size(), false);
    for(const ReadGraphEdge& edge: readGraph.edges) {
        keepAlignment[edge.alignmentId] = true;
    }

    cout << timestamp << "Finding bridges in the read graph." << endl;
    cout << "The read graph uses " <<
        count(keepAlignment.begin(), keepAlignment.end(), true) <<
        " alignments out of " << alignmentData.size() << endl;

    // Unflag alignments corresponding to read graph bridges.
    readGraph.findBridges(keepAlignment, maxDistance);

    // Recreate the read graph using the surviving alignments.
    readGraph.edges.remove();
    readGraph.connectivity.remove();
    createReadGraphUsingSelectedAlignments(keepAlignment);

    cout << timestamp << "After removing bridges, the read graph uses " <<
        count(keepAlignment.begin(), keepAlignment.end(), true) <<
        " alignments out of " << alignmentData.size() << endl;
}



void Assembler::analyzeReadGraph()
{
    // Check that we have what we need.
    SHASTA_ASSERT(readGraph.edges.isOpen);
    SHASTA_ASSERT(readGraph.connectivity.isOpen());

    ofstream csv("AnalyzeReadGraph.csv");
    csv <<
        "OrientedReadId0,OrientedReadId1,"
        "Neighbors0,Neighbors1,"
        "ExclusiveNeighbors0,ExclusiveNeighbors1,"
        "CommonNeighbors,"
        "ExclusiveEdges,"
        "\n";


    // Loop over read graph edges.
    for(uint64_t edgeId=0; edgeId<readGraph.edges.size(); edgeId++) {
        const ReadGraphEdge& edge = readGraph.edges[edgeId];

        // Find neighbors of the two vertices of this edge.
        array<vector<OrientedReadId>, 2> neighbors;
        for(uint64_t i=0; i<2; i++) {
            const OrientedReadId orientedReadId = edge.orientedReadIds[i];
            readGraph.findNeighbors(orientedReadId, neighbors[i]);
        }

        // For each of the two vertex, find neighbors that are not:
        // - The other vertex.
        // - A neighbor of the other vertex.
        array<vector<OrientedReadId>, 2> exclusiveNeighbors;
        for(uint64_t i=0; i<2; i++) {
            const vector<OrientedReadId>& x = neighbors[i];
            const vector<OrientedReadId>& y = neighbors[1-i];
            for(const OrientedReadId orientedReadId: x) {
                if(orientedReadId != edge.orientedReadIds[1-i] and
                    not binary_search(y.begin(), y.end(), orientedReadId)) {
                    exclusiveNeighbors[i].push_back(orientedReadId);
                }
            }
        }

        // Find common neighbors.
        vector<OrientedReadId> commonNeighbors;
        set_intersection(
            neighbors[0].begin(), neighbors[0].end(),
            neighbors[1].begin(), neighbors[1].end(),
            back_inserter(commonNeighbors));


        // Count edges between exclusive neighbors.
        uint64_t exclusiveEdgeCount = 0;
        for(const OrientedReadId orientedReadId0: exclusiveNeighbors[0]) {
            for(const uint32_t edgeId: readGraph.connectivity[orientedReadId0.getValue()]) {
                const ReadGraphEdge& edge = readGraph.edges[edgeId];
                const OrientedReadId orientedReadId1 = edge.getOther(orientedReadId0);
                if(binary_search(exclusiveNeighbors[1].begin(), exclusiveNeighbors[1].end(), orientedReadId1)) {
                    ++exclusiveEdgeCount;
                }
            }
        }



        for(uint64_t i=0; i<2; i++) {
            csv << edge.orientedReadIds[i] << ",";
        }
        for(uint64_t i=0; i<2; i++) {
            csv << neighbors[i].size()-1 << ",";
        }
        for(uint64_t i=0; i<2; i++) {
            csv << exclusiveNeighbors[i].size() << ",";
        }
        csv << commonNeighbors.size() << ",";
        csv << exclusiveEdgeCount << ",";
        csv << "\n";
    }

}


#if 0
// This version runs the clustering once.
void Assembler::readGraphClustering()
{
    SHASTA_ASSERT(readGraph.edges.isOpen);
    SHASTA_ASSERT(readGraph.connectivity.isOpen());

    const uint32_t seed = 231;
    std::mt19937 randomSource(seed);

    vector<ReadId> cluster;
    const bool debug = true;
    readGraph.clustering(randomSource, cluster, debug);
}
#endif


// This version runs the clustering many times.
void Assembler::readGraphClustering()
{
    SHASTA_ASSERT(readGraph.edges.isOpen);
    SHASTA_ASSERT(readGraph.connectivity.isOpen());

    const uint32_t seed = 231;
    std::mt19937 randomSource(seed);

    vector<ReadId> cluster;
    const bool debug = false;

    // Vector to count, for each edge, how many times
    // the two vertices belong to the same cluster.
    vector<uint64_t> isSameClusterEdge(readGraph.edges.size(), 0);

    // Do the clustering many times.
    const uint64_t iterationCount = 1000;
    for(uint64_t iteration=0; iteration<iterationCount; iteration++) {
        readGraph.clustering(randomSource, cluster, debug);

        // Increment isSameClusterEdge counters for each edge.
        for(uint64_t edgeId=0; edgeId<readGraph.edges.size(); edgeId++) {
            const ReadGraphEdge& edge = readGraph.edges[edgeId];
            const OrientedReadId orientedReadId0 = edge.orientedReadIds[0];
            const OrientedReadId orientedReadId1 = edge.orientedReadIds[1];
            const ReadId cluster0 = cluster[orientedReadId0.getValue()];
            const ReadId cluster1 = cluster[orientedReadId1.getValue()];
            if(cluster0 == cluster1) {
                ++isSameClusterEdge[edgeId];
            }
        }
    }


    // Histogram isSameClusterEdge.
    vector<uint64_t> histogram(iterationCount, 0);
    for(uint64_t edgeId=0; edgeId<readGraph.edges.size(); edgeId++) {
        histogram[isSameClusterEdge[edgeId]]++;
    }
    ofstream csv("Histogram.csv");
    for(uint64_t i=0; i<histogram.size(); i++) {
        csv << i << "," << histogram[i] << "\n";
    }



    // Write the read graph in graphviz format, with edge thickness
    // proportional to isSameClusterEdge.
    ofstream graphOut("ReadGraph.dot");
    graphOut << "graph ReadGraph {\n"
        "tooltip=\" \"";
    const ReadId vertexCount = ReadId(readGraph.connectivity.size());
    for(ReadId vertexId=0; vertexId<vertexCount; vertexId++) {
        const OrientedReadId orientedReadId = OrientedReadId(vertexId);
        graphOut << "\"" << orientedReadId << "\"[" <<
            " tooltip=\"" << orientedReadId << "\""
            "];\n";
    }
    for(uint64_t edgeId=0; edgeId<readGraph.edges.size(); edgeId++) {
        const ReadGraphEdge& edge = readGraph.edges[edgeId];
        graphOut << "\"" << edge.orientedReadIds[0] << "\"--\"" <<
            edge.orientedReadIds[1] << "\"";
        /*
        const ReadId vertexId0 = edge.orientedReadIds[0].getValue();
        const ReadId vertexId1 = edge.orientedReadIds[1].getValue();
        const ReadId cluster0 = cluster[vertexId0];
        const ReadId cluster1 = cluster[vertexId1];
        if(cluster0 != cluster1) {
            graphOut << " [color=red]";
        }
        */
        if(isSameClusterEdge[edgeId] < uint64_t(0.1 * double(iterationCount))) {
            graphOut << "[color=red]";
        }
        graphOut << ";\n";
    }
    graphOut << "}\n";
}



// Singular value decomposition analysis of the local read graph.

// Call x the vector of estimated center positions for the oriented reads
// corresponding to the vertices of the LocalReadGraph.

// Each edge U--V corresponds to an alignment that gives an equation
// xU - xV = offset, the average offset between the centers.

// The linear system consisting of these equations is overdetermined
// and can only be solved in a least square sense.
// The least square solution is determined up to a constant,
// because we can add a constant to all the x coordinates
// without affecting the alignments.

// We use a singular value decomposition to compute the least square
// solution, as described for eexample here:
// https://www2.math.uconn.edu/~leykekhman/courses/MATH3795/Lectures/Lecture_9_Linear_least_squares_SVD.pdf

// The number of rows, M, of the constraint matrix A is equal to the
// number of edges, because each edge contributes one equation.
// The number of columns, N, is equal to the number of vertices.

void Assembler::leastSquareAnalysis(
    LocalReadGraph& graph,
    vector<double>& S) const
{
    using vertex_descriptor = LocalReadGraph::vertex_descriptor;
    using edge_descriptor = LocalReadGraph::edge_descriptor;

    const double svThreshold = 1.e-3;   // EXPOSE WHEN CODE STABILIZES **********

    // Count vertices and edges. Use int's because this is what Lapack wants.
    const int N = int(num_vertices(graph));
    SHASTA_ASSERT(N > 0);
    const int M = int(num_edges(graph));
    if(M == 0) {
        throw runtime_error("The local read graph has no edges.");
    }

    // Map the vertices to integers in [0, N).
    vector<vertex_descriptor> vertexTable;
    std::map<vertex_descriptor, int> vertexMap;
    // cout << "Vertices:" << endl;
    int j = 0;
    BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
        // cout << j << " " << graph[v].orientedReadId << endl;
        vertexTable.push_back(v);
        vertexMap.insert(make_pair(v, j++));
    }



    // Map the edges to integers in [0, M] and construct the constraint equations.
    // Each equation is of the form xU - xV = offset.
    vector<edge_descriptor> edgeTable;
    vector< tuple<int, int, double> > equations;
    // cout << "Edges:" << endl;
    BGL_FORALL_EDGES(e, graph, LocalReadGraph) {
        edgeTable.push_back(e);

        // Get the vertices of this edge.
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        // Find the corresponding integer indexes.
        const auto it0 = vertexMap.find(v0);
        const auto it1 = vertexMap.find(v1);
        SHASTA_ASSERT(it0 != vertexMap.end());
        SHASTA_ASSERT(it1 != vertexMap.end());
        const int j0 = it0->second;
        const int j1 = it1->second;

        /// Find the corresponding OrientedReadId's.
        const OrientedReadId orientedReadId0 = graph[v0].orientedReadId;
        const OrientedReadId orientedReadId1 = graph[v1].orientedReadId;

        // Get the alignment data.
        const uint64_t globalEdgeId = graph[e].globalEdgeId;
        const ReadGraphEdge& globalEdge = readGraph.edges[globalEdgeId];
        const uint64_t alignmentId = globalEdge.alignmentId;
        const AlignmentInfo alignmentInfo = alignmentData[alignmentId].orient(orientedReadId0, orientedReadId1);

        // Get the offset.
        const double offset = - alignmentInfo.averageOrdinalOffset;
        graph[e].averageAlignmentOffset = offset;

        // Store this equation.
        /*
        cout << equations.size() << " " << orientedReadId0 << " " << orientedReadId1 << " " <<
            j0 << " " << j1 << " " << offset << endl;
        */
        equations.push_back(make_tuple(j0, j1, offset));

    }
    SHASTA_ASSERT(int(equations.size()) == M);



    // Fill in the constraint matrix A and the right-hand side B.
    vector<double> A(M*N, 0.);  // Constraint matrix.
    vector<double> B(M, 0.);    // Right-hand side.
    for(int i=0; i<M; i++) {
        const auto& equation = equations[i];
        const int j0 = get<0>(equation);
        const int j1 = get<1>(equation);
        const double offset = get<2>(equation);

        // Fill in this equation.
        // Use Fortran matrix storage by columns!
        A[j0*M + i] = -1.;
        A[j1*M + i] = +1.;
        B[i] = offset;
    }



    // Compute the SVD.
    const string JOBU = "A";
    const string JOBVT = "A";
    const int LDA = M;
    S.resize(min(M, N));
    vector<double> U(M*M);
    const int LDU = M;
    vector<double> VT(N*N);
    const int LDVT = N;
    const int LWORK = 10 * max(M, N);
    vector<double> WORK(LWORK);
    int INFO = 0;
    dgesvd_(
        JOBU.data(), JOBVT.data(),
        M, N,
        &A[0], LDA, &S[0], &U[0], LDU, &VT[0], LDVT, &WORK[0], LWORK, INFO);
    if(INFO != 0) {
        throw runtime_error("Error " + to_string(INFO) +
            " computing SVD decomposition of local read graph.");
    }
    /*
    cout << "Singular values: " << endl;
    for(const double v: S) {
        cout << v << endl;
    }
    */

    // We know that there must be a zero singular value because
    // the least square solution is defined up to a constant.
    SHASTA_ASSERT(abs(S[N-1]) < 1.e-14);



    // Now that we have the SVD we can compute the minimum norm solution.
    // See page 11 of
    // https://www2.math.uconn.edu/~leykekhman/courses/MATH3795/Lectures/Lecture_9_Linear_least_squares_SVD.pdf

    // Compute BB = UT * B, the right hand side in the space transformed
    // according to the SVD.
    vector<double> BB(M, 0.);
    const string TRANS = "T";
    dgemv_(
        TRANS.data(), M, M,
        1.,
        &U[0], M,
        &B[0], 1,
        0.,
        &BB[0], 1);

    // Solve for XX = VT * X, the solution vector in the space transformed
    // according to the SVD. In this space, the solution is trivial.
    vector<double> XX(N, 0.);
    for(int i=0; i<N-1; i++) {
        const double sv = S[i];
        SHASTA_ASSERT(sv >= 0.);
        if(sv > svThreshold) {
            XX[i] = BB[i] / sv; // Otherwise it stays at zero.
        }
    }

    // Compute the least square solution vector in the untransformed space,
    // X = V * XX
    vector<double> X(N, 0.);
    dgemv_(
        TRANS.data(), N, N,
        1.,
        &VT[0], N,
        &XX[0], 1,
        0.,
        &X[0], 1);

    // cout << "Least square solution vector:" << endl;
    for(int j=0; j<N; j++) {
        // cout << j << " " << graph[vertexTable[j]].orientedReadId << " " << X[j] << endl;
        graph[vertexTable[j]].leastSquarePosition = X[j];
    }

    /*
    cout << "Residuals at edges: " << endl;
    double maxResidual = 0.;
    for(int i=0; i<M; i++) {
        const auto& equation = equations[i];
        const int j0 = get<0>(equation);
        const int j1 = get<1>(equation);
        const double offset = get<2>(equation);

        const double leastSquareOffset = X[j1] - X[j0];
        const double residual = leastSquareOffset - offset;
        maxResidual = max(maxResidual, abs(residual));

        cout << i << " " <<
            graph[vertexTable[j0]].orientedReadId << " " <<
            graph[vertexTable[j1]].orientedReadId << " " <<
            j0 << " " <<
            j1 << " " <<
            offset << " " <<
            leastSquareOffset << " " <<
            residual << endl;
    }
    cout << "Maximum residual absolute value " << maxResidual << endl;
    */
}



// Triangle analysis of the local read graph.
// Returns a vector of triangles and their alignment residuals,
// sorted by decreasing residual.
void Assembler::triangleAnalysis(
    LocalReadGraph& graph,
    vector< pair<array<LocalReadGraph::edge_descriptor, 3>, int32_t> >& triangles) const
{
    using vertex_descriptor = LocalReadGraph::vertex_descriptor;
    // using edge_descriptor = LocalReadGraph::edge_descriptor;

    // Loop over triangles.
    // To make sure each triangle only gets seen once,
    // only consider it if orientedReadId0<orientedReadId1<orientedReadId2.
    triangles.clear();
    BGL_FORALL_VERTICES(v0, graph, LocalReadGraph) {
        const OrientedReadId orientedReadId0 = graph[v0].orientedReadId;
        BGL_FORALL_OUTEDGES(v0, e01, graph, LocalReadGraph) {
            const vertex_descriptor v1 = target(e01, graph);
            const OrientedReadId orientedReadId1 = graph[v1].orientedReadId;
            if(orientedReadId1 <= orientedReadId0) {
                continue;
            }

            // Get the offset of orientedReadId1 relative to orientedReadId0.
            const uint64_t globalEdgeId01 = graph[e01].globalEdgeId;
            const ReadGraphEdge& globalEdge01 = readGraph.edges[globalEdgeId01];
            const uint64_t alignmentId01 = globalEdge01.alignmentId;
            const AlignmentInfo alignmentInfo01 = alignmentData[alignmentId01].orient(orientedReadId0, orientedReadId1);
            const int32_t offset01 = - alignmentInfo01.averageOrdinalOffset;

            BGL_FORALL_OUTEDGES(v1, e12, graph, LocalReadGraph) {
                const vertex_descriptor v2 = target(e12, graph);
                const OrientedReadId orientedReadId2 = graph[v2].orientedReadId;
                if(orientedReadId2 <= orientedReadId1) {
                    continue;
                }

                // Get the offset of orientedReadId2 relative to orientedReadId1.
                const uint64_t globalEdgeId12 = graph[e12].globalEdgeId;
                const ReadGraphEdge& globalEdge12 = readGraph.edges[globalEdgeId12];
                const uint64_t alignmentId12 = globalEdge12.alignmentId;
                const AlignmentInfo alignmentInfo12 = alignmentData[alignmentId12].orient(orientedReadId1, orientedReadId2);
                const int32_t offset12 = - alignmentInfo12.averageOrdinalOffset;

                // Get the offset of orientedReadId2 relative to orientedReadId0.
                const int32_t offset02 = offset01 + offset12;

                BGL_FORALL_OUTEDGES(v2, e23, graph, LocalReadGraph) {
                    const vertex_descriptor v3 = target(e23, graph);
                    if(v3 != v0) {
                        continue;
                    }

                    // Get the offset of orientedReadId0 relative to orientedReadId2.
                    const uint64_t globalEdgeId23 = graph[e23].globalEdgeId;
                    const ReadGraphEdge& globalEdge23 = readGraph.edges[globalEdgeId23];
                    const uint64_t alignmentId23 = globalEdge23.alignmentId;
                    const AlignmentInfo alignmentInfo23 = alignmentData[alignmentId23].orient(orientedReadId2, orientedReadId0);
                    const int32_t offset23 = - alignmentInfo23.averageOrdinalOffset;

                    // Since v3 is the same as v0, the total offset should be zero.
                    const int32_t offsetError = offset02 + offset23;

                    // Store this triangle.
                    triangles.push_back(make_pair(
                        array<LocalReadGraph::edge_descriptor, 3>({e01, e12, e23}), offsetError));
                }
            }
        }
    }

    // Sort them by decreasing absolute value of offset error.
    using Triangle = pair<array<LocalReadGraph::edge_descriptor, 3>, int32_t>;
    sort(triangles.begin(), triangles.end(),
        [](const Triangle& x, const Triangle& y)
        {
            return abs(x.second) > abs(y.second);
        });

}



void Assembler::flagInconsistentAlignments(
    uint64_t triangleErrorThreshold,
    uint64_t leastSquareErrorThreshold,
    uint64_t leastSquareMaxDistance,
    size_t threadCount)
{
    // Check that we have what we need.
    SHASTA_ASSERT(alignmentData.isOpenWithWriteAccess);
    SHASTA_ASSERT(readGraph.edges.isOpenWithWriteAccess);
    SHASTA_ASSERT(readGraph.connectivity.isOpen());

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store the arguments so they are visible to the threads.
    flagInconsistentAlignmentsData.triangleErrorThreshold = triangleErrorThreshold;
    flagInconsistentAlignmentsData.leastSquareErrorThreshold = leastSquareErrorThreshold;
    flagInconsistentAlignmentsData.leastSquareMaxDistance = leastSquareMaxDistance;

    // Compute the alignment offset for each edge of the read graph,
    // oriented with the lowest OrientedReadId first.
    flagInconsistentAlignmentsData.edgeOffset.createNew(
        largeDataName("tmp-FlagInconsistentAlignmentsDataOffset"), largeDataPageSize);
    flagInconsistentAlignmentsData.edgeOffset.resize(readGraph.edges.size());
    setupLoadBalancing(readGraph.edges.size(), 1000);
    runThreads(&Assembler::flagInconsistentAlignmentsThreadFunction1, threadCount);

    // Loop over triangles in the read graph.
    flagInconsistentAlignmentsData.threadEdgeIds.clear();
    flagInconsistentAlignmentsData.threadEdgeIds.resize(threadCount);
    const uint64_t readCount = readGraph.connectivity.size() / 2;
    setupLoadBalancing(readCount, 100);
    runThreads(&Assembler::flagInconsistentAlignmentsThreadFunction2, threadCount);

    // We no longer need the offsets.
    flagInconsistentAlignmentsData.edgeOffset.remove();

    // Gather the inconsistent edge ids found by all threads.
    vector<uint64_t> edgeIds;
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        vector<uint64_t>& threadEdgeIds = flagInconsistentAlignmentsData.threadEdgeIds[threadId];
        copy(threadEdgeIds.begin(), threadEdgeIds.end(),
            back_inserter(edgeIds));
    }
    deduplicate(edgeIds);

    // Flag these edges as inconsistent and the corresponding alignments
    // as not in the read graph.
    cout << "Flagged " << edgeIds.size() << " read graph edges as inconsistent." << endl;
    for(const uint64_t edgeId: edgeIds) {
        ReadGraphEdge& edge = readGraph.edges[edgeId];
        edge.hasInconsistentAlignment = 1;
        alignmentData[edge.alignmentId].info.isInReadGraph = 0;
        cout << edge.orientedReadIds[0] << " " <<
            edge.orientedReadIds[1] << " " << edge.alignmentId << endl;
    }
}



void Assembler::flagInconsistentAlignmentsThreadFunction1(size_t threadId)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all read graph edges assigned to this batch.
        for(uint64_t edgeId=begin; edgeId!=end; edgeId++) {
            const ReadGraphEdge& edge = readGraph.edges[edgeId];

            // Make sure the lowest OrientedRead is first.
            array<OrientedReadId, 2> orientedReadIds = edge.orientedReadIds;
            if(orientedReadIds[1] < orientedReadIds[0]) {
                swap(orientedReadIds[0], orientedReadIds[1]);
            }
            SHASTA_ASSERT(orientedReadIds[0] < orientedReadIds[1]);

            // Orient the alignment accordingly.
            const AlignmentData& ad = alignmentData[edge.alignmentId];
            const AlignmentInfo alignmentInfo = ad.orient(orientedReadIds[0], orientedReadIds[1]);

            // Store the offset.
            flagInconsistentAlignmentsData.edgeOffset[edgeId] = alignmentInfo.averageOrdinalOffset;
        }
    }

}



// Here we loop over triangles.
// We only consider triangles with oriented read ids 012 where:
// - orientedReadId0 is on strand 0.
// - orientedReadId0<orientedReadId1<orientedReadId2.
// This way each pair of reverse complemented triangles gets looked at exactly once.
// This code is written with a triple loop for each start vertex,
// but could be made faster if necessary.
// In the triple loop, we exclude vertices corresponding to chimeric reads
// and edges marked as cross-strand edges.

void Assembler::flagInconsistentAlignmentsThreadFunction2(size_t threadId)
{
    using vertex_descriptor = LocalReadGraph::vertex_descriptor;
    using edge_descriptor = LocalReadGraph::edge_descriptor;

    const bool debug = false;
    ofstream out;
    if(debug) {
        out.open("flagInconsistentAlignments-" + to_string(threadId) + ".log");
    }

    const uint64_t triangleErrorThreshold = flagInconsistentAlignmentsData.triangleErrorThreshold;
    const double leastSquareErrorThreshold = double(flagInconsistentAlignmentsData.leastSquareErrorThreshold);
    const uint64_t leastSquareMaxDistance = flagInconsistentAlignmentsData.leastSquareMaxDistance;
    vector<uint64_t>& inconsistentEdgeIds = flagInconsistentAlignmentsData.threadEdgeIds[threadId];

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.
        for(ReadId readId0=ReadId(begin); readId0!=ReadId(end); readId0++) {
            if(reads->getFlags(readId0).isChimeric) {
                continue;
            }
            const OrientedReadId orientedReadId0(readId0, 0);

            // Loop over edges of orientedReadId0.
            const span<uint32_t> edgeIds0 = readGraph.connectivity[orientedReadId0.getValue()];
            for(uint32_t edgeId01: edgeIds0){
                const ReadGraphEdge edge01 = readGraph.edges[edgeId01];
                const OrientedReadId orientedReadId1 = edge01.getOther(orientedReadId0);
                if(orientedReadId1 < orientedReadId0) {
                    continue;
                }
                if(reads->getFlags(orientedReadId1.getReadId()).isChimeric) {
                    continue;
                }
                if(edge01.crossesStrands) {
                    continue;
                }
                if(edge01.hasInconsistentAlignment) {
                    continue;
                }
                const int32_t offset01 = flagInconsistentAlignmentsData.edgeOffset[edgeId01];

                // Loop over edges of orientedReadId1.
                const span<uint32_t> edgeIds1 = readGraph.connectivity[orientedReadId1.getValue()];
                for(uint32_t edgeId12: edgeIds1){
                    const ReadGraphEdge edge12 = readGraph.edges[edgeId12];
                    const OrientedReadId orientedReadId2 = edge12.getOther(orientedReadId1);
                    if(orientedReadId2 < orientedReadId1) {
                        continue;
                    }
                    if(reads->getFlags(orientedReadId2.getReadId()).isChimeric) {
                        continue;
                    }
                    if(edge12.crossesStrands) {
                        continue;
                    }
                    if(edge12.hasInconsistentAlignment) {
                        continue;
                    }
                    const int32_t offset12 = flagInconsistentAlignmentsData.edgeOffset[edgeId12];
                    const int32_t offset02 = offset01 + offset12;

                    // Loop over edges of orientedReadId2.
                    const span<uint32_t> edgeIds2 = readGraph.connectivity[orientedReadId2.getValue()];
                    for(uint32_t edgeId20: edgeIds2){
                        const ReadGraphEdge edge20 = readGraph.edges[edgeId20];
                        if(edge20.crossesStrands) {
                            continue;
                        }
                        if(edge20.hasInconsistentAlignment) {
                            continue;
                        }
                        if(edge20.getOther(orientedReadId2) != orientedReadId0) {
                            continue;
                        }

                        // We found a triangle.
                        const int32_t offset20 = -flagInconsistentAlignmentsData.edgeOffset[edgeId20];
                        const int32_t offsetError = offset02 + offset20;

                        // If the error is small, don't do anything.
                        if(abs(offsetError) < triangleErrorThreshold) {
                            continue;
                        }

                        if(debug) {
                            out << "Working on triangle ";
                            out << orientedReadId0 << " ";
                            out << orientedReadId1 << " ";
                            out << orientedReadId2 << " ";
                            out << offset01 << " ";
                            out << offset12 << " ";
                            out << offset20 << " ";
                            out << offsetError << "\n";
                        }

                        // Construct a local read graph around this triangle.
                        LocalReadGraph graph;
                        const vector<OrientedReadId> orientedReadIds =
                            {orientedReadId0, orientedReadId1, orientedReadId2};
                        createLocalReadGraph(orientedReadIds,
                            uint32_t(leastSquareMaxDistance), false, false, 0., graph);

                        // Iterate, removing one edge at a time
                        // until all residuals are small.
                        while(true) {

                            // Perform least square analysis.
                            vector<double> singularValues;
                            leastSquareAnalysis(graph, singularValues);

                            // Find the edge with the worst residual absolute value.
                            double maxResidual = -1.;
                            LocalReadGraph::edge_iterator it, end, itWorst;
                            tie(it, end) = edges(graph);
                            for(; it!=end; ++it) {
                                const edge_descriptor e = *it;
                                const vertex_descriptor v0 = source(e, graph);
                                const vertex_descriptor v1 = target(e, graph);
                                const double x0 = graph[v0].leastSquarePosition;
                                const double x1 = graph[v1].leastSquarePosition;
                                const double residual = abs((x1 - x0) - graph[e].averageAlignmentOffset);
                                if(residual > maxResidual) {
                                    maxResidual = residual;
                                    itWorst = it;
                                }
                            }
                            const edge_descriptor eWorst = *itWorst;
                            const uint64_t globalEdgeId = graph[eWorst].globalEdgeId;
                            if(debug) {
                                 out << "Edge with worst residual " <<
                                    graph[source(eWorst, graph)].orientedReadId << " " <<
                                    graph[target(eWorst, graph)].orientedReadId << " " << maxResidual << endl;
                            }

                            // If the residual is small, end the iteration.
                            if(maxResidual < leastSquareErrorThreshold) {
                                break;
                            }

                            // Remove the edge with the worst residual and its
                            // reverse complement.
                            inconsistentEdgeIds.push_back(globalEdgeId);
                            inconsistentEdgeIds.push_back(readGraph.getReverseComplementEdgeId(globalEdgeId));
                            if(debug) {
                                const ReadGraphEdge& globalEdge = readGraph.edges[globalEdgeId];
                                const AlignmentData& ad = alignmentData[globalEdge.alignmentId];
                                out << "Alignment " << globalEdge.alignmentId << " " <<
                                    ad.readIds[0] << " " << ad.readIds[1] << " " << int(ad.isSameStrand) <<
                                    " flagged as inconsistent." << endl;
                                }
                            remove_edge(eWorst, graph);
                        }
                    }
                }
            }

        }
    }
    deduplicate(inconsistentEdgeIds);
}

