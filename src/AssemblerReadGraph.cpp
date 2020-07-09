
// Shasta.
#include "Assembler.hpp"
#include "LocalReadGraph.hpp"
#include "orderPairs.hpp"
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
        if(!keepAlignment[alignmentId]) {
            continue;
        }
        const AlignmentData& alignment = alignmentData[alignmentId];

        // Create the edge corresponding to this alignment.
        ReadGraphEdge edge;
        edge.alignmentId = alignmentId & 0x7fff'ffff'ffff'ffff;
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
    readGraph.connectivity.beginPass1(2 * readCount());
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
    for(ReadId readId=0; readId<readCount(); readId++) {
        const OrientedReadId orientedReadId(readId, 0);
        const uint64_t neighborCount = readGraph.connectivity.size(orientedReadId.getValue());
        if(neighborCount > 0) {
            continue;
        }
        ++isolatedReadCount;
        isolatedReadBaseCount += getReadRawSequenceLength(readId);
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
        size_t maxTrim,                 // Used to define containment.
        double timeout,                 // Or 0 for no timeout.
        LocalReadGraph& graph)
{
    vector<OrientedReadId> starts = {start};
    bool success = createLocalReadGraph(
            starts,
            maxDistance,           // How far to go from starting oriented read.
            allowChimericReads,
            allowCrossStrandEdges,
            maxTrim,                 // Used to define containment.
            timeout,                 // Or 0 for no timeout.
            graph
    );

    return success;
}


// Create a local subgraph of the global read graph,
// starting at a given vertex and extending out to a specified
// distance (number of edges).
bool Assembler::createLocalReadGraph(
    vector<OrientedReadId>& starts,
    uint32_t maxDistance,           // How far to go from starting oriented read.
    bool allowChimericReads,
    bool allowCrossStrandEdges,
    size_t maxTrim,                 // Used to define containment.
    double timeout,                 // Or 0 for no timeout.
    LocalReadGraph& graph)
{
    const auto startTime = steady_clock::now();

    // Initialize a BFS starting at the start vertex.
    std::queue<OrientedReadId> q;

    for (auto& start: starts) {
        // If the starting read is chimeric and we don't allow chimeric reads, do nothing.
        if (!allowChimericReads && readFlags[start.getReadId()].isChimeric) {
            continue;
        }

        // Add the starting vertex.
        graph.addVertex(start, uint32_t(markers[start.getValue()].size()),
                        readFlags[start.getReadId()].isChimeric, 0);

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
            if (!allowChimericReads && readFlags[orientedReadId1.getReadId()].isChimeric) {
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
            const AlignmentType alignmentType = alignmentInfo.classify(uint32_t(maxTrim));
            const uint32_t markerCount = alignmentInfo.markerCount;

            // Update our BFS.
            // Note that we are pushing to the queue vertices at maxDistance,
            // so we can find all of their edges to other vertices at maxDistance.
            if (distance0 < maxDistance) {
                if (!graph.vertexExists(orientedReadId1)) {
                    graph.addVertex(orientedReadId1,
                                    uint32_t(markers[orientedReadId1.getValue()].size()),
                                    readFlags[orientedReadId1.getReadId()].isChimeric, distance1);
                    q.push(orientedReadId1);
                }
                graph.addEdge(
                        orientedReadId0,
                        orientedReadId1,
                        markerCount,
                        alignmentType,
                        globalEdge.crossesStrands == 1);
            } else {
                SHASTA_ASSERT(distance0 == maxDistance);
                if (graph.vertexExists(orientedReadId1)) {
                    graph.addEdge(
                            orientedReadId0,
                            orientedReadId1,
                            markerCount,
                            alignmentType,
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
            readFlags[readId].isChimeric = 0;
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
        if(readFlags[readId].isChimeric) {
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
            readFlags[startReadId].isChimeric = 0;



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
                        readFlags[startReadId].isChimeric = 1;
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
    SHASTA_ASSERT(readFlags.isOpenWithWriteAccess);
    checkReadGraphIsOpen();
    const size_t readCount = reads.size();
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
        if(readFlags[readId0].isChimeric) {
            continue;
        }
        if(readFlags[readId1].isChimeric) {
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
    for(ReadFlags& f: readFlags) {
        f.isInSmallComponent = 0;
        f.strand = 0;
    }



    // Strand separation. Process the connected components one at a time.
    for(ReadId componentId=0; componentId<components.size(); componentId++) {
        const vector<OrientedReadId>& component = components[componentId];

        // If this component is small, set the isInSmallComponent flag for all
        // the reads it contains.
        if(component.size() < minComponentSize) {
            for(const OrientedReadId orientedReadId: component) {
                const ReadId readId = orientedReadId.getReadId();
                readFlags[readId].isInSmallComponent = 1;
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
                    readFlags[readId].strand = strand & 1;
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
    for(const ReadFlags& flags: readFlags) {
        if(flags.isChimeric) {
            SHASTA_ASSERT(flags.isInSmallComponent);
        }
    }
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
        std::numeric_limits<size_t>::max(),
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
        const auto readName = readNames[readId];
        fasta << ">" << readId << " ";
        copy(readName.begin(), readName.end(), ostream_iterator<char>(fasta));
        const auto metaData = readMetaData[readId];
        if(metaData.size() > 0) {
            fasta << " ";
            copy(metaData.begin(), metaData.end(), ostream_iterator<char>(fasta));
        }
        fasta << "\n";

        // Write the sequence.
        const auto& sequence = reads[readId];
        const auto& counts = readRepeatCounts[readId];
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
    const size_t readCount = reads.size();
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
                    readGraph.edges[edgeId].crossesStrands = 1;
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
    const size_t readCount = reads.size();
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


