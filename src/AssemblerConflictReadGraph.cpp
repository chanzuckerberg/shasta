#include "Assembler.hpp"
#include "orderPairs.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>



void Assembler::createConflictReadGraph(
    uint64_t threadCount,
    uint32_t maxOffsetSigma,
    uint32_t maxTrim,
    uint32_t maxSkip)
{
    cout << timestamp << "createConflictReadGraph begins." << endl;

    // Check that we have what we need.
    // The code as written only supports the directed read graph.
    SHASTA_ASSERT(directedReadGraph.edges.isOpen);
    SHASTA_ASSERT(directedReadGraph.edgesBySource.isOpen());
    SHASTA_ASSERT(directedReadGraph.edgesByTarget.isOpen());
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();

    // Adjust the number of threads.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store criteria for the induced alignment so all threads can see them.
    createConflictReadGraphData.inducedAlignmentCriteria.maxOffsetSigma = maxOffsetSigma;
    createConflictReadGraphData.inducedAlignmentCriteria.maxTrim = maxTrim;
    createConflictReadGraphData.inducedAlignmentCriteria.maxSkip = maxSkip;

    // Initialize the conflict read graph.
    conflictReadGraph.createNew(largeDataName("ConflictReadGraph"), largeDataPageSize);
    conflictReadGraph.createVertices(readCount());

    // Add edges.
    conflictReadGraph.edges.reserve(10 * readCount());
    setupLoadBalancing(readCount(), 1);
    runThreads(&Assembler::createConflictReadGraphThreadFunction, threadCount);
    conflictReadGraph.edges.unreserve();
    conflictReadGraph.computeConnectivity();


    cout << "The conflict read graph has " <<
        conflictReadGraph.vertices.size() << " vertices and " <<
        conflictReadGraph.edges.size() << " edges" << endl;
    cout << timestamp << "createConflictReadGraph ends." << endl;
}



void Assembler::accessConflictReadGraph()
{
    conflictReadGraph.accessExistingReadWrite(largeDataName("ConflictReadGraph"));

}



void Assembler::createConflictReadGraphThreadFunction(size_t threadId)
{
    const InducedAlignmentCriteria inducedAlignmentCriteria =
        createConflictReadGraphData.inducedAlignmentCriteria;

    // Work areas for addConflictGraphEdges.
    vector<OrientedReadId> conflictCandidates;
    vector<OrientedReadId> conflictingOrientedReads;
    vector<InducedAlignment> inducedAlignments;
    vector<bool> work0;
    vector<bool> work1;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads in this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); ++readId) {
            // cout << "Working on read " << readId << endl;
            addConflictGraphEdges(
                readId,
                inducedAlignmentCriteria,
                conflictCandidates,
                conflictingOrientedReads,
                inducedAlignments,
                work0,
                work1);
        }
    }

}



// This creates edges of the conflict read graph where
// the lower number read is readId0.
// It add the edges to the conflict read graph directly, under mutex protection.
// This should not create significant contention as adding edges to the
// graph is most of the times much faster than computing them.
void Assembler::addConflictGraphEdges(
    ReadId readId0,
    const InducedAlignmentCriteria& inducedAlignmentCriteria,

    // Work areas.
    vector<OrientedReadId>& conflictCandidates,
    vector<OrientedReadId>& conflictingOrientedReads,
    vector<InducedAlignment>& inducedAlignments,
    vector<bool>& work0,
    vector<bool>& work1)
{

    // Put this read on strand 0.
    // When adding edges to the conflict read graph, we will make sure
    // to also add the reverse complemented edge.
    const OrientedReadId orientedReadId0(readId0, 0);



    // Find conflict candidates for orientedReadId0.
    // These are OrientedReadId's that share at least one marker graph vertex
    // with orientedReadId0.
    // To do this, we loop over all markers of orientedReadId0.
    conflictCandidates.clear();
    const MarkerId firstMarkerId = markers.begin(orientedReadId0.getValue()) - markers.begin();
    const uint32_t markerCount = uint32_t(markers.size(orientedReadId0.getValue()));
    for(uint32_t ordinal=0; ordinal<markerCount; ordinal++) {
        const MarkerId markerId0 = firstMarkerId + ordinal;

        // Find the vertex that this marker is on.
        const MarkerGraph::CompressedVertexId compressedVertexId =
            markerGraph.vertexTable[markerId0];

        // If this marker is not on a marker graph vertex, skip.
        if(compressedVertexId == MarkerGraph::invalidCompressedVertexId) {
            continue;
        }

        // Loop over all markers on this vertex.
        const span<MarkerId> vertexMarkers =
            markerGraph.vertices[compressedVertexId];
        for(const MarkerId markerId1: vertexMarkers) {

            // Skip the marker that we started from.
            if(markerId1 == markerId0) {
                continue;
            }

            // Find the oriented read on this marker.
            OrientedReadId orientedReadId1;
            uint32_t ordinal1;
            tie(orientedReadId1, ordinal1) = findMarkerId(markerId1);

            // Only consider markers with readId1>readId0.
            if(orientedReadId1.getReadId() <= readId0) {
                continue;
            }

            // Add this oriented read to our conflict candidates.
           conflictCandidates.push_back(orientedReadId1);
        }
    }
    deduplicate(conflictCandidates);



    // Remove conflict candidates that correspond to an edge of the read graph.
    // For those we already have a good alignment.
    auto itA = conflictCandidates.begin();
    auto itB = itA;
    for(; itA!=conflictCandidates.end(); ++itA) {
        const OrientedReadId orientedReadId1 = *itA;
        const DirectedReadGraph::VertexId v0 = orientedReadId0.getValue();
        const DirectedReadGraph::VertexId v1 = orientedReadId1.getValue();
        const bool forwardExists = directedReadGraph.findEdge(v0, v1)
            != DirectedReadGraph::invalidEdgeId;
        const bool backwardExists = directedReadGraph.findEdge(v1, v0)
            != DirectedReadGraph::invalidEdgeId;

        if(forwardExists or backwardExists) {
            continue;
        } else {
            *itB++ = orientedReadId1;
        }
    }
    conflictCandidates.resize(itB - conflictCandidates.begin());



    // Compute induced alignments between orientedReadId0 and these conflict candidates.
    computeInducedAlignments(
        orientedReadId0,
        conflictCandidates,
        inducedAlignments
    );
    SHASTA_ASSERT(inducedAlignments.size() == conflictCandidates.size());

    // Find which of the induced alignments are bad.
    conflictingOrientedReads.clear();
    const uint32_t markerCount0 = uint32_t(markers.size(orientedReadId0.getValue()));
    for(uint64_t i=0;i<inducedAlignments.size(); i++) {
        const OrientedReadId orientedReadId1 = conflictCandidates[i];
        const uint32_t markerCount1 = uint32_t(markers.size(orientedReadId1.getValue()));
        // cout << "Checking induced alignment of " << orientedReadId0 << " " << orientedReadId1 << endl;

        if(not inducedAlignments[i].evaluate(
            markerCount0,
            markerCount1,
            inducedAlignmentCriteria)) {
            conflictingOrientedReads.push_back(orientedReadId1);
        }
#if 0
        // This also takes into account the presence or absence of marker graph vertices.
        if(not evaluateInducedAlignment(
            orientedReadId0,
            orientedReadId1,
            inducedAlignments[i],
            inducedAlignmentCriteria,
            work0,
            work1)) {
            conflictingOrientedReads.push_back(orientedReadId1);
        }
#endif
    }
    // cout << "Counts: " << conflictCandidates.size () << " " << conflictingOrientedReads.size() << endl;


    // Add edges to the conflict graph.
    {
        // Find the vertices corresponding to the first read.
        using VertexId = ConflictReadGraph::VertexId;
        const VertexId vertexId0 =
            ConflictReadGraph::getVertexId(orientedReadId0);
        OrientedReadId orientedReadId0ReverseComplement = orientedReadId0;
        orientedReadId0ReverseComplement.flipStrand();
        const VertexId vertexId0ReverseComplement =
            ConflictReadGraph::getVertexId(orientedReadId0ReverseComplement);

        std::lock_guard<std::mutex> lock(mutex);
        for(const OrientedReadId orientedReadId1: conflictingOrientedReads) {

            // Find the vertices corresponding to the second read.
            const VertexId vertexId1 =
                ConflictReadGraph::getVertexId(orientedReadId1);
            OrientedReadId orientedReadId1ReverseComplement = orientedReadId1;
            orientedReadId1ReverseComplement.flipStrand();
            const VertexId vertexId1ReverseComplement =
                ConflictReadGraph::getVertexId(orientedReadId1ReverseComplement);

            // Add the edges.
            conflictReadGraph.addEdge(
                vertexId0,
                vertexId1,
                ConflictReadGraphEdge());
            conflictReadGraph.addEdge(
                vertexId0ReverseComplement,
                vertexId1ReverseComplement,
                ConflictReadGraphEdge());
        }

    }
}


#if 0
// Use a greedy coloring algorithm to color each non-trivial
// connected component of the conflict read graph.
// The coloring algorithm is in ConflictReadGraph::colorConnectedComponent
void Assembler::colorConflictReadGraph()
{
    SHASTA_ASSERT(conflictReadGraph.isOpen());
    using VertexId = ConflictReadGraph::VertexId;
    using EdgeId = ConflictReadGraph::EdgeId;
    const auto invalid = ConflictReadGraphVertex::invalid;

    // Initialize all vertex components and colors to invalid.
    const VertexId n = conflictReadGraph.vertices.size();
    for(VertexId vertexId=0; vertexId<n; vertexId++) {
        auto& vertex = conflictReadGraph.getVertex(vertexId);
        vertex.componentId = invalid;
        vertex.color = invalid;
    }


    // Compute connected components.
    vector<VertexId> rank(n);
    vector<VertexId> parent(n);
    boost::disjoint_sets<VertexId*, VertexId*> disjointSets(&rank[0], &parent[0]);
    for(VertexId vertexId=0; vertexId<n; vertexId++) {
        disjointSets.make_set(vertexId);
    }
    for(EdgeId edgeId=0; edgeId<conflictReadGraph.edges.size(); edgeId++) {
        disjointSets.union_set(
            conflictReadGraph.v0(edgeId),
            conflictReadGraph.v1(edgeId));
    }



    // Gather the vertices of each connected component.
    std::map<VertexId, vector<VertexId> > componentMap;
    for(VertexId vertexId=0; vertexId<n; vertexId++) {
        const VertexId componentId = disjointSets.find_set(vertexId);
        componentMap[componentId].push_back(vertexId);
    }




    // Sort the components by decreasing size.
    // Discard trivial connected components containing only one isolated vertex.
    // componentTable contains pairs(size, componentId as key in componentMap).
    vector< pair<VertexId, VertexId> > componentTable;
    for(const auto& p: componentMap) {
        const VertexId componentId = p.first;
        const vector<VertexId>& component = p.second;
        const VertexId componentSize = component.size();
        if(componentSize > 1) {
            componentTable.push_back(make_pair(componentSize, componentId));
        }
    }
    sort(componentTable.begin(), componentTable.end(), std::greater<pair<VertexId, VertexId>>());



    // Store components in this order of decreasing size.
    vector< vector<VertexId> > components;
    VertexId nonTrivialComponentVertexCount = 0;
    for(const auto& p: componentTable) {
        const vector<VertexId>& component = componentMap[p.second];
        components.push_back(component);
        nonTrivialComponentVertexCount += component.size();
    }
    cout << "Found " << components.size() <<
        " non-trivial connected components with a total " <<
        nonTrivialComponentVertexCount <<  " vertices out " << n <<
        " vertices in the conflict read graph." << endl;



    // Color each connected component separately.
    for(VertexId componentId=0; componentId<components.size(); componentId++) {
        const vector<VertexId>& component = components[componentId];

        // Store the component id in each vertex.
        for(const VertexId vertexId: component) {
            conflictReadGraph.getVertex(vertexId).componentId = componentId;
        }

        cout << "Coloring component " << componentId << " with " <<
            component.size() << " vertices." << endl;
        conflictReadGraph.colorConnectedComponent(component);
    }



    // Write graphviz file containing each of the non-trivial connected component.
    for(VertexId componentId=0; componentId<components.size(); componentId++) {
        const vector<VertexId>& component = components[componentId];
        ofstream graphOut("ConflictReadGaph-" + to_string(componentId) + ".dot");
        graphOut << "graph component" << componentId << " {\n";

        // Write the vertices.
        for(const VertexId vertexId: component) {
            const ConflictReadGraphVertex& vertex = conflictReadGraph.getVertex(vertexId);
            graphOut << "\"" << conflictReadGraph.getOrientedReadId(vertexId) << "\" [";

            // Tooltip.
            graphOut << " tooltip=\"" << conflictReadGraph.getOrientedReadId(vertexId);
            if(vertex.componentId != invalid) {
                SHASTA_ASSERT(vertex.color != invalid);
                graphOut << " component " << vertex.componentId << " color " <<
                    vertex.color;
            }
            graphOut << "\"";

            // Colors.
            if(vertex.componentId != invalid) {
                SHASTA_ASSERT(vertex.color != invalid);
                graphOut << " color=\"/set18/" << (vertex.color % 8) + 1 << "\"";
                graphOut << " fillcolor=\"/set18/" << (vertex.color % 8) + 1 << "\"";

            }
            graphOut << "];\n";
        }



        // Write the edges.
        for(const VertexId vertexId: component) {
            const span<EdgeId> incidentEdges =
                conflictReadGraph.incidentEdges(vertexId);
            for(const EdgeId edgeId: incidentEdges) {
                const VertexId v0 = conflictReadGraph.v0(edgeId);
                if(v0 == vertexId) {
                    const VertexId v1 = conflictReadGraph.v1(edgeId);
                    graphOut << "\"" << conflictReadGraph.getOrientedReadId(v0) << "\"--\"" <<
                        conflictReadGraph.getOrientedReadId(v1) << "\";\n";

                    // Also check that the other vertex has the same
                    // component and a different color.
                    const ConflictReadGraphVertex& vertex0 = conflictReadGraph.getVertex(v0);
                    const ConflictReadGraphVertex& vertex1 = conflictReadGraph.getVertex(v1);
                    SHASTA_ASSERT(vertex0.componentId == vertex1.componentId);
                    SHASTA_ASSERT(vertex0.color != vertex1.color);
                }

            }
        }

        graphOut << "}\n";
    }




}
#endif



// This colors the ConflictReadGraph by walking the
// DirectedReadGraph.
void Assembler::colorConflictReadGraph()
{
    const bool debug = false;

    // Check that we have what we need.
    SHASTA_ASSERT(directedReadGraph.isOpen());
    SHASTA_ASSERT(conflictReadGraph.isOpen());

    // Types for vertices and edges of the two graphs we will use.
    using VertexId = DirectedReadGraph::VertexId;
    using EdgeId   = DirectedReadGraph::EdgeId;

    // We are assuming the two graphs use the same VertexId and EdgeId.
    static_assert(std::is_same<VertexId, ConflictReadGraph::VertexId>::value,
        "Unexpected VertexId discrepancy.");
    static_assert(std::is_same<VertexId, ConflictReadGraph::VertexId>::value,
        "Unexpected VertexId discrepancy.");

    // Initialize all vertex components and colors to invalid.
    const auto invalid = ConflictReadGraphVertex::invalid;
    const VertexId n = conflictReadGraph.vertices.size();
    for(VertexId vertexId=0; vertexId<n; vertexId++) {
        auto& vertex = conflictReadGraph.getVertex(vertexId);
        vertex.componentId = invalid;
        vertex.color = invalid;
    }



    // Find clusters of vertices in the DirectedReadGraph that
    // have no conflicts. These will be used as starting points
    // for each iteration in the coloring process.
    vector<VertexId> rank(n);
    vector<VertexId> parent(n);
    boost::disjoint_sets<VertexId*, VertexId*> disjointSets(&rank[0], &parent[0]);
    for(VertexId v=0; v<n; v++) {
        disjointSets.make_set(v);
    }
    for(EdgeId edgeId=0; edgeId<directedReadGraph.edges.size(); edgeId++) {
        const VertexId v0 = directedReadGraph.source(edgeId);
        if(conflictReadGraph.degree(v0) > 0) {
            continue;   // v0 has conflict. Skip.
        }
        const VertexId v1 = directedReadGraph.target(edgeId);
        if(conflictReadGraph.degree(v1) > 0) {
            continue;   // v1 has conflict. Skip.
        }
        disjointSets.union_set(v0, v1);
    }



    // Gather the vertices of each cluster.
    std::map<VertexId, vector<VertexId> > clusterMap;
    for(VertexId vertexId=0; vertexId<n; vertexId++) {
        const VertexId clusterId = disjointSets.find_set(vertexId);
        clusterMap[clusterId].push_back(vertexId);
    }




    // Sort the clusters by decreasing size.
    // clusterTable contains pairs(size, clusterId as key in clusterMap).
    vector< pair<VertexId, VertexId> > clusterTable;
    for(const auto& p: clusterMap) {
        const VertexId clusterId = p.first;
        const vector<VertexId>& cluster = p.second;
        const VertexId clusterSize = cluster.size();
        SHASTA_ASSERT(clusterSize > 0);
        if(clusterSize==1) {
            // The cluster consists of a single vertex.
            const VertexId v = cluster.front();
            if(conflictReadGraph.degree(v) > 0) {
                // That single vertex has conflicts, so don't consider it.
            } else {
                // That single vertex has no conflicts, so store it.
                clusterTable.push_back(make_pair(clusterSize, clusterId));
            }
        } else {
            // The cluster consists of two or more vertices.
            clusterTable.push_back(make_pair(clusterSize, clusterId));
        }
    }
    sort(clusterTable.begin(), clusterTable.end(), std::greater<pair<VertexId, VertexId>>());



    // Store clusters in this order of decreasing size.
    vector< vector<VertexId> > clusters;
    VertexId clustersVertexCount = 0;
    for(const auto& p: clusterTable) {
        const vector<VertexId>& cluster = clusterMap[p.second];
        clusters.push_back(cluster);
        clustersVertexCount += cluster.size();
    }
    if(debug) {
        cout << "Found " << clusters.size() <<
            " clusters of vertices without conflict for a total " <<
            clustersVertexCount <<  " vertices out " << n <<
            " vertices in the read graph." << endl;
    }

    // Data structures to keep track, at each iteration, of the vertices
    // that conflict with vertices we already encountered.
    vector<VertexId> forbiddenVertices;
    vector<bool> isForbidden(directedReadGraph.vertices.size(), false);

    // Other data structures used below.
    vector<VertexId> adjacentVertices;
    vector< pair<VertexId, uint64_t > > adjacentVerticesSortedByConflictCount;



    // Iterate over clusters in this order of decreasing size.
    for(uint64_t iteration=0; iteration<clusters.size(); iteration++) {
        const vector<VertexId>& cluster = clusters[iteration];
        if(debug) {
            cout << "Coloring iteration " << iteration <<
                " starting from a cluster with " << cluster.size() << " vertices." << endl;
            for(const VertexId vertexId: cluster) {
                cout << ConflictReadGraph::getOrientedReadId(vertexId) << " ";
            }
            cout << endl;
        }


        // We use a process similar to a BFS starting at the vertices
        // of this cluster, with a few differences:
        // - When encountering a vertex that conflicts with another vertex
        //   we already encountered at this iteration, we skip it.
        // - When encountering a vertex that was already colored at a previous
        //   iteration (not just one colored at the current iteration), we skip it.
        // - When we enqueue neighbors of a vertex, we enqueue them in order
        //   of increasing number of conflicting vertices.


        // Initialize the BFS.
        std::queue<VertexId> q;
        for(const VertexId vertexId: cluster) {
            q.push(vertexId);
        }


        // BFS loop.
        while(not q.empty()) {

            // Dequeue a vertex.
            const VertexId v0 = q.front();
            const OrientedReadId orientedReadId0 = ConflictReadGraph::getOrientedReadId(v0);
            if(debug) {
                cout << "Queue size " << q.size() << ", dequeued " << orientedReadId0 << endl;
            }
            q.pop();

            // If v0 is now forbidden, skip it.
            if(isForbidden[v0]) {
                // v0 was not forbidden when we enqueud it, but it is forbidden now.
                if(debug) {
                    cout << orientedReadId0 << " skipped because it is now forbidden." << endl;
                }
                continue;
            }

            // Give it a color equal to this iteration.
            conflictReadGraph.getVertex(v0).color = iteration;
            if(debug) {
                cout<< orientedReadId0 << " being colored " << iteration << endl;
            }

            // Gather adjacent vertices.
            directedReadGraph.findKeptAdjacent(v0, adjacentVertices);

            // Sort adjacent vertices by increasing number of conflicts.
            // Skip the ones that are conflicting with vertices we
            // already encountered at this iteration or that have
            // already been colored at a previous iteration.
            adjacentVerticesSortedByConflictCount.clear();
            for(const VertexId v1: adjacentVertices) {
                if(isForbidden[v1]) {
                    if(debug) {
                        cout << ConflictReadGraph::getOrientedReadId(v1) << " forbidden" << endl;
                    }
                    continue;
                }
                if(conflictReadGraph.vertices[v1].color != invalid) {
                    if(debug) {
                        cout << ConflictReadGraph::getOrientedReadId(v1) << " already colored" << endl;
                    }
                    continue;
                }
                if(conflictReadGraph.vertices[v1].componentId != invalid) {
                    if(debug) {
                        cout << ConflictReadGraph::getOrientedReadId(v1) << " already enqueued" << endl;
                    }
                    continue;
                }
                adjacentVerticesSortedByConflictCount.push_back(
                    make_pair(v1, conflictReadGraph.degree(v1)));
            }
            sort(
                adjacentVerticesSortedByConflictCount.begin(),
                adjacentVerticesSortedByConflictCount.end(),
                OrderPairsBySecondOnly<VertexId, uint64_t>());

            // Loop over adjacent vertices, in this order of increasing number of conflicts.
            for(const auto& p: adjacentVerticesSortedByConflictCount) {
                const VertexId v1 = p.first;

                // We know this vertex is not forbidden and was not already
                // color at this or the previous iteration, so we can enqueue it now.
                // Set the component id to indicate that it was enqueued.
                conflictReadGraph.vertices[v1].componentId = iteration;
                q.push(v1);
                if(debug) {
                    cout << "Enqueued " << conflictReadGraph.getOrientedReadId(v1) << endl;
                }

                // Mark as forbidden the vertices that conflict with v1.
                for(const EdgeId e12: conflictReadGraph.edgesByVertex[v1]) {
                    const VertexId v2 = conflictReadGraph.otherVertex(e12, v1);
                    if(not isForbidden[v2]) {
                        isForbidden[v2] = true;
                        forbiddenVertices.push_back(v2);
                        if(debug) {
                            cout << "Marked as forbidden " <<
                                conflictReadGraph.getOrientedReadId(v2) << endl;
                        }
                    }
                }
            }
        }



        // Clean up data structures to prepare them for the next iteration.
        for(const VertexId v: forbiddenVertices) {
            isForbidden[v] = false;
        }
        forbiddenVertices.clear();

    }

}



void Assembler::markDirectedReadGraphConflictEdges()
{
    // Check that we have what we need.
    SHASTA_ASSERT(directedReadGraph.isOpen());
    SHASTA_ASSERT(conflictReadGraph.isOpen());

    // Loop over all edges of the directed read graph.
    const auto invalid = ConflictReadGraphVertex::invalid;
    uint64_t invalidEdgeCount = 0;
    for(DirectedReadGraph::EdgeId edgeId=0; edgeId<directedReadGraph.edges.size(); edgeId++) {
        DirectedReadGraphEdge& edge = directedReadGraph.getEdge(edgeId);

        // Get the vertices of the DirectedReadGraph..
        const DirectedReadGraph::VertexId v0 = directedReadGraph.source(edgeId);
        const DirectedReadGraph::VertexId v1 = directedReadGraph.target(edgeId);

        // Get the corresponding OrientedReadId's.
        const OrientedReadId orientedReadId0 = OrientedReadId(OrientedReadId::Int(v0));
        const OrientedReadId orientedReadId1 = OrientedReadId(OrientedReadId::Int(v1));

        // Get the corresponding vertices of the ConclictReadGraph.
        const ConflictReadGraph::VertexId u0 = ConflictReadGraph::getVertexId(orientedReadId0);
        const ConflictReadGraph::VertexId u1 = ConflictReadGraph::getVertexId(orientedReadId1);
        const ConflictReadGraphVertex& cVertex0 = conflictReadGraph.getVertex(u0);
        const ConflictReadGraphVertex& cVertex1 = conflictReadGraph.getVertex(u1);

        // With current numbering, the vertex ids should be the same.
        SHASTA_ASSERT(u0 == v0);
        SHASTA_ASSERT(u1 == v1);



        // Figure out if this a conflict edge.
        if(cVertex0.componentId == invalid or cVertex1.componentId==invalid) {

            // One or both vertices are isolated in the conflict graph.
            // Therefore there is no conflict.
            edge.isConflict = 0;

        } else {

            // Both vertices belong to non-trivial connected components of the
            // conflict read graph. Check that the color was set.
            SHASTA_ASSERT(cVertex0.color != invalid);
            SHASTA_ASSERT(cVertex1.color != invalid);

            if(cVertex0.componentId == cVertex1.componentId) {

                // The two vertices belong to the same connected components
                // of the conflict read graph. Therefore there is conflict
                // if their colors are different.
                edge.isConflict = (cVertex0.color != cVertex1.color);

            } else {

                // The two vertices belong to different connected components
                // of the conflict read graph. Therefore there is no conflict.
                edge.isConflict = 0;
            }


        }
        if(edge.isConflict) {
            ++invalidEdgeCount;
        }
    }

    cout << "Marked as conflict edges " << invalidEdgeCount << " edges out of " <<
        directedReadGraph.edges.size() <<
        " edges in the directed read graph." << endl;
}
