#include "Assembler.hpp"
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
                inducedAlignments);
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
    vector<InducedAlignment>& inducedAlignments)
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
        const MemoryAsContainer<MarkerId> vertexMarkers =
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



// Use a greedy coloring algorithm to color each non-trivial
// connected component of the conflict read graph.
void Assembler::colorConflictReadGraph()
{
    SHASTA_ASSERT(conflictReadGraph.isOpen());
    using VertexId = ConflictReadGraph::VertexId;
    using EdgeId = ConflictReadGraph::EdgeId;

    // Initialize all vertex components and colors to invalid.
    const VertexId n = conflictReadGraph.vertices.size();
    for(VertexId vertexId=0; vertexId<n; vertexId++) {
        auto& vertex = conflictReadGraph.getVertex(vertexId);
        vertex.componentId = ConflictReadGraphVertex::invalid;
        vertex.color = ConflictReadGraphVertex::invalid;
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



    // Write graphviz file containing each of the non-trivial connected component.
    for(VertexId componentId=0; componentId<components.size(); componentId++) {
        const vector<VertexId>& component = components[componentId];
        ofstream graphOut("ConflictReadGaph-" + to_string(componentId) + ".dot");
        graphOut << "graph component" << componentId << " {\n";

        // Write the vertices.
        for(const VertexId vertexId: component) {
            graphOut << "\"" << conflictReadGraph.getOrientedReadId(vertexId) << "\";\n";
        }

        // Write the edges.
        for(const VertexId vertexId: component) {
            const MemoryAsContainer<EdgeId> incidentEdges =
                conflictReadGraph.incidentEdges(vertexId);
            for(const EdgeId edgeId: incidentEdges) {
                const VertexId v0 = conflictReadGraph.v0(edgeId);
                if(v0 == vertexId) {
                    const VertexId v1 = conflictReadGraph.v1(edgeId);
                    graphOut << "\"" << conflictReadGraph.getOrientedReadId(v0) << "\"--\"" <<
                        conflictReadGraph.getOrientedReadId(v1) << "\";\n";
                }

            }
        }

        graphOut << "}\n";
    }




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

}
