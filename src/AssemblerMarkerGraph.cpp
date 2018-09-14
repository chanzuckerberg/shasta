// shasta.
#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "LocalMarkerGraph2.hpp"
#include "LocalReadGraph.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Standard library.
#include "chrono.hpp"
#include <queue>



void Assembler::accessMarkerGraphVertices()
{
    globalMarkerGraphVertex.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphVertex"));

    globalMarkerGraphVertices.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphVertices"));
}



void Assembler::checkMarkerGraphVerticesAreAvailable()
{
    if(!globalMarkerGraphVertices.isOpen() || !globalMarkerGraphVertex.isOpen) {
        throw runtime_error("Vertices of the marker graph are not accessible.");
    }
}



// Find the vertex of the global marker graph that contains a given marker.
GlobalMarkerGraphVertexId Assembler::getGlobalMarkerGraphVertex(
    ReadId readId,
    Strand strand,
    uint32_t ordinal) const
{
    return getGlobalMarkerGraphVertex(OrientedReadId(readId, strand), ordinal);

}
GlobalMarkerGraphVertexId Assembler::getGlobalMarkerGraphVertex(
    OrientedReadId orientedReadId,
    uint32_t ordinal) const
{
    const MarkerId markerId =  getMarkerId(orientedReadId, ordinal);
    return globalMarkerGraphVertex[markerId];
}


// Find the markers contained in a given vertex of the global marker graph.
// Returns the markers as tuples(read id, strand, ordinal).
vector< tuple<ReadId, Strand, uint32_t> >
    Assembler::getGlobalMarkerGraphVertexMarkers(
        GlobalMarkerGraphVertexId globalMarkerGraphVertexId) const
{
    // Call the lower level function.
    vector< pair<OrientedReadId, uint32_t> > markers;
    getGlobalMarkerGraphVertexMarkers(globalMarkerGraphVertexId, markers);

    // Create the return vector.
    vector< tuple<ReadId, Strand, uint32_t> > returnVector;
    for(const auto& marker: markers) {
        const OrientedReadId orientedReadId = marker.first;
        const uint32_t ordinal = marker.second;
        returnVector.push_back(make_tuple(orientedReadId.getReadId(), orientedReadId.getStrand(), ordinal));
    }
    return returnVector;
}
void Assembler::getGlobalMarkerGraphVertexMarkers(
    GlobalMarkerGraphVertexId vertexId,
    vector< pair<OrientedReadId, uint32_t> >& markers) const
{
    markers.clear();
    for(const MarkerId markerId: globalMarkerGraphVertices[vertexId]) {
        OrientedReadId orientedReadId;
        uint32_t ordinal;
        tie(orientedReadId, ordinal) = findMarkerId(markerId);
        markers.push_back(make_pair(orientedReadId, ordinal));
    }
}



// Find the children of a vertex of the global marker graph.
vector<GlobalMarkerGraphVertexId>
    Assembler::getGlobalMarkerGraphVertexChildren(
    GlobalMarkerGraphVertexId vertexId) const
{
    vector<GlobalMarkerGraphVertexId> children;
    getGlobalMarkerGraphVertexChildren(vertexId, children);
    return children;
}
void Assembler::getGlobalMarkerGraphVertexChildren(
    GlobalMarkerGraphVertexId vertexId,
    vector<GlobalMarkerGraphVertexId>& children,
    bool append
    ) const
{

    if(!append) {
        children.clear();
    }
    if(isBadMarkerGraphVertex(vertexId)) {
        return;
    }

    // Loop over the markers of this vertex.
    for(const MarkerId markerId: globalMarkerGraphVertices[vertexId]) {

        // Find the OrientedReadId and ordinal.
        OrientedReadId orientedReadId;
        uint32_t ordinal;
        tie(orientedReadId, ordinal) = findMarkerId(markerId);

        // Find the next marker that is contained in a vertex.
        ++ordinal;
        for(; ordinal<markers.size(orientedReadId.getValue()); ++ordinal) {

            // Find the vertex id.
            const MarkerId childMarkerId =  getMarkerId(orientedReadId, ordinal);
            const GlobalMarkerGraphVertexId childVertexId =
                globalMarkerGraphVertex[childMarkerId];

            // If this marker correspond to a vertex, add it to our list.
            if(childVertexId != invalidCompressedGlobalMarkerGraphVertexId &&
                !isBadMarkerGraphVertex(childVertexId)) {
                children.push_back(childVertexId);
                break;
            }
        }

    }



    // Deduplicate.
    sort(children.begin(), children.end());
    children.resize(std::unique(children.begin(), children.end()) - children.begin());
}



// This version also returns the oriented read ids and ordinals
// that caused a child to be marked as such.
void Assembler::getGlobalMarkerGraphVertexChildren(
    GlobalMarkerGraphVertexId vertexId,
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > >& children,
    vector< pair<GlobalMarkerGraphVertexId, MarkerGraphNeighborInfo> >& workArea
    ) const
{
    children.clear();
    workArea.clear();

    if(isBadMarkerGraphVertex(vertexId)) {
        return;
    }

    // Loop over the markers of this vertex.
    for(const MarkerId markerId: globalMarkerGraphVertices[vertexId]) {

        // Find the OrientedReadId and ordinal.
        MarkerGraphNeighborInfo info;
        tie(info.orientedReadId, info.ordinal0) = findMarkerId(markerId);

        // Find the next marker that is contained in a vertex.
        const auto markerCount = markers.size(info.orientedReadId.getValue());
        for(info.ordinal1=info.ordinal0+1; info.ordinal1<markerCount; ++info.ordinal1) {

            // Find the vertex id.
            const MarkerId childMarkerId =  getMarkerId(info.orientedReadId, info.ordinal1);
            const GlobalMarkerGraphVertexId childVertexId =
                globalMarkerGraphVertex[childMarkerId];

            // If this marker correspond to a vertex, add it to our list.
            if( childVertexId!=invalidCompressedGlobalMarkerGraphVertexId &&
                !isBadMarkerGraphVertex(childVertexId)) {
                workArea.push_back(make_pair(childVertexId, info));
                break;
            }
        }

    }
    sort(workArea.begin(), workArea.end());



    // Now construct the children by gathering streaks of workArea entries
    // with the same child vertex id.
    for(auto streakBegin=workArea.begin(); streakBegin!=workArea.end(); ) {
        auto streakEnd = streakBegin + 1;
        for(;
            streakEnd!=workArea.end() && streakEnd->first==streakBegin->first;
            streakEnd++) {
        }
        children.resize(children.size() + 1);
        children.back().first = streakBegin->first;
        auto& v = children.back().second;
        for(auto it=streakBegin; it!=streakEnd; it++) {
            v.push_back(it->second);
        }

        // Process the next streak.
        streakBegin = streakEnd;
    }
}



// Find the parents of a vertex of the global marker graph.
vector<GlobalMarkerGraphVertexId>
    Assembler::getGlobalMarkerGraphVertexParents(
    GlobalMarkerGraphVertexId vertexId) const
{
    vector<GlobalMarkerGraphVertexId> parents;
    getGlobalMarkerGraphVertexParents(vertexId, parents);
    return parents;
}
void Assembler::getGlobalMarkerGraphVertexParents(
    GlobalMarkerGraphVertexId vertexId,
    vector<GlobalMarkerGraphVertexId>& parents,
    bool append
    ) const
{

    if(!append) {
        parents.clear();
    }
    if(isBadMarkerGraphVertex(vertexId)) {
        return;
    }

    // Loop over the markers of this vertex.
    for(const MarkerId markerId: globalMarkerGraphVertices[vertexId]) {

        // Find the OrientedReadId and ordinal.
        OrientedReadId orientedReadId;
        uint32_t ordinal;
        tie(orientedReadId, ordinal) = findMarkerId(markerId);

        // Find the previous marker that is contained in a vertex.
        if(ordinal == 0) {
            continue;
        }
        --ordinal;
        for(; ; --ordinal) {

            // Find the vertex id.
            const MarkerId parentMarkerId =  getMarkerId(orientedReadId, ordinal);
            const GlobalMarkerGraphVertexId parentVertexId =
                globalMarkerGraphVertex[parentMarkerId];

            // If this marker correspond to a vertex, add it to our list.
            if(parentVertexId != invalidCompressedGlobalMarkerGraphVertexId &&
                !isBadMarkerGraphVertex(parentVertexId)) {
                parents.push_back(parentVertexId);
                break;
            }

            if(ordinal == 0) {
                break;
            }
        }
    }

    // Deduplicate.
    sort(parents.begin(), parents.end());
    parents.resize(std::unique(parents.begin(), parents.end()) - parents.begin());
}



// This version also returns the oriented read ids and ordinals
// that caused a parent to be marked as such.
void Assembler::getGlobalMarkerGraphVertexParents(
    GlobalMarkerGraphVertexId vertexId,
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > >& parents,
    vector< pair<GlobalMarkerGraphVertexId, MarkerGraphNeighborInfo> >& workArea
    ) const
{
    parents.clear();
    workArea.clear();

    if(isBadMarkerGraphVertex(vertexId)) {
        return;
    }

    // Loop over the markers of this vertex.
    for(const MarkerId markerId: globalMarkerGraphVertices[vertexId]) {

        // Find the OrientedReadId and ordinal.
        MarkerGraphNeighborInfo info;
        tie(info.orientedReadId, info.ordinal0) = findMarkerId(markerId);
        if(info.ordinal0 == 0) {
            continue;
        }

        // Find the previous marker that is contained in a vertex.
        for(info.ordinal1=info.ordinal0-1; ; --info.ordinal1) {

            // Find the vertex id.
            const MarkerId parentMarkerId =  getMarkerId(info.orientedReadId, info.ordinal1);
            const GlobalMarkerGraphVertexId parentVertexId =
                globalMarkerGraphVertex[parentMarkerId];

            // If this marker correspond to a vertex, add it to our list.
            if( parentVertexId!=invalidCompressedGlobalMarkerGraphVertexId &&
                !isBadMarkerGraphVertex(parentVertexId)) {
                workArea.push_back(make_pair(parentVertexId, info));
                break;
            }

            if(info.ordinal1  == 0) {
                break;
            }
        }

    }
    sort(workArea.begin(), workArea.end());



    // Now construct the parents by gathering streaks of workArea entries
    // with the same child vertex id.
    for(auto streakBegin=workArea.begin(); streakBegin!=workArea.end(); ) {
        auto streakEnd = streakBegin + 1;
        for(;
            streakEnd!=workArea.end() && streakEnd->first==streakBegin->first;
            streakEnd++) {
        }
        parents.resize(parents.size() + 1);
        parents.back().first = streakBegin->first;
        auto& v = parents.back().second;
        for(auto it=streakBegin; it!=streakEnd; it++) {
            v.push_back(it->second);
        }

        // Process the next streak.
        streakBegin = streakEnd;
    }
}



// Python-callable function to get information about an edge of the
// global marker graph. Returns an empty vector if the specified
// edge does not exist.
vector<Assembler::GlobalMarkerGraphEdgeInformation> Assembler::getGlobalMarkerGraphEdgeInformation(
    GlobalMarkerGraphVertexId vertexId0,
    GlobalMarkerGraphVertexId vertexId1
    )
{
    const uint32_t k = uint32_t(assemblerInfo->k);

    // Find the children of vertexId0.
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > > children;
    vector< pair<GlobalMarkerGraphVertexId, MarkerGraphNeighborInfo> > workArea;
    getGlobalMarkerGraphVertexChildren(vertexId0, children, workArea);

    // Find vertexId1 in the children.
    vector<GlobalMarkerGraphEdgeInformation> v;
    for(const auto& child: children) {
        if(child.first != vertexId1) {
            continue;
        }

        // We found vertexId1. Loop over its MarkerGraphNeighborInfo.
        const auto& childInfos = child.second;
        v.resize(childInfos.size());
        for(size_t i=0; i<v.size(); i++) {
            const auto& childInfo = childInfos[i];
            auto& info = v[i];
            info.readId = childInfo.orientedReadId.getReadId();
            info.strand = childInfo.orientedReadId.getStrand();
            info.ordinal0 = childInfo.ordinal0;
            info.ordinal1 = childInfo.ordinal1;

            // Get the positions.
            const MarkerId markerId0 = getMarkerId(childInfo.orientedReadId, info.ordinal0);
            const MarkerId markerId1 = getMarkerId(childInfo.orientedReadId, info.ordinal1);
            const auto& marker0 = markers.begin()[markerId0];
            const auto& marker1 = markers.begin()[markerId1];
            info.position0 = marker0.position;
            info.position1 = marker1.position;

            // Construct the sequence.
            if(info.position1 <= info.position0+k) {
                // The marker overlap.
                info.overlappingBaseCount = info.position0+k - info.position1;
            } else {
                // The markers don't overlap.
                info.overlappingBaseCount = 0;
                for(uint32_t position=info.position0+k; position!=info.position1; position++) {
                    const Base base = getOrientedReadBase(childInfo.orientedReadId, position);
                    info.sequence.push_back(base.character());
                }
            }
        }
    }

    // If getting here, vertexId1 was not found, and we return
    // and empty veector.
    return v;
}



// Return true if a vertex of the global marker graph has more than
// one marker for at least one oriented read id.
bool Assembler::isBadMarkerGraphVertex(GlobalMarkerGraphVertexId vertexId) const
{
    // Get the markers of this vertex.
    const auto& vertexMarkerIds = globalMarkerGraphVertices[vertexId];

    // The markers are sorted by OrientedReadId, so we can just check each
    // consecutive pairs.
    for(size_t i=1; i<vertexMarkerIds.size(); i++) {
        const MarkerId markerId0 = vertexMarkerIds[i-1];
        const MarkerId markerId1 = vertexMarkerIds[i];
        OrientedReadId orientedReadId0;
        OrientedReadId orientedReadId1;
        tie(orientedReadId0, ignore) = findMarkerId(markerId0);
        tie(orientedReadId1, ignore) = findMarkerId(markerId1);
        if(orientedReadId0 == orientedReadId1) {
            return true;
        }
    }
    return false;
}



void Assembler::extractLocalMarkerGraph(

    // The ReadId, Strand, and ordinal that identify the
    // marker corresponding to the start vertex
    // for the local marker graph to be created.
    ReadId readId,
    Strand strand,
    uint32_t ordinal,

    // Maximum distance from the start vertex (number of edges in the global marker graph).
    int distance,

    // Minimum coverage for a strong vertex or edge (affects coloring).
    size_t minCoverage

    )
{
    // Create the local marker graph.
    LocalMarkerGraph2 graph(
        uint32_t(assemblerInfo->k),
        reads,
        assemblerInfo->useRunLengthReads,
        readRepeatCounts,
        markers,
        globalMarkerGraphVertex);
    extractLocalMarkerGraph(OrientedReadId(readId, strand), ordinal, distance, 0., graph);

    cout << "The local marker graph has " << num_vertices(graph);
    cout << " vertices and " << num_edges(graph) << " edges." << endl;

    // Write it out.
    graph.write("MarkerGraph.dot", minCoverage, distance, false, true);
    graph.write("DetailedMarkerGraph.dot", minCoverage, distance, true, true);

}



bool Assembler::extractLocalMarkerGraph(
    OrientedReadId orientedReadId,
    uint32_t ordinal,
    int distance,
    double timeout,                 // Or 0 for no timeout.
    LocalMarkerGraph2& graph
    )
{
    const GlobalMarkerGraphVertexId startVertexId =
        getGlobalMarkerGraphVertex(orientedReadId, ordinal);
    return extractLocalMarkerGraph(startVertexId, distance, timeout, graph);

}



bool Assembler::extractLocalMarkerGraph(
    GlobalMarkerGraphVertexId startVertexId,
    int distance,
    double timeout,                 // Or 0 for no timeout.
    LocalMarkerGraph2& graph
    )
{


    using vertex_descriptor = LocalMarkerGraph2::vertex_descriptor;
    using edge_descriptor = LocalMarkerGraph2::edge_descriptor;
    const auto startTime = steady_clock::now();

    // Add the start vertex.
    if(startVertexId == invalidCompressedGlobalMarkerGraphVertexId) {
        return true;    // Because no timeout occurred.
    }
    const vertex_descriptor vStart = graph.addVertex(startVertexId, 0, globalMarkerGraphVertices[startVertexId]);

    // Some vectors used inside the BFS.
    // Define them here to reduce memory allocation activity.
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > > children;
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > > parents;
    vector< pair<GlobalMarkerGraphVertexId, MarkerGraphNeighborInfo> > workArea;
    vector<LocalMarkerGraph2Edge::Info> infoVector;



    // Do the BFS.
    std::queue<vertex_descriptor> q;
    if(distance > 0) {
        q.push(vStart);
    }
    while(!q.empty()) {

        // See if we exceeded the timeout.
        if(timeout>0. && seconds(steady_clock::now() - startTime) > timeout) {
            graph.clear();
            return false;
        }

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const LocalMarkerGraph2Vertex& vertex0 = graph[v0];
        const GlobalMarkerGraphVertexId vertexId0 = vertex0.vertexId;
        const int distance0 = vertex0.distance;
        const int distance1 = distance0 + 1;

        // Get the children and parents.
        getGlobalMarkerGraphVertexChildren(vertexId0, children, workArea);
        getGlobalMarkerGraphVertexParents (vertexId0, parents, workArea);

        // Loop over the children.
        for(const auto& p: children) {
            const GlobalMarkerGraphVertexId vertexId1 = p.first;
            bool vertexExists;

            // Find the vertex corresponding to this child, creating it if necessary.
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);
            if(!vertexExists) {
                v1 = graph.addVertex(
                    vertexId1, distance1, globalMarkerGraphVertices[vertexId1]);
                if(distance1 < distance) {
                    q.push(v1);
                }
            }

            // Create the edge v0->v1, if it does not already exist.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            if(!edgeExists) {
                tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
                CZI_ASSERT(edgeExists);

                // Fill in edge information.
                infoVector.clear();
                const auto& v = p.second;
                for(const MarkerGraphNeighborInfo& x: v) {
                    infoVector.push_back(LocalMarkerGraph2Edge::Info(
                        x.orientedReadId, x.ordinal0, x.ordinal1));
                }
                graph.storeEdgeInfo(e, infoVector);
            }
        }


        // Loop over the parents.
        for(const auto& p: parents) {
            const GlobalMarkerGraphVertexId vertexId1 = p.first;
            bool vertexExists;

            // Find the vertex corresponding to this parent, creating it if necessary.
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);
            if(!vertexExists) {
                v1 = graph.addVertex(
                    vertexId1, distance1, globalMarkerGraphVertices[vertexId1]);
                if(distance1 < distance) {
                    q.push(v1);
                }
            }

            // Create the edge v1->v0, if it does not already exist.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v1, v0, graph);
            if(!edgeExists) {
                tie(e, edgeExists) = boost::add_edge(v1, v0, graph);
                CZI_ASSERT(edgeExists);

                // Fill in edge information.
                infoVector.clear();
                const auto& v = p.second;
                for(const MarkerGraphNeighborInfo& x: v) {
                    infoVector.push_back(LocalMarkerGraph2Edge::Info(
                        x.orientedReadId, x.ordinal1, x.ordinal0));
                }
                graph.storeEdgeInfo(e, infoVector);
            }
        }

    }



    // The BFS process did not create edges between vertices at maximum distance.
    // Do it now.
    // Loop over all vertices at maximum distance.
    BGL_FORALL_VERTICES(v0, graph, LocalMarkerGraph2) {
        const LocalMarkerGraph2Vertex& vertex0 = graph[v0];
        if(vertex0.distance != distance) {
            continue;
        }

        // Loop over the children that exist in the local marker graph
        // and are also at maximum distance.
        getGlobalMarkerGraphVertexChildren(vertex0.vertexId, children, workArea);
        for(const auto& p: children) {
            const GlobalMarkerGraphVertexId vertexId1 = p.first;
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);

            // If it does not exist in the local marker graph, skip.
            if(!vertexExists) {
                continue;
            }

            // If it is not at maximum distance, skip.
            const LocalMarkerGraph2Vertex& vertex1 = graph[v1];
            if(vertex1.distance != distance) {
                continue;
            }

            // There is no way we already created this edge.
            // Check that this is the case.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            CZI_ASSERT(!edgeExists);

            // Add the edge.
            tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
            CZI_ASSERT(edgeExists);

            // Fill in edge information.
            infoVector.clear();
            const auto& v = p.second;
            for(const MarkerGraphNeighborInfo& x: v) {
                infoVector.push_back(LocalMarkerGraph2Edge::Info(
                    x.orientedReadId, x.ordinal0, x.ordinal1));
            }
            graph.storeEdgeInfo(e, infoVector);
        }
    }

    // Fill in the oriented read ids represented in the graph.
    graph.findOrientedReadIds();

    // If using the run-length representation of reads,
    // also compute SeqAn alignments for all edges.
    // This is work in progress and not yet activated.
    // graph.computeSeqanAlignments();
    return true;
}



// Create a local marker graph and return its local assembly path.
// The local marker graph is specified by its start vertex
// and maximum distance (number of edges) form the start vertex.
vector<GlobalMarkerGraphVertexId> Assembler::getLocalAssemblyPath(
    GlobalMarkerGraphVertexId startVertexId,
    int maxDistance
    )
{

    // Create the local marker graph.
    LocalMarkerGraph2 graph(
        uint32_t(assemblerInfo->k),
        reads,
        assemblerInfo->useRunLengthReads,
        readRepeatCounts,
        markers,
        globalMarkerGraphVertex);
    extractLocalMarkerGraph(startVertexId, maxDistance, 0., graph);

    // Construct the local assembly path.
    graph.approximateTopologicalSort();
    graph.computeOptimalSpanningTree();
    graph.computeOptimalSpanningTreeBestPath();
    graph.computeLocalAssemblyPath(maxDistance);

    // Get the vertex ids in the assembly path.
    vector<GlobalMarkerGraphVertexId> path;
    if(!graph.localAssemblyPath.empty()) {
        LocalMarkerGraph2::edge_descriptor e = graph.localAssemblyPath.front();
        const LocalMarkerGraph2::vertex_descriptor v = source(e, graph);
        path.push_back(graph[v].vertexId);
        for(LocalMarkerGraph2::edge_descriptor e: graph.localAssemblyPath) {
            const LocalMarkerGraph2::vertex_descriptor v = target(e, graph);
            path.push_back(graph[v].vertexId);
        }
    }
    return path;
}



// Compute connectivity of the global marker graph.
// Vertices with more than markerCountOverflow are skipped.
void Assembler::createMarkerGraphConnectivity(
    size_t threadCount,
    size_t markerCountOverflow
    )
{
    cout << timestamp << "createMarkerGraphConnectivity begins." << endl;

    // Check that we have what we need.
    checkMarkerGraphVerticesAreAvailable();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Store markerCountOverflow so all threads can see it.
    markerGraphConnectivity.markerCountOverflow = markerCountOverflow;

    // Each thread stores the edges it finds in a separate vector.
    markerGraphConnectivity.threadEdges.resize(threadCount);
    setupLoadBalancing(globalMarkerGraphVertices.size(), 100000);
    runThreads(&Assembler::createMarkerGraphConnectivityThreadFunction0, threadCount,
        "threadLogs/createMarkerGraphConnectivity0");

    // Combine the edges found by each thread.
    markerGraphConnectivity.edges.createNew(
            largeDataName("GlobalMarkerGraphEdges"),
            largeDataPageSize);
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        auto& thisThreadEdges = *markerGraphConnectivity.threadEdges[threadId];
        for(const auto& edge: thisThreadEdges) {
            markerGraphConnectivity.edges.push_back(edge);
        }
        thisThreadEdges.remove();
    }
    cout << timestamp << "Found " << markerGraphConnectivity.edges.size();
    cout << " edges for " << globalMarkerGraphVertices.size() << " vertices." << endl;



    // Now we need to create edgesBySource and edgesByTarget.
    markerGraphConnectivity.edgesBySource.createNew(
        largeDataName("GlobalMarkerGraphEdgesBySource"),
        largeDataPageSize);
    markerGraphConnectivity.edgesByTarget.createNew(
        largeDataName("GlobalMarkerGraphEdgesByTarget"),
        largeDataPageSize);

    cout << timestamp << "Creating connectivity: pass 1 begins." << endl;
    markerGraphConnectivity.edgesBySource.beginPass1(globalMarkerGraphVertices.size());
    markerGraphConnectivity.edgesByTarget.beginPass1(globalMarkerGraphVertices.size());
    setupLoadBalancing(markerGraphConnectivity.edges.size(), 100000);
    runThreads(&Assembler::createMarkerGraphConnectivityThreadFunction1, threadCount);

    cout << timestamp << "Creating connectivity: pass 2 begins." << endl;
    markerGraphConnectivity.edgesBySource.beginPass2();
    markerGraphConnectivity.edgesByTarget.beginPass2();
    setupLoadBalancing(markerGraphConnectivity.edges.size(), 100000);
    runThreads(&Assembler::createMarkerGraphConnectivityThreadFunction2, threadCount);
    markerGraphConnectivity.edgesBySource.endPass2();
    markerGraphConnectivity.edgesByTarget.endPass2();

    cout << timestamp << "createMarkerGraphConnectivity ends." << endl;
}



void Assembler::createMarkerGraphConnectivityThreadFunction0(size_t threadId)
{
    ostream& out = getLog(threadId);

    // Create the vector to contain the edges found by this thread.
    using std::shared_ptr;
    using std::make_shared;
    shared_ptr< MemoryMapped::Vector<MarkerGraphConnectivity::Edge> > thisThreadEdgesPointer =
        make_shared< MemoryMapped::Vector<MarkerGraphConnectivity::Edge> >();
    markerGraphConnectivity.threadEdges[threadId] = thisThreadEdgesPointer;
    MemoryMapped::Vector<MarkerGraphConnectivity::Edge>& thisThreadEdges = *thisThreadEdgesPointer;
    thisThreadEdges.createNew(
            largeDataName("tmp-ThreadGlobalMarkerGraphEdges-" + to_string(threadId)),
            largeDataPageSize);

    // Some things used inside the loop but defined here for performance.
    vector<GlobalMarkerGraphVertexId> children;
    MarkerGraphConnectivity::Edge edge;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << begin << endl;

        // Loop over all marker graph vertices assigned to this batch.
        for(GlobalMarkerGraphVertexId vertex0=begin; vertex0!=end; ++vertex0) {
            // out << timestamp << vertex0 << " " << globalMarkerGraphVertices.size(vertex0) << endl;
            edge.source = vertex0;
            const auto markerIds0 = globalMarkerGraphVertices[vertex0];

            // Skip it is it has too many markers.
            if(markerIds0.size() > markerGraphConnectivity.markerCountOverflow) {
                continue;
            }

            // Loop over children of this vertex.
            getGlobalMarkerGraphVertexChildren(vertex0, children);
            for(const GlobalMarkerGraphVertexId vertex1: children) {
                edge.target = vertex1;
                const auto markerIds1 = globalMarkerGraphVertices[vertex1];

                // Skip it is it has too many markers.
                if(markerIds1.size() > markerGraphConnectivity.markerCountOverflow) {
                    continue;
                }

                // Joint loop over markers to compute coverage.
                // This code is similar to LocalMarkerGraph2::storeEdgeInfo.
                uint32_t coverage = 0;

                // Find pairs of markers for the same oriented read in the two vertices.
                // We exploit the fact that the markers in each
                // of the vertices are sorted.
                auto it0 = markerIds0.begin();
                auto it1 = markerIds1.begin();
                const auto end0 = markerIds0.end();
                const auto end1 = markerIds1.end();
                while(it0!=end0 && it1!=end1) {
                    const MarkerId markerId0 = *it0;
                    const MarkerId markerId1 = *it1;

                    // Find the oriented read ids and ordinals for these markers.
                    // This uses findMarkerId which requires binary searches
                    // and therefore could be expensive.
                    // If this becomes a performance problem we can
                    // store the OrientedReadId in each marker (4 extra bytes per marker).
                    OrientedReadId orientedReadId0;
                    uint32_t ordinal0;
                    tie(orientedReadId0, ordinal0) = findMarkerId(markerId0);
                    OrientedReadId orientedReadId1;
                    uint32_t ordinal1;
                    tie(orientedReadId1, ordinal1) = findMarkerId(markerId1);

                    if(orientedReadId0 < orientedReadId1) {
                        ++it0;
                        continue;
                    }
                    if(orientedReadId1 < orientedReadId0) {
                        ++it1;
                        continue;
                    }

                    // If getting here, the two oriented read ids are the same.
                    CZI_ASSERT(orientedReadId0 == orientedReadId1);
                    const OrientedReadId orientedReadId = orientedReadId0;

                    // Find the range of marker ids that correspond to this orientedReadId.
                    const auto thisOrientedReadMarkers = markers[orientedReadId.getValue()];
                    const MarkerId markerIdEnd   = thisOrientedReadMarkers.end()   - markers.begin();


                    // Find the streaks of markers for the same oriented readId.
                    auto it0StreakEnd = it0;
                    while(it0StreakEnd!=end0 && *it0StreakEnd<markerIdEnd) {
                        ++it0StreakEnd;
                    }
                    auto it1StreakEnd = it1;
                    while(it1StreakEnd!=end1 && *it1StreakEnd<markerIdEnd) {
                        ++it1StreakEnd;
                    }


                    // Only do it if both streaks contain one marker,
                    // the ordinal for the source vertex
                    // is less than the ordinal for the target vertex,
                    // and there are no intervening markers that also belong to a
                    // vertex of the marker graph.
                    if(it0StreakEnd-it0==1 && it1StreakEnd-it1==1 && ordinal0<ordinal1) {

                        // Check that there are no intervening markers that also belong to a
                        // vertex of the marker graph.
                        bool interveningVertexFound = false;
                        for(MarkerId markerId=markerId0+1; markerId!=markerId1; markerId++) {
                            if(globalMarkerGraphVertex[markerId] != invalidCompressedGlobalMarkerGraphVertexId) {
                                interveningVertexFound = true;
                                break;
                            }

                        }
                        if(!interveningVertexFound) {
                            ++coverage;
                        }
                    }

                    // Update the iterators to point to the end of the streaks.
                    it0 = it0StreakEnd;
                    it1 = it1StreakEnd;
                }


                // Store this edge.
                if(coverage >= 255) {
                    edge.coverage = 255;
                } else {
                    edge.coverage = uint8_t(coverage);
                }
                thisThreadEdges.push_back(edge);

            }

        }
    }

}



void Assembler::createMarkerGraphConnectivityThreadFunction1(size_t threadId)
{
    createMarkerGraphConnectivityThreadFunction12(threadId, 1);
}
void Assembler::createMarkerGraphConnectivityThreadFunction2(size_t threadId)
{
    createMarkerGraphConnectivityThreadFunction12(threadId, 2);
}
void Assembler::createMarkerGraphConnectivityThreadFunction12(size_t threadId, size_t pass)
{
    CZI_ASSERT(pass==1 || pass==2);

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all marker graph edges assigned to this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            const auto& edge = markerGraphConnectivity.edges[i];
            if(pass == 1) {
                markerGraphConnectivity.edgesBySource.incrementCountMultithreaded(edge.source);
                markerGraphConnectivity.edgesByTarget.incrementCountMultithreaded(edge.target);
            } else {
                markerGraphConnectivity.edgesBySource.storeMultithreaded(edge.source, Uint40(i));
                markerGraphConnectivity.edgesByTarget.storeMultithreaded(edge.target, Uint40(i));
            }
        }
    }

}


void Assembler::accessMarkerGraphConnectivity(bool accessEdgesReadWrite)
{
    if(accessEdgesReadWrite) {
        markerGraphConnectivity.edges.accessExistingReadWrite(
            largeDataName("GlobalMarkerGraphEdges"));
    } else {
        markerGraphConnectivity.edges.accessExistingReadOnly(
            largeDataName("GlobalMarkerGraphEdges"));
    }
    markerGraphConnectivity.edgesBySource.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphEdgesBySource"));
    markerGraphConnectivity.edgesByTarget.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphEdgesByTarget"));
}



// Locate the edge given the vertices.
const Assembler::MarkerGraphConnectivity::Edge*
    Assembler::MarkerGraphConnectivity::findEdge(Uint40 source, Uint40 target) const
{
    const auto edgesWithThisSource = edgesBySource[source];
    for(const uint64_t i: edgesWithThisSource) {
        const Edge& edge = edges[i];
        if(edge.target == target) {
            return &edge;
        }
    }
    return 0;

}



// Flag as not good a marker graph edge if:
// - It has coverage<minCoverage, AND
// - A path of length <= maxPathLength edges exists that:
//    * Starts at the source vertex of the edge.
//    * Ends at the target vertex of the edge.
//    * Only uses edges with coverage>=minCoverage.
void Assembler::flagMarkerGraphEdges(
    size_t threadCount,
    size_t minCoverage,
    size_t maxPathLength)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Store the parameters so all threads can see them.
    flagMarkerGraphEdgesData.minCoverage = minCoverage;
    flagMarkerGraphEdgesData.maxPathLength = maxPathLength;

    // Do it in parallel.
    const size_t vertexCount = markerGraphConnectivity.edgesBySource.size();
    setupLoadBalancing(vertexCount, 100000);
    runThreads(
        &Assembler::flagMarkerGraphEdgesThreadFunction,
        threadCount,
        "threadLogs/flagMarkerGraphEdges");
}



void Assembler::flagMarkerGraphEdgesThreadFunction(size_t threadId)
{
    ostream& out = getLog(threadId);

    const size_t minCoverage = flagMarkerGraphEdgesData.minCoverage;
    const size_t maxPathLength = flagMarkerGraphEdgesData.maxPathLength;

    vector< vector<GlobalMarkerGraphVertexId> > verticesByDistance(maxPathLength+1);
    verticesByDistance.front().resize(1);
    vector<GlobalMarkerGraphVertexId> vertices;

    // Loop over all batches assigned to this thread.
    size_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << begin << endl;

        // Loop over vertices assigned to this batch.
        for(GlobalMarkerGraphVertexId startVertexId=begin; startVertexId!=end; ++startVertexId) {

            // Find all vertices within maxPathLength of this vertex,
            // computed using only edges with coverage >= minCoverage.
            // This uses a very simple algorithm that should work well
            // because the number of vertices involved is usually very small.
            verticesByDistance.front().front() = startVertexId;
            vertices.clear();
            vertices.push_back(startVertexId);
            for(size_t distance=1; distance<=maxPathLength; distance++) {
                auto& verticesAtThisDistance = verticesByDistance[distance];
                verticesAtThisDistance.clear();
                for(const GlobalMarkerGraphVertexId vertexId0: verticesByDistance[distance-1]) {
                    for(const uint64_t i: markerGraphConnectivity.edgesBySource[vertexId0]) {
                        const auto& edge = markerGraphConnectivity.edges[i];
                        if(edge.coverage < minCoverage) {
                            continue;
                        }
                        CZI_ASSERT(edge.source == vertexId0);
                        const GlobalMarkerGraphVertexId vertexId1 = edge.target;
                        verticesAtThisDistance.push_back(vertexId1);
                        vertices.push_back(vertexId1);
                    }
                }
            }
            sort(vertices.begin(), vertices.end());
            vertices.resize(unique(vertices.begin(), vertices.end()) - vertices.begin());

            // Loop over low coverage edges with source startVertexId.
            // If the target is one of the vertices we found, flag the edge as not good.
            for(const uint64_t i: markerGraphConnectivity.edgesBySource[startVertexId]) {
                auto& edge = markerGraphConnectivity.edges[i];
                if(edge.coverage >= minCoverage) {
                    continue;
                }
                CZI_ASSERT(edge.source == startVertexId);
                if(std::binary_search(vertices.begin(), vertices.end(), edge.target)) {
                    edge.isGood = false;
                }
            }

        }

    }

}

