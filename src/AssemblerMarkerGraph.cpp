// shasta.
#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "LocalMarkerGraph.hpp"
#include "LocalMarkerGraph2.hpp"
#include "LocalReadGraph.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Standard library.
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
    bool append) const
{
    if(!append) {
        children.clear();
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
            if(childVertexId != invalidCompressedGlobalMarkerGraphVertexId) {
                children.push_back(childVertexId);
                break;
            }
        }

    }



    // Deduplicate.
    sort(children.begin(), children.end());
    children.resize(std::unique(children.begin(), children.end()) - children.begin());
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
    bool append) const
{
    if(!append) {
        parents.clear();
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
            if(parentVertexId != invalidCompressedGlobalMarkerGraphVertexId) {
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



void Assembler::extractLocalMarkerGraph(

    // The ReadId, Strand, and ordinal that identify the
    // marker corresponding to the start vertex
    // for the local marker graph to be created.
    ReadId readId,
    Strand strand,
    uint32_t ordinal,

    // Maximum distance from the start vertex (number of edges in the global marker graph).
    int distance,

    // Minimum coverage for a strong vertex.
    size_t minCoverage,

    // Minimum consensus for a strong edge.
    size_t minConsensus
    )
{
    // Create the local marker graph.
    LocalMarkerGraph2 graph(uint32_t(assemblerInfo->k), reads, markers, globalMarkerGraphVertex);
    extractLocalMarkerGraph(OrientedReadId(readId, strand), ordinal, distance, graph);

    cout << "The local marker graph has " << num_vertices(graph);
    cout << " vertices and " << num_edges(graph) << " edges." << endl;

    // Write it out.
    graph.write("MarkerGraph.dot", minCoverage, distance, false, true);
    graph.write("DetailedMarkerGraph.dot", minCoverage, distance, true, true);

}



void Assembler::extractLocalMarkerGraph(
    OrientedReadId orientedReadId,
    uint32_t ordinal,
    int distance,
    LocalMarkerGraph2& graph
    )
{
    using vertex_descriptor = LocalMarkerGraph2::vertex_descriptor;

    // Add the start vertex.
    const GlobalMarkerGraphVertexId startVertexId =
        getGlobalMarkerGraphVertex(orientedReadId, ordinal);
    if(startVertexId == invalidCompressedGlobalMarkerGraphVertexId) {
        return;
    }
    const vertex_descriptor vStart = graph.addVertex(startVertexId, 0, globalMarkerGraphVertices[startVertexId]);

    // Do the BFS. Do not add the edges now.
    // We will add the edges later.
    std::queue<vertex_descriptor> q;
    if(distance > 0) {
        q.push(vStart);
    }
    vector<GlobalMarkerGraphVertexId> neighbors;
    while(!q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const LocalMarkerGraph2Vertex& vertex0 = graph[v0];
        const GlobalMarkerGraphVertexId vertexId0 = vertex0.vertexId;
        const int distance0 = vertex0.distance;
        const int distance1 = distance0 + 1;

        // Get the neighbors.
        neighbors.clear();
        getGlobalMarkerGraphVertexChildren(vertexId0, neighbors, true);
        getGlobalMarkerGraphVertexParents (vertexId0, neighbors, true);

        // Loop over the neighbors.
        for(const GlobalMarkerGraphVertexId vertexId1: neighbors) {
            bool vertexExists;
            tie(vertexExists, ignore) = graph.findVertex(vertexId1);
            if(!vertexExists) {
                const vertex_descriptor v1 = graph.addVertex(
                    vertexId1, distance1, globalMarkerGraphVertices[vertexId1]);
                if(distance1 < distance) {
                    q.push(v1);
                }
            }
        }
    }


    // Now we can add the edges.
    BGL_FORALL_VERTICES(v0, graph, LocalMarkerGraph2) {
        getGlobalMarkerGraphVertexChildren(graph[v0].vertexId, neighbors);
        for(const GlobalMarkerGraphVertexId vertexId1: neighbors) {
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);
            if(vertexExists) {
                LocalMarkerGraph2::edge_descriptor  e;
                bool edgeWasAdded = false;
                tie(e, edgeWasAdded) = add_edge(v0, v1, graph);
                CZI_ASSERT(edgeWasAdded);
                graph.storeEdgeInfo(e);
            }
        }
    }

}



// Connectivity of the global marker graph.
void Assembler::createMarkerGraphConnectivity(size_t threadCount)
{
    cout << timestamp << "createMarkerGraphConnectivity begins." << endl;

    // Check that we have what we need.
    checkMarkerGraphVerticesAreAvailable();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Each thread stores the edges it finds in a separate vector.
    markerGraphConnectivity.threadEdges.resize(threadCount);
    setupLoadBalancing(globalMarkerGraphVertices.size(), 100000);
    runThreads(&Assembler::createMarkerGraphConnectivityThreadFunction0, threadCount);

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

        // Loop over all marker graph vertices assigned to this batch.
        for(GlobalMarkerGraphVertexId vertex0=begin; vertex0!=end; ++vertex0) {
            edge.source = vertex0;
            const auto markerIds0 = globalMarkerGraphVertices[vertex0];

            // Loop over children of this vertex.
            getGlobalMarkerGraphVertexChildren(vertex0, children);
            for(const GlobalMarkerGraphVertexId vertex1: children) {
                edge.target = vertex1;
                const auto markerIds1 = globalMarkerGraphVertices[vertex1];

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


void Assembler::accessMarkerGraphConnectivity()
{
    CZI_ASSERT(0);
}
