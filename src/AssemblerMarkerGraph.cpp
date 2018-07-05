// shasta.
#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "LocalMarkerGraph.hpp"
#include "LocalMarkerGraph2.hpp"
#include "LocalReadGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Standard library.
#include <queue>



void Assembler::accessGlobalMarkerGraph()
{
    globalMarkerGraphVertex.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphVertex"));

    globalMarkerGraphVertices.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphVertices"));
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
    graph.write("MarkerGraph.dot", minCoverage, minConsensus, distance, false);
    graph.write("DetailedMarkerGraph.dot", minCoverage, minConsensus, distance, true);

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

