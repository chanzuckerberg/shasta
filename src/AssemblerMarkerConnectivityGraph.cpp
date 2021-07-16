// Shasta.
#include "Assembler.hpp"
#include "MarkerConnectivityGraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <map>
#include <queue>



// Create the marker connectivity graph starting with a given marker.
void Assembler::createMarkerConnectivityGraph(
    OrientedReadId orientedReadId,
    uint32_t ordinal,
    bool useReadGraphAlignmentsOnly,
    MarkerConnectivityGraph& graph) const
{
    MarkerConnectivityGraphVertexMap vertexMap;
    createMarkerConnectivityGraph(orientedReadId, ordinal, useReadGraphAlignmentsOnly,
        graph, vertexMap);

}



void Assembler::createMarkerConnectivityGraph(
    OrientedReadId orientedReadId,
    uint32_t ordinal,
    bool useReadGraphAlignmentsOnly,
    MarkerConnectivityGraph& graph,
    MarkerConnectivityGraphVertexMap& vertexMap) const
{


    using vertex_descriptor = MarkerConnectivityGraph::vertex_descriptor;
    vertexMap.clear();

    // Initialize a BFS in the space of aligned markers.
    const vertex_descriptor v = add_vertex(MarkerDescriptor(orientedReadId, ordinal), graph);
    vertexMap[MarkerDescriptor(orientedReadId, ordinal)] = v;
    std::queue<vertex_descriptor> q;
    q.push(v);



    // BFS loop.
    vector<MarkerDescriptor> alignedMarkers;
    while(not q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const MarkerDescriptor markerDescriptor0 =graph[v0];
        SHASTA_ASSERT(vertexMap[markerDescriptor0] == v0);

        // Find aligned markers for this vertex.
        findAlignedMarkers(markerDescriptor0.first, markerDescriptor0.second,
            useReadGraphAlignmentsOnly, alignedMarkers);

        for(const MarkerDescriptor& markerDescriptor1: alignedMarkers) {

            // If there is no vertex for markerDescriptor1, add it and enqueue it.
            auto it1 = vertexMap.find(markerDescriptor1);
            if(it1 == vertexMap.end()) {
                const vertex_descriptor v1 = add_vertex(markerDescriptor1, graph);
                tie(it1, ignore) = vertexMap.insert(make_pair(markerDescriptor1, v1));
                q.push(v1);
            }
            const vertex_descriptor v1 = it1->second;

            // If there is no edge between v0 and v1, create one.
            bool edgeExists;
            tie(ignore, edgeExists) = edge(v0, v1, graph);
            if(not edgeExists) {
                add_edge(v0, v1, graph);
            }
        }
    }

}
