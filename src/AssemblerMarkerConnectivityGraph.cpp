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
    using vertex_descriptor = MarkerConnectivityGraph::vertex_descriptor;
    std::map<MarkerPair, vertex_descriptor> vertexMap;

    // Initialize a BFS in the space of aligned markers.
    const vertex_descriptor v = add_vertex(MarkerPair(orientedReadId, ordinal), graph);
    vertexMap[MarkerPair(orientedReadId, ordinal)] = v;
    std::queue<vertex_descriptor> q;
    q.push(v);



    // BFS loop.
    vector<MarkerPair> alignedMarkers;
    while(not q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const MarkerPair markerPair0 =graph[v0];
        SHASTA_ASSERT(vertexMap[markerPair0] == v0);

        // Find aligned markers for this vertex.
        findAlignedMarkers(markerPair0.first, markerPair0.second,
            useReadGraphAlignmentsOnly, alignedMarkers);

        for(const MarkerPair& markerPair1: alignedMarkers) {

            // If there is no vertex for markerPair1, add it and enqueue it.
            auto it1 = vertexMap.find(markerPair1);
            if(it1 == vertexMap.end()) {
                const vertex_descriptor v1 = add_vertex(markerPair1, graph);
                tie(it1, ignore) = vertexMap.insert(make_pair(markerPair1, v1));
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
