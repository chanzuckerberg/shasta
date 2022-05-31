// Shasta.
#include "mode3-PathGraph.hpp"
#include "orderPairs.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "iostream.hpp"



// Create the PathGraph from the AssemblyGraph.
// Start with a single segment for each vertex
// (that is, paths of length 1).
PathGraph::PathGraph(const AssemblyGraph& assemblyGraph) :
    MultithreadedObject<PathGraph>(*this),
    assemblyGraph(assemblyGraph)
{
    const uint64_t minCoverage = 3; // EXPOSE WHEN CODE STABILIZES

    PathGraph& pathGraph = *this;

    createVertices();
    createEdges(minCoverage);
    cout << "The initial path graph has " << num_vertices(pathGraph) <<
        " vertices and " << num_edges(pathGraph) << " edges." << endl;
}



// Initial creation of the vertices.
// Start with a single segment for each vertex
// (that is, paths of length 1).
void PathGraph::createVertices() {

    PathGraph& pathGraph = *this;


    // Create a vertex for each segment in the AssemblyGraph.
    for(uint64_t segmentId=0; segmentId<assemblyGraph.paths.size(); segmentId++) {

        // Create the vertex.
        const vertex_descriptor v = add_vertex(pathGraph);
        PathGraphVertex& vertex = pathGraph[v];
        vertex.id = nextVertexId++;

        // Store the path.
        vertex.path.push_back(segmentId);

        // Store the AssemblyGraphJourneyInterval's.
        const span<const pair<OrientedReadId, uint64_t> > journeyInfos =
            assemblyGraph.assemblyGraphJourneyInfos[segmentId];
        for(const pair<OrientedReadId, uint64_t>& p: journeyInfos) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t position = p.second;
            AssemblyGraphJourneyInterval interval;
            interval.orientedReadId = orientedReadId;
            interval.first = position;
            interval.last = position;
            vertex.journeyIntervals.push_back(interval);
        }
    }

}



// Recreate all edges from scratch, using only the
// information stored in the vertices.
void PathGraph::createEdges(uint64_t minCoverage)
{
    PathGraph& pathGraph = *this;

    // Gather AssemblyGraphJourneyInterval's for all oriented reads.
    vector< vector<pair<AssemblyGraphJourneyInterval, vertex_descriptor> > >
        journeyIntervals(2 * assemblyGraph.readCount());
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        for(const AssemblyGraphJourneyInterval& interval: pathGraph[v].journeyIntervals) {
            journeyIntervals[interval.orientedReadId.getValue()].push_back(
                make_pair(interval, v));
        }
    }
    for(auto& v: journeyIntervals) {
        sort(v.begin(), v.end(),
            OrderPairsByFirstOnly<AssemblyGraphJourneyInterval, vertex_descriptor>());
    }


    // Create the edges.
    for(const auto& orientedReadJourneyIntervals: journeyIntervals) {

        for(uint64_t i=1; i<orientedReadJourneyIntervals.size(); i++) {
            const vertex_descriptor v0 = orientedReadJourneyIntervals[i-1].second;
            const vertex_descriptor v1 = orientedReadJourneyIntervals[i  ].second;

            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(v0, v1, pathGraph);
            if(not edgeExists) {
                tie(e, edgeExists) = add_edge(v0, v1, pathGraph);
                SHASTA_ASSERT(edgeExists);
            }
            ++pathGraph[e].coverage;
        }
    }



    // Remove the low coverage edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        if(pathGraph[e].coverage < minCoverage) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, pathGraph);
    }
}

