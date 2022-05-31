#include "mode3-PathGraph.hpp"
using namespace shasta;
using namespace mode3;

#include "iostream.hpp"



// Create the PathGraph from the AssemblyGraph.
// Start with a single segment for each vertex
// (that is, paths of length 1).
PathGraph::PathGraph(const AssemblyGraph& assemblyGraph) :
    MultithreadedObject<PathGraph>(*this),
    assemblyGraph(assemblyGraph)
{
    PathGraph& pathGraph = *this;

    createVertices();
    cout << "The initial path graph has " << num_vertices(pathGraph) <<
        " vertices." << endl;
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
