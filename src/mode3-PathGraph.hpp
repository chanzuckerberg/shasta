#ifndef SHASTA_MODE3_PATH_GRAPH_HPP
#define SHASTA_MODE3_PATH_GRAPH_HPP

/*******************************************************************************

The mode3::PathGraph is a directed graph in which each vertex represents
a path in the mode3::AssemblyGraph.

*******************************************************************************/

// Shasta.
#include "mode3.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standara libraries.
#include <limits>
#include "vector.hpp"

namespace shasta {
    namespace mode3 {
        class PathGraph;
        class PathGraphVertex;
        class PathGraphEdge;

        using PathGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            PathGraphVertex, PathGraphEdge>;

    }
}



// Each vertex of the PathGraph describes a path
// in the mode3::AssemblyGraph.
class shasta::mode3::PathGraphVertex {
public:

    // The segment ids of the mode3::AssemblyGraph path
    // that this vertex describes.
    vector<uint64_t> path;

    // We also store the portions of the assembly graph journeys
    // that visit this path.
    // Note that an oriented read can have more than one
    // (e. g. if it goes around in a cycle).
    vector<AssemblyGraphJourneyInterval> journeyIntervals;

    // The vertex id is only used to help keep track of vertices
    // for testing and debugging.
    uint64_t id;

    // The partition this vertex was assigned to.
    uint64_t subgraphId = std::numeric_limits<uint64_t>::max();

    // Distance from the start vertex of the BFS.
    // Only used during the BFS.
    uint64_t distance = 0;
};



class shasta::mode3::PathGraphEdge {
public:
    uint64_t coverage = 0;
};



class shasta::mode3::PathGraph :
    public PathGraphBaseClass,
    public MultithreadedObject<PathGraph> {
public:

    // Create the PathGraph from the AssemblyGraph.
    PathGraph(const AssemblyGraph&);

    void writeGfa(const string& baseName) const;

private:

    // The AssemblyGraph this PathGraph refers to.
    const AssemblyGraph& assemblyGraph;

    // Initial creation of the vertices.
    // Start with a single segment for each vertex
    // (that is, paths of length 1).
    void createVertices();

    // Recreate all edges from scratch, using only the
    // information stored in the vertices.
    void createEdges(uint64_t minCoverage);

    // The id of the next vertex to be added.
    // Vertex ids are only used to help keep track of vertices
    // for testing and debugging.
    uint64_t nextVertexId = 0;

    // Partition the PathGraph into subgraphs.
    void partition(
        uint64_t maxDistance,
        uint64_t minSubgraphSize);
    static const uint64_t noSubgraph = std::numeric_limits<uint64_t>::max();

    // Gather subgraphs using the subgraphId stored in each vertex.
    // A subgraph can have size 0, and in that case it should be ignored.
    void gatherSubgraphs();
    void histogramSubgraphs();
    vector< vector<vertex_descriptor> > subgraphs;

    // A partition iteration does a single BFS starting at v.
    // It moves forward from v, avoiding vertices already
    // assigned to a subgraph, and up to maxDistance from v.
    // It also returns the boundaryVertices, that is the
    // vertices found in the process that are at distance maxDistance+1
    // from v and are nto yet assigned to a subgraph.
    // These can then used as starting points new partition iterations.
    void partitionIteration(
        vertex_descriptor v,
        uint64_t maxDistance,
        uint64_t subgraphId,
        vector<vertex_descriptor>& boundaryVertices);
};

#endif
