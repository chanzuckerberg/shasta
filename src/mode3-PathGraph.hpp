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
    // (e. g. if ti goes around in a cycle).
    vector<AssemblyGraphJourneyInterval> journeyIntervals;

    // The vertex id is only used to help keep track of vertices
    // for testing and debugging.
    uint64_t id;
};



class shasta::mode3::PathGraphEdge {
public:
};



class shasta::mode3::PathGraph :
    public PathGraphBaseClass,
    public MultithreadedObject<PathGraph> {
public:

    // Create the PathGraph from the AssemblyGraph.
    PathGraph(const AssemblyGraph&);

    // Initial creation of the vertices.
    // Start with a single segment for each vertex
    // (that is, paths of length 1).
    void createVertices();

    const AssemblyGraph& assemblyGraph;

    uint64_t nextVertexId = 0;
};

#endif
