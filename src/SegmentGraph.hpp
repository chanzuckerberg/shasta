#ifndef SHASTA_READ_LOADER_HPP
#define SHASTA_READ_LOADER_HPP

// The segment graph is used to detangle and analyze reachability
// in the assembly graph.
// See AssemblerSegmentGraph.cpp for more information.


#include "AssemblyGraph.hpp"
#include <boost/graph/adjacency_list.hpp>

namespace shasta {
    class SegmentGraph;
    class SegmentGraphVertex;
    class SegmentGraphEdge;

    using SegmentGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        SegmentGraphVertex,
        SegmentGraphEdge
        >;

}



class shasta::SegmentGraphVertex {
public:

    // The id of the key assembly graph edge (segment)
    // corresponding to this vertex of the segment graph.
    AssemblyGraph::EdgeId assemblyGraphEdgeId;

    SegmentGraphVertex(AssemblyGraph::EdgeId assemblyGraphEdgeId) :
        assemblyGraphEdgeId(assemblyGraphEdgeId) {}
};



class shasta::SegmentGraphEdge {
public:

    // The number of oriented reads supporting this edge.
    uint64_t coverage = 0;
};



class shasta::SegmentGraph : public SegmentGraphBaseClass {
public:

    void removeLowCoverageEdges(uint64_t minCoverage);

    // Approximate transitive reduction.
    void transitiveReduction(uint64_t maxDistance);


    void writeGraphviz(const string& fileName) const;

    // Find chains in the segment graph.
    // This assumes that no vertex has in-degree or out-degree
    // greater than one.
    void findChains();
    class Chain {
    public:
        vector<vertex_descriptor> vertices;
        bool isCircular;
    };
    vector<Chain> chains;

};

#endif
