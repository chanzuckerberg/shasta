#ifndef SHASTA_MODE3_LOCAL_ASSEMBLY_GRAPH_HPP
#define SHASTA_MODE3_LOCAL_ASSEMBLY_GRAPH_HPP

#include "mode3.hpp"

namespace shasta {
    namespace mode3 {

        class LocalAssemblyGraph;
        class LocalAssemblyGraphEdge;
        class LocalAssemblyGraphVertex;
    }

}


// Classes used to display in the http server a local portion of the AssemblyGraph.
class shasta::mode3::LocalAssemblyGraphVertex {
public:
    uint64_t segmentId;
    uint64_t distance;  // From the start vertex.
    vector<MarkerGraphEdgeInfo> path;
    array<double, 2> position;
    LocalAssemblyGraphVertex(
        uint64_t segmentId,
        uint64_t distance,
        const span<const MarkerGraphEdgeInfo> path);
    LocalAssemblyGraphVertex();
};



class shasta::mode3::LocalAssemblyGraphEdge {
public:
    uint64_t coverage;
    LocalAssemblyGraphEdge(uint64_t coverage = 0) : coverage(coverage) {}
};



class shasta::mode3::LocalAssemblyGraph :
    public boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
    LocalAssemblyGraphVertex, LocalAssemblyGraphEdge> {
public:

    LocalAssemblyGraph(
        const AssemblyGraph&,
        uint64_t startSegmentId,
        uint64_t maxDistance);
    uint64_t maxDistance;

    vertex_descriptor addVertex(
        uint64_t segmentId,
        uint64_t distance,
        const span<const MarkerGraphEdgeInfo> path);

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

    void writeSvg1(const string& fileName, uint64_t sizePixels);
    void writeSvg1(ostream&, uint64_t sizePixels);
};
#endif

