#ifndef SHASTA_MODE3_HPP
#define SHASTA_MODE3_HPP

// Classes used for Mode 3 assembly.

// Shasta.
#include "MemoryMappedVectorOfVectors.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "array.hpp"
#include "vector.hpp"



namespace shasta {
    namespace mode3 {

        class DynamicAssemblyGraph;
        class DynamicAssemblyGraphVertex;
        class DynamicAssemblyGraphEdge;
        using DynamicAssemblyGraphBaseClass =
            boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
            DynamicAssemblyGraphVertex, DynamicAssemblyGraphEdge>;

        class AssemblyGraph;

        class VirtualMarkerGraphEdge;
    }

    class CompressedMarker;
    class MarkerGraph;
    class ReadFlags;
}



// A VirtualMarkerGraphEdge described a marker graph edge
// that actually does not exist in the marker graph.
class shasta::mode3::VirtualMarkerGraphEdge {
public:
    array<MarkerGraphVertexId, 2> vertices;
};



// Each vertex corresponds to a gfa Segment of the assembly graph -
// a linear sequence of marker graph edges, possibly containing
// "virtual marker graph edges".
class shasta::mode3::DynamicAssemblyGraphVertex {
public:

    DynamicAssemblyGraphVertex(const vector<MarkerGraphEdgeId>&);

    // A small class that can describe both a real and a virtual
    // marker graph edge.
    // If virtual, the edgeId is an index into the
    // virtualMarkerGraphEdges vector.
    class MarkerGraphEdgeInfo {
    public:
        uint64_t isVirtual : 1;
        MarkerGraphEdgeId edgeId: 63;
        MarkerGraphEdgeInfo(MarkerGraphEdgeId, bool isVirtual);
    };

    // The marker graph edges (real or virtual) that
    // make up the gfa Segment corresponding to this vertex.
    // The target vertex of each marker graph edge (real or virtual)
    // is always equal to the source vertex of the next marker graph edge.
    vector<MarkerGraphEdgeInfo> markerGraphEdges;
};



// Each edge corresponds to a gfa Link in the assembly graph.
class shasta::mode3::DynamicAssemblyGraphEdge {
public:
};



// The DynamicAssemblyGraph is used for initial creation
// of the Mode 3 assembly graph. It is implemented using
// the Boost Graph library.
// Note that the roles of vertices and edges are swapped
// relative to what happens in the shasta::AssemnblyGraph,
// where eahc gfa segment corresponds to an edge.
class shasta::mode3::DynamicAssemblyGraph :
    public DynamicAssemblyGraphBaseClass {
public:

    DynamicAssemblyGraph(
        const MemoryMapped::Vector<ReadFlags>&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&);
    void createVertices(
        const MemoryMapped::Vector<ReadFlags>&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&);

    vector<VirtualMarkerGraphEdge> virtualMarkerGraphEdges;
};



// The AssemblyGraph is used to store the Mode 3 assembly graph,
// when it no longer needs to be changed,
// in memory mapped data structures.
class shasta::mode3::AssemblyGraph {
public:
};

#endif

