#ifndef SHASTA_MODE3_HPP
#define SHASTA_MODE3_HPP

// Classes used for Mode 3 assembly.

// Shasta.
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"
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
        class MarkerGraphEdgeInfo;
        class PseudoPathEntry;
        class Transition;
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



// An entry of the pseudo-path of an oriented read in the DynamicAssemblyGraph.
class shasta::mode3::PseudoPathEntry {
public:
    DynamicAssemblyGraphBaseClass::vertex_descriptor v;
    uint32_t position;
    array<uint32_t, 2> ordinals;

    bool operator<(const PseudoPathEntry& that) const
    {
        return ordinals[0] < that.ordinals[0];
    }
};



// A Transition occurs when the pseudopath of an oriented read
// moves from a Segment to a different segment.
// Transitions are used to create edges (gfa links).
class shasta::mode3::Transition : public array<PseudoPathEntry, 2> {
public:
    Transition(const array<PseudoPathEntry, 2>& x) : array<PseudoPathEntry, 2>(x) {}
};



// A small class that can describe both a real and a virtual
// marker graph edge.
// If virtual, the edgeId is an index into the
// virtualMarkerGraphEdges vector.
class shasta::mode3::MarkerGraphEdgeInfo {
public:
    uint64_t isVirtual : 1;
    MarkerGraphEdgeId edgeId: 63;
    MarkerGraphEdgeInfo(MarkerGraphEdgeId, bool isVirtual);
};



// Each vertex corresponds to a gfa Segment of the assembly graph -
// a linear sequence of marker graph edges, possibly containing
// "virtual marker graph edges".
class shasta::mode3::DynamicAssemblyGraphVertex {
public:

    DynamicAssemblyGraphVertex(const vector<MarkerGraphEdgeId>&, uint64_t vertexId);

    uint64_t vertexId;

    // The marker graph path that
    // makes up the gfa Segment corresponding to this vertex.
    // This can contain virtual marker graph edges.
    vector<MarkerGraphEdgeInfo> path;

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
    public DynamicAssemblyGraphBaseClass,
    public MultithreadedObject<DynamicAssemblyGraph> {
public:
    using G = DynamicAssemblyGraph;

    DynamicAssemblyGraph(
        const MemoryMapped::Vector<ReadFlags>&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&,
        const string& largeDataFileNamePrefix,
        size_t largeDataPageSize,
        size_t threadCount);
    ~DynamicAssemblyGraph();

    // Store some information passed in to the constructor.
    const MemoryMapped::Vector<ReadFlags>& readFlags;
    const MarkerGraph& markerGraph;
    const string& largeDataFileNamePrefix;
    size_t largeDataPageSize;
    size_t threadCount;

    uint64_t nextVertexId = 0;

    // Initial creation of vertices of the DynamicAssemblyGraph.
    // Each vertex (gfa segment) corresponds to a linear sequence
    // of edges (a path) in the marker graph.
    void createVertices(
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers);

    // For each marker graph edge, store in the marker graph edge table
    // the corresponding DynamicAssemblyGraph vertex (segment)
    // and position in the path, if any.
    // This is needed when computing pseudopaths.
    MemoryMapped::Vector< pair<vertex_descriptor, uint32_t> > markerGraphEdgeTable;
    void computeMarkerGraphEdgeTable();
    void computeMarkerGraphEdgeTableThreadFunction(size_t threadId);
    class MarkerGraphEdgeTableData {
    public:
        vector<vertex_descriptor> allVertices;
    };
    MarkerGraphEdgeTableData markerGraphEdgeTableData;


    // Compute pseudopaths for all oriented reads.
    // The pseudopath of an oriented read is the
    // sequence of MarkerIntervals it encounters.
    // This is indexed by OrientedReadId::getValue();
    MemoryMapped::VectorOfVectors<PseudoPathEntry, uint64_t> pseudoPaths;
    void computePseudoPaths();
    void computePseudoPathsPass1(size_t threadId);
    void computePseudoPathsPass2(size_t threadId);
    void computePseudoPathsPass12(uint64_t pass);
    void sortPseudoPaths(size_t threadId);
    class ComputePseudoPathsData {
    public:
    };
    ComputePseudoPathsData computePseudoPathsData;

    void writePseudoPaths(const string& fileName) const;
    void writePseudoPaths(ostream&) const;



    // Use transitions in pseudopaths to create edges (gfa links).
    void createEdges();


    // The virtual marker graph edges.
    // These don't exists in the marker graph and we
    // define them here.
    vector<VirtualMarkerGraphEdge> virtualMarkerGraphEdges;
};



// The AssemblyGraph is used to store the Mode 3 assembly graph,
// when it no longer needs to be changed,
// in memory mapped data structures.
class shasta::mode3::AssemblyGraph {
public:
};

#endif

