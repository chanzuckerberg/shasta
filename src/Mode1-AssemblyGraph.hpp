#ifndef SHASTA_MODE1_ASSEMBLY_GRAPH_HPP
#define SHASTA_MODE1_ASSEMBLY_GRAPH_HPP

#include "MarkerGraph.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

#include <boost/graph/adjacency_list.hpp>

namespace shasta {
    namespace Mode1 {

        class AssemblyGraph;
        class AssemblyGraphVertex;
        class AssemblyGraphEdge;

        using AssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            AssemblyGraphVertex,
            AssemblyGraphEdge
            >;


    }

    class CompressedMarker;
    class MarkerGraph;
}


/*******************************************************************************

Class Mode1::AssemblyGraph is the assembly graph class used in mode 1
assembly (--Assembly.mode 1).

Here, and in contrast to shasta::AssemblyGraph, each vertex corresponds
to a path in the marker graph and, eventually, to an assembled segment.

For development flexibility, this class is currently implemented
using the Boost Graph library. and
without multithreading and without usage of MemoryMapped containers.
This can result in performance degradation.

*******************************************************************************/



class shasta::Mode1::AssemblyGraphVertex {
public:

};



class shasta::Mode1::AssemblyGraphEdge {
public:

};



class shasta::Mode1::AssemblyGraph : public AssemblyGraphBaseClass {
public:

    AssemblyGraph(
        uint64_t minEdgeCoverage,
        uint64_t minEdgeCoveragePerStrand,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&);

private:

    void createVertices(
        uint64_t minEdgeCoverage,
        uint64_t minEdgeCoveragePerStrand);

    // Constructor parameters.
    uint64_t minEdgeCoverage;
    uint64_t minEdgeCoveragePerStrand;

    // References to data owned by the Assembler.
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MarkerGraph& markerGraph;

    // Return true if a given marker graph edge has sufficient coverage
    // (total and for each strand).
    bool markerGraphEdgeHasSufficientCoverage(MarkerGraph::EdgeId) const;

    // Return the out-degree of a marker graph vertex,
    // counting only edges with sufficient coverage.
    uint64_t markerGraphVertexOutDegree(MarkerGraph::VertexId) const;

    // Return the in-degree of a marker graph vertex,
    // counting only edges with sufficient coverage.
    uint64_t markerGraphVertexInDegree(MarkerGraph::VertexId) const;

    // Return the unique next/previous marker graph edge for a given marker graph edge,
    // or MarkerGraph::invalidEdgeId if there are none or more than one.
    // The next/previous marker graph edge is chosen among the
    // ones with sufficient coverage.
    MarkerGraph::EdgeId getMarkerGraphUniqueNextEdge(MarkerGraph::EdgeId) const;
    MarkerGraph::EdgeId getMarkerGraphUniquePreviousEdge(MarkerGraph::EdgeId) const;
};

#endif
