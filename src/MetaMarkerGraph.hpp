#ifndef SHASTA_META_MARKER_GRAPH_HPP
#define SHASTA_META_MARKER_GRAPH_HPP

// Experimental. See AssemblerAnalyzePaths.cpp for more information.

// Shasta.
#include "AssemblyGraph.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {
    class MetaMarkerGraph;
    class MetaMarkerGraphVertex;
    class MetaMarkerGraphEdge;

    using MetaMarkerGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        MetaMarkerGraphVertex,
        MetaMarkerGraphEdge
        >;

}



class shasta::MetaMarkerGraphVertex {
public:

    uint64_t vertexId;

    // The id of the assembly graph edge (segment)
    // corresponding to this vertex of the MetaMarkerGraph.
    AssemblyGraph::EdgeId segmentId;

    // The number of markers on that segment,
    uint64_t markerCount;

    // The OrientedReadIds and meta-ordinals.
    // The meta-ordinal is the position of this segment
    // in the oriented read pseudo-path.
    vector< pair<OrientedReadId, uint64_t> > orientedReads;

    MetaMarkerGraphVertex(
        uint64_t vertexId,
        AssemblyGraph::EdgeId segmentId,
        uint64_t markerCount,
        const vector< pair<OrientedReadId, uint64_t> >& orientedReads) :
        vertexId(vertexId),
        segmentId(segmentId),
        markerCount(markerCount),
        orientedReads(orientedReads)
    {
    }
};



class shasta::MetaMarkerGraphEdge {
public:

    // The OrientedReadIds and meta-ordinals.
    // The meta-ordinal is the position of the source vertex of this edge
    // in the oriented read pseudo-path.
    vector< pair<OrientedReadId, uint64_t> > orientedReads;
};



class shasta::MetaMarkerGraph : public MetaMarkerGraphBaseClass {
public:
    void createEdges();
    void writeGraphviz(const string& fileName) const;
    void writeGfa(const string& fileName) const;
    void writeVerticesCsv(const string& fileName) const;
    void writeEdgesCsv(const string& fileName) const;
};

#endif
