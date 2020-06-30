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

    // The id of the assembly graph edge (segment)
    // corresponding to this vertex of the MetaMarkerGraph.
    AssemblyGraph::EdgeId segmentId;

    // The sequence number among vertices corresponding to the same segmentId.
    uint64_t sequenceNumber = std::numeric_limits<uint64_t>::max();

    // The number of markers on that segment.
    uint64_t markerCount;

    // The OrientedReadIds and meta-ordinals.
    // The meta-ordinal is the position of this segment
    // in the oriented read pseudo-path.
    vector< pair<OrientedReadId, uint64_t> > orientedReads;

    MetaMarkerGraphVertex(
        AssemblyGraph::EdgeId segmentId,
        uint64_t markerCount,
        const vector< pair<OrientedReadId, uint64_t> >& orientedReads) :
        segmentId(segmentId),
        markerCount(markerCount),
        orientedReads(orientedReads)
    {
    }
    MetaMarkerGraphVertex() {}

    string gfaId() const
    {
        return to_string(segmentId) + "_" + to_string(sequenceNumber);
    }

    uint64_t coverage() const
    {
        return orientedReads.size();
    }

    // Color the vertex by coverage.
    string color() const
    {
        const auto coverage = orientedReads.size();
        SHASTA_ASSERT(coverage > 0);
        if(coverage == 1) {
            return "red";
        } else if(coverage == 2) {
            return "orange";
        } else if(coverage == 3) {
            return "green";
        } else {
            return "black";
        }
    }
};



class shasta::MetaMarkerGraphEdge {
public:

    // The OrientedReadIds and meta-ordinals.
    // The meta-ordinal is the position of the source vertex of this edge
    // in the oriented read pseudo-path.
    vector< pair<OrientedReadId, uint64_t> > orientedReads;

    // Color the edge by coverage.
    string color() const
    {
        const auto coverage = orientedReads.size();
        SHASTA_ASSERT(coverage > 0);
        if(coverage == 1) {
            return "red";
        } else if(coverage == 2) {
            return "orange";
        } else if(coverage == 3) {
            return "green";
        } else {
            return "black";
        }
    }
};



class shasta::MetaMarkerGraph : public MetaMarkerGraphBaseClass {
public:
    using SegmentId = AssemblyGraph::EdgeId;

    void addVertex(
        SegmentId segmentId,
        uint64_t markerCount,
        const vector< pair<OrientedReadId, uint64_t> >& orientedReads);

    void createEdges();
    void transitiveReduction();
    void generateSequenceNumbers();
    void writeGraphviz(const string& fileName,
        SegmentId startSegmentId = std::numeric_limits<SegmentId>::max()) const;
    void writeGfa(const string& fileName) const;
    void writeVerticesCsv(const string& fileName) const;
    void writeVerticesDetailCsv(const string& fileName) const;
    void writeEdgesCsv(const string& fileName) const;

    void findLinearChain(
        SegmentId segmentId,
        vector<SegmentId>& chain) const;

    void findForwardChokePoints(
        SegmentId segmentId,
        vector<vertex_descriptor>&);
    void findBackwardChokePoints(
        SegmentId segmentId,
        vector<vertex_descriptor>&);

    // Recursively prune all leafs with coverage less than minCoverage.
    void prune(uint64_t minCoverage);

private:

    // Find the vertex corresponding to a given segment id,
    // or null_vertex if not exactly one found.
    vertex_descriptor findVertex(SegmentId) const;

};

#endif
