#ifndef SHASTA_LOCAL_ALIGNMENT_CANDIDATE_GRAPH_HPP
#define SHASTA_LOCAL_ALIGNMENT_CANDIDATE_GRAPH_HPP



/*******************************************************************************

The local candidate graph is a local subgraph of the global candidate pairs,
starting at a given oriented read and going out up to a specified distance.
It is an undirected graph in which each vertex corresponds to a read.
An edge is created between two vertices if the corresponding
oriented reads have a stored alignment.

It has been augmented from the LocalAlignmentGraph to include additional
attributes that describe whether each edge exists in the AlignmentGraph,
the ReadGraph, or (optionally) whether an overlap is inferred to exist in a
reference alignment of reads.
*******************************************************************************/

// Shasta.
#include "Alignment.hpp"
#include "computeLayout.hpp"
#include "OrientedReadPair.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard libraries.
#include <map>

namespace shasta {

    class LocalAlignmentCandidateGraphVertex;
    class LocalAlignmentCandidateGraphEdge;
    class LocalAlignmentCandidateGraph;
    using LocalAlignmentCandidateGraphBaseClass = boost::adjacency_list<
        boost::setS,
        boost::listS,
        boost::undirectedS,
        LocalAlignmentCandidateGraphVertex,
        LocalAlignmentCandidateGraphEdge
        >;

}


class shasta::LocalAlignmentCandidateGraphVertex {
public:

    // The OrientedReadId that this vertex corresponds to.
    // We store it as OrientedRead::Int so we can also use
    // it as a vertex id for graphviz output.
    OrientedReadId::Int orientedReadId;

    // The number of bases in this read.
    uint32_t baseCount;

    // Layout position of this vertex, for display.
    array<double, 2> position;

    // The distance of this vertex from the starting vertex.
    uint32_t distance;

    uint8_t getSvgOrdering() const;

    LocalAlignmentCandidateGraphVertex(
        OrientedReadId orientedReadId,
        uint32_t baseCount,
        uint32_t distance) :
        orientedReadId(orientedReadId.getValue()),
        baseCount(baseCount),
        distance(distance)
        {}

};



class shasta::LocalAlignmentCandidateGraphEdge {
public:

    bool inAlignments;
    bool inReadGraph;
    bool inReferenceAlignments;

    uint8_t getSvgOrdering() const;

    string getSvgClassName() const;

    LocalAlignmentCandidateGraphEdge(
        bool inAlignments,
        bool inReadgraph,
        bool inReferenceAlignments) :
            inAlignments(inAlignments),
            inReadGraph(inReadgraph),
            inReferenceAlignments(inReferenceAlignments)
    {};

    LocalAlignmentCandidateGraphEdge():
            inAlignments(false),
            inReadGraph(false),
            inReferenceAlignments(false)
    {};
};



class shasta::LocalAlignmentCandidateGraph :
    public LocalAlignmentCandidateGraphBaseClass {
public:

    void addVertex(
        OrientedReadId orientedReadId,
        uint32_t baseCount,
        uint32_t distance);

    void addEdge(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        bool inAlignments,
        bool inReadgraph,
        bool inReferenceAlignments);

    // Find out if a vertex with a given OrientedId exists.
    bool vertexExists(OrientedReadId) const;

    // Get the distance of an existing vertex from the start vertex.
    uint32_t getDistance(OrientedReadId) const;

    // Compute sfdp layout using graphviz and store the results
    // in the vertex positions.
    ComputeLayoutReturnCode computeLayout(
        const string& layoutMethod,
        double timeout);

    // Write directly to svg, without using Graphviz rendering.
    // This assumes that the layout was already computed
    // and stored in the vertices.
    void writeSvg(
        const string& svgId,
        uint64_t width,
        uint64_t height,
        double vertexScalingFactor,
        double edgeThicknessScalingFactor,
        uint64_t maxDistance,
        ostream& svg) const;

    // Write in Graphviz format.
    void write(ostream&, uint32_t maxDistance) const;
    void write(const string& fileName, uint32_t maxDistance) const;

    uint8_t getEdgeOrdering(const LocalAlignmentCandidateGraphEdge& e);
    uint8_t getVertexOrdering(const LocalAlignmentCandidateGraphVertex& v);

private:

    // Map that gives the vertex corresponding to an OrientedReadId.
    std::map<OrientedReadId, vertex_descriptor> vertexMap;

    // Graphviz writer.
    class Writer {
    public:
        Writer(const LocalAlignmentCandidateGraph&, uint32_t maxDistance);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const LocalAlignmentCandidateGraph& graph;
        uint32_t maxDistance;
    };
};



#endif
