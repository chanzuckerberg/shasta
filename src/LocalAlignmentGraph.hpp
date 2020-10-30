#ifndef SHASTA_LOCAL_ALIGNMENT_GRAPH_HPP
#define SHASTA_LOCAL_ALIGNMENT_GRAPH_HPP



/*******************************************************************************

The local alignment graph is a local subgraph of the global alignment graph,
starting at a given oriented read and going out up to a specified distance.
It is an undirected graph in which each vertex corresponds to a read.
An edge is created between two vertices if the corresponding
oriented reads have a stored alignment.

The alignment graph is being phased out in favor of the read graph
(see AssemblerReadGraph.cpp).

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

    class LocalAlignmentGraphVertex;
    class LocalAlignmentGraphEdge;
    class LocalAlignmentGraph;
    using LocalAlignmentGraphBaseClass = boost::adjacency_list<
        boost::setS,
        boost::listS,
        boost::undirectedS,
        LocalAlignmentGraphVertex,
        LocalAlignmentGraphEdge
        >;

}


class shasta::LocalAlignmentGraphVertex {
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

    LocalAlignmentGraphVertex(
        OrientedReadId orientedReadId,
        uint32_t baseCount,
        uint32_t distance) :
        orientedReadId(orientedReadId.getValue()),
        baseCount(baseCount),
        distance(distance)
        {}

};



class shasta::LocalAlignmentGraphEdge {
public:

    // Copies of the AlignmentInfo that caused this edge to be created.
    AlignmentInfo alignmentInfo;

    LocalAlignmentGraphEdge(
        const AlignmentInfo& alignmentInfo) :
        alignmentInfo(alignmentInfo)
        {}
};



class shasta::LocalAlignmentGraph :
    public LocalAlignmentGraphBaseClass {
public:

    void addVertex(
        OrientedReadId orientedReadId,
        uint32_t baseCount,
        uint32_t distance);

    void addEdge(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        const AlignmentInfo&);

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

private:

    // Map that gives the vertex corresponding to an OrientedReadId.
    std::map<OrientedReadId, vertex_descriptor> vertexMap;

    // Graphviz writer.
    class Writer {
    public:
        Writer(const LocalAlignmentGraph&, uint32_t maxDistance);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const LocalAlignmentGraph& graph;
        uint32_t maxDistance;
    };
};



#endif
