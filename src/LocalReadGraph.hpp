#ifndef SHASTA_LOCAL_READ_GRAPH_HPP
#define SHASTA_LOCAL_READ_GRAPH_HPP



/*******************************************************************************

The local read graph is a local subgraph of the global read graph,
starting at a given read and going out up to a specified distance.
It is an undirected graph in which each vertex corresponds to a read.

See comments at the beginning of AssemblerReadGraph.cpp for
more information on the global read graph.

*******************************************************************************/

// Shasta.
#include <computeLayout.hpp>
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard libraries.
#include "array.hpp"
#include <map>

namespace shasta {

    // Forward declaration of types declared in this file.
    class LocalReadGraphVertex;
    class LocalReadGraphEdge;
    class LocalReadGraph;
    using LocalReadGraphBaseClass = boost::adjacency_list<
        boost::setS,
        boost::listS,
        boost::undirectedS,
        LocalReadGraphVertex,
        LocalReadGraphEdge
        >;

    enum class AlignmentType;
}


class shasta::LocalReadGraphVertex {
public:

    OrientedReadId orientedReadId;
    uint32_t orientedReadIdValue;

    // The number of markers in this read.
    uint32_t markerCount;

    // Flag that indicates whether the associated read is chimeric.
    bool isChimeric;

    // The distance of this vertex from the starting vertex.
    uint32_t distance;

    // Used for Blast annotations.
    string additionalToolTipText;

    // The least square position of this vertex,
    // computed by Assembler::analyzeLocalReadGraph.
    double leastSquarePosition;

    LocalReadGraphVertex(
        OrientedReadId orientedReadId,
        uint32_t markerCount,
        bool isChimeric,
        uint32_t distance) :
        orientedReadId(orientedReadId),
        orientedReadIdValue(orientedReadId.getValue()),
        markerCount(markerCount),
        isChimeric(isChimeric),
        distance(distance)
        {}

    // The position of this vertex in the layout for display,
    //  as computed by computeSfdpLayout.
    array<double, 2> position;

};



class shasta::LocalReadGraphEdge {
public:

    // The number of alignment markers in the alignment
    // that created this edge.
    uint32_t markerCount;

    // Flag that indicates this edge jumps across strands.
    bool crossesStrands;

    // The corresponding global read graph edge.
    uint64_t globalEdgeId;

    // The average alignment offset.
    double averageAlignmentOffset = 0.;

    string color;

    LocalReadGraphEdge(
        uint32_t markerCount,
        bool crossesStrands,
        uint64_t globalEdgeId) :
        markerCount(markerCount),
        crossesStrands(crossesStrands),
        globalEdgeId(globalEdgeId)
        {}
};



class shasta::LocalReadGraph :
    public LocalReadGraphBaseClass {
public:

    void addVertex(
        OrientedReadId,
        uint32_t baseCount,
        bool isChimeric,
        uint32_t distance);

    void addEdge(
        OrientedReadId,
        OrientedReadId,
        uint32_t markerCount,
        uint64_t globalEdgeId,
        bool crossesStrands);

    // Find out if a vertex with a given OrientedReadId exists.
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
    void write(
            ostream&,
            const string& layoutMethod,
            uint32_t maxDistance,
            double vertexScalingFactor,
            double edgeThicknessScalingFactor) const;
    void write(
            const string& fileName,
            const string& layoutMethod,
            uint32_t maxDistance,
            double vertexScalingFactor,
            double edgeThicknessScalingFactor) const;

private:

    // Map that gives the vertex corresponding to an OrientedReadId.
    std::map<OrientedReadId, vertex_descriptor> vertexMap;

    // Graphviz writer.
    class Writer {
    public:
        Writer(
                const LocalReadGraph&,
                const string& layoutMethod,
                uint32_t maxDistance,
                double vertexScalingFactor,
                double edgeThicknessScalingFactor);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const LocalReadGraph& graph;
        const string& layoutMethod;
        uint32_t maxDistance;
        double vertexScalingFactor;
        double edgeThicknessScalingFactor;
        double edgeArrowScalingFactor;
    };
};



#endif
