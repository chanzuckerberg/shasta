#ifndef SHASTA_LOCAL_DIRECTED_READ_GRAPH_HPP
#define SHASTA_LOCAL_DIRECTED_READ_GRAPH_HPP




// Shasta.
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard libraries.
#include <map>

namespace shasta {

    // Forward declaration of types declared in this file.
    class LocalDirectedReadGraphVertex;
    class LocalDirectedReadGraphEdge;
    class LocalDirectedReadGraph;
    using LocalDirectedReadGraphBaseClass = boost::adjacency_list<
        boost::setS,
        boost::listS,
        boost::bidirectionalS,
        LocalDirectedReadGraphVertex,
        LocalDirectedReadGraphEdge
        >;
}


class shasta::LocalDirectedReadGraphVertex {
public:

    OrientedReadId orientedReadId;
    OrientedReadId::Int orientedReadIdValue;   // Used as vertex id.

    // The number of bases and markers in this read.
    uint64_t baseCount;
    uint64_t markerCount;

    // The distance of this vertex from the starting vertex.
    uint64_t distance;

    // Used for Blast annotations.
    string additionalToolTipText;

    LocalDirectedReadGraphVertex(
        OrientedReadId orientedReadId,
        uint64_t baseCount,
        uint64_t markerCount,
        uint64_t distance) :
        orientedReadId(orientedReadId),
        orientedReadIdValue(orientedReadId.getValue()),
        baseCount(baseCount),
        markerCount(markerCount),
        distance(distance)
        {}

};



class shasta::LocalDirectedReadGraphEdge {
public:

    // Twice the offset between centers.
    int twiceOffsetAtCenter;

    // The number of alignment markers in the alignment
    // that created this edge.
    uint64_t markerCount;

    bool wasRemovedByTransitiveReduction;

    LocalDirectedReadGraphEdge(
        int twiceOffsetAtCenter,
        uint64_t markerCount,
        bool wasRemovedByTransitiveReduction) :
        twiceOffsetAtCenter(twiceOffsetAtCenter),
        markerCount(markerCount),
        wasRemovedByTransitiveReduction(wasRemovedByTransitiveReduction)
        {}
};



class shasta::LocalDirectedReadGraph :
    public LocalDirectedReadGraphBaseClass {
public:

    void addVertex(
        OrientedReadId,
        uint64_t baseCount,
        uint64_t markerCount,
        uint64_t distance);

    void addEdge(
        OrientedReadId,
        OrientedReadId,
        int twiceOffsetAtCenter,
        uint64_t markerCount,
        bool wasRemovedByTransitiveReduction);

    // Find out if a vertex with a given OrientedReadId exists.
    bool vertexExists(OrientedReadId) const;

    // Get the distance of an existing vertex from the start vertex.
    uint64_t getDistance(OrientedReadId) const;

    // Write in Graphviz format.
    void write(
        ostream&,
        uint64_t maxDistance,
        double vertexScalingFactor,
        double edgeThicknessScalingFactor,
        double edgeArrowScalingFactor) const;
    void write(
        const string& fileName,
        uint64_t maxDistance,
        double vertexScalingFactor,
        double edgeThicknessScalingFactor,
        double edgeArrowScalingFactor) const;

    // Return the vertex corresponding to a given OrientedReadId,
    // or null_vertex() if none.
    vertex_descriptor getVertex(OrientedReadId orientedReadId) const
    {
        const auto it = vertexMap.find(orientedReadId);
        if(it == vertexMap.end()) {
            return null_vertex();
        } else {
            return it->second;
        }
    }

private:

    // Map that gives the vertex corresponding to an OrientedReadId.
    std::map<OrientedReadId, vertex_descriptor> vertexMap;

    // Graphviz writer.
    class Writer {
    public:
        Writer(
            const LocalDirectedReadGraph&,
            uint64_t maxDistance,
            double vertexScalingFactor,
            double edgeThicknessScalingFactor,
            double edgeArrowScalingFactor);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const LocalDirectedReadGraph& graph;
        uint64_t maxDistance;
        double vertexScalingFactor;
        double edgeThicknessScalingFactor;
        double edgeArrowScalingFactor;
    };
};



#endif
