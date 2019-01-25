#ifndef CZI_SHASTA_LOCAL_READ_GRAPH_HPP
#define CZI_SHASTA_LOCAL_READ_GRAPH_HPP



/*******************************************************************************

The local read graph is a local subgraph of the global read graph,
starting at a given read and going out up to a specified distance.
It is a bidirected graph in which each vertex corresponds to a read.

See comments at the beginning of AssemblerReadGraph.cpp for
more information on the global read graph.

*******************************************************************************/

// Shasta.
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard libraries.
#include <map>

namespace ChanZuckerberg {
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

    }
}


class ChanZuckerberg::shasta::LocalReadGraphVertex {
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

};



class ChanZuckerberg::shasta::LocalReadGraphEdge {
public:

    // The id of the global read graph edge that corresponds to this edge.
    // This is an index into ReadGraph::edges.
    size_t globalEdgeId;

    LocalReadGraphEdge(size_t globalEdgeId) :
        globalEdgeId(globalEdgeId) {}
};



class ChanZuckerberg::shasta::LocalReadGraph :
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
        size_t globalEdgeId);

    // Find out if a vertex with a given OrientedReadId exists.
    bool vertexExists(OrientedReadId) const;

    // Get the distance of an existing vertex from the start vertex.
    uint32_t getDistance(OrientedReadId) const;

    // Write in Graphviz format.
    void write(ostream&, uint32_t maxDistance) const;
    void write(const string& fileName, uint32_t maxDistance) const;

private:

    // Map that gives the vertex corresponding to an OrientedReadId.
    std::map<OrientedReadId, vertex_descriptor> vertexMap;

    // Graphviz writer.
    class Writer {
    public:
        Writer(const LocalReadGraph&, uint32_t maxDistance);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const LocalReadGraph& graph;
        uint32_t maxDistance;
    };
};



#endif
