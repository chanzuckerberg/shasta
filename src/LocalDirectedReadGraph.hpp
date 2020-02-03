#ifndef SHASTA_LOCAL_DIRECTED_READ_GRAPH_HPP
#define SHASTA_LOCAL_DIRECTED_READ_GRAPH_HPP




// Shasta.
#include "Alignment.hpp"
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
    uint32_t baseCount;
    uint32_t markerCount;

    // The distance of this vertex from the starting vertex.
    uint32_t distance;

    bool isContained;

    // Used for Blast annotations.
    string additionalToolTipText;

    LocalDirectedReadGraphVertex(
        OrientedReadId orientedReadId,
        uint32_t baseCount,
        uint32_t markerCount,
        uint32_t distance,
        bool isContained) :
        orientedReadId(orientedReadId),
        orientedReadIdValue(orientedReadId.getValue()),
        baseCount(baseCount),
        markerCount(markerCount),
        distance(distance),
        isContained(isContained)
        {}

    // Information from the global conflict read graph.
    bool isConflict = false;
    uint64_t clusterId = std::numeric_limits<uint64_t>::max();
    uint64_t conflictCount = 0;
    bool hasLongGap = false;

    bool isConflictingGreen = false;
    bool isConflictingRed = false;

};



class shasta::LocalDirectedReadGraphEdge {
public:

    AlignmentInfo alignmentInfo;

    bool involvesTwoContainedVertices;
    bool involvesOneContainedVertex;
    bool keep;
    bool isConflict;

    uint32_t commonNeighborCount;

    LocalDirectedReadGraphEdge(
        const AlignmentInfo& alignmentInfo,
        bool involvesTwoContainedVertices,
        bool involvesOneContainedVertex,
        bool keep,
        bool isConflict,
        uint32_t commonNeighborCount):
        alignmentInfo(alignmentInfo),
        involvesTwoContainedVertices(involvesTwoContainedVertices),
        involvesOneContainedVertex(involvesOneContainedVertex),
        keep(keep),
        isConflict(isConflict),
        commonNeighborCount(commonNeighborCount)
        {}
};



class shasta::LocalDirectedReadGraph :
    public LocalDirectedReadGraphBaseClass {
public:

    void addVertex(
        OrientedReadId,
        uint32_t baseCount,
        uint32_t markerCount,
        uint32_t distance,
        bool isContained);

    void addEdge(
        OrientedReadId,
        OrientedReadId,
        const AlignmentInfo& alignmentInfo,
        bool involvesTwoContainedVertices,
        bool involvesOneContainedVertex,
        bool keep,
        bool isConflict,
        uint32_t commonNeighborCount);

    // Find out if a vertex with a given OrientedReadId exists.
    bool vertexExists(OrientedReadId) const;

    // Get the distance of an existing vertex from the start vertex.
    uint32_t getDistance(OrientedReadId) const;

    // Write in Graphviz format.
    enum class VertexColoringMethod {
        None, ByConflictCount, ByCluster
    };
    void write(
        ostream&,
        uint32_t maxDistance,
        double vertexScalingFactor,
        double edgeThicknessScalingFactor,
        double edgeArrowScalingFactor,
        bool colorEdgeArrows,
        VertexColoringMethod) const;
    void write(
        const string& fileName,
        uint32_t maxDistance,
        double vertexScalingFactor,
        double edgeThicknessScalingFactor,
        double edgeArrowScalingFactor,
        bool colorEdgeArrows,
        VertexColoringMethod) const;

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
            uint32_t maxDistance,
            double vertexScalingFactor,
            double edgeThicknessScalingFactor,
            double edgeArrowScalingFactor,
            bool colorEdgeArrows,
            VertexColoringMethod);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const LocalDirectedReadGraph& graph;
        uint32_t maxDistance;
        double vertexScalingFactor;
        double edgeThicknessScalingFactor;
        double edgeArrowScalingFactor;
        bool colorEdgeArrows;
        VertexColoringMethod vertexColoringMethod;
    };
};



#endif
