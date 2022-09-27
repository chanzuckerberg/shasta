#ifndef SHASTA_MODE3_JACCARD_GRAPH_HPP
#define SHASTA_MODE3_JACCARD_GRAPH_HPP

/*******************************************************************************

The mode3::JaccardGraph is a directed graph in which each vertex represents
a segment in the mode3::AssemblyGraph.

A directed edge S0->S1 is created if S0 and S1 have:
- A sufficient number of common reads.
- High Jaccard similarity.
- Low unexplained fractions.
(The above quantities defined as computed by
mode3::AssemblyGraph::analyzeSegmentPair).
For the edge to be created, we also require one of the following:
1. S1 is the first primary segment encountered starting from S0,
   and performing a forward path search using the algorithm defined by
   mode3::AssemblyGraph::createAssemblyPath3.
2. S0 is the first primary segment encountered starting from S1,
   and performing a backward path search using the algorithm defined by
   mode3::AssemblyGraph::createAssemblyPath3.

*******************************************************************************/

// Shasta.
#include "mode3-SegmentPairInformation.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "cstdint.hpp"
#include <map>
#include "tuple.hpp"
#include "vector.hpp"



namespace shasta {
    namespace mode3 {
        class JaccardGraph;
        class JaccardGraphEdge;
        class JaccardGraphEdgeInfo;
        class JaccardGraphVertex;

        using JaccardGraphBaseClass = boost::adjacency_list<
            boost::listS, boost::listS, boost::bidirectionalS,
            JaccardGraphVertex, JaccardGraphEdge>;

    }
}



class shasta::mode3::JaccardGraphVertex {
public:

    // The assembly graph segment corresponding to this vertex.
    uint64_t segmentId;
};



class shasta::mode3::JaccardGraphEdge {
public:

    // The SegmentPairInformation computed by
    // mode3::AssemblyGraph::analyzeSegmentPair
    // when called for (segmentId0, segmentId1), in this order.
    SegmentPairInformation segmentPairInformation;

    // Flags for the directions in which this edge was found
    // (0=forward, 1=backward).
    array<bool, 2> wasFoundInDirection = {false, false};
    bool wasFoundInBothDirections() const
    {
        return wasFoundInDirection[0] and wasFoundInDirection[1];
    }

    JaccardGraphEdge(
        const SegmentPairInformation& segmentPairInformation,
        uint64_t direction) :
        segmentPairInformation(segmentPairInformation)
        {
            wasFoundInDirection[direction] = true;
        }
};



// This is only used during parallel creation of the edges.
class shasta::mode3::JaccardGraphEdgeInfo {
public:
    uint64_t segmentId0;
    uint64_t segmentId1;
    uint64_t direction;
    SegmentPairInformation segmentPairInformation;
};



class shasta::mode3::JaccardGraph : public JaccardGraphBaseClass {
public:

    // Create a JaccardGraph with the given number of vertices
    // (one for each segment) and no edges.
    JaccardGraph(uint64_t segmentCount);

    // Map segment ids to vertices.
    // If  vertex is removed, the corresponding entry will be null_vertex().
    vector<vertex_descriptor> vertexTable;

    // The edges found by each thread.
    // Only used during edge creation.
    vector< vector<JaccardGraphEdgeInfo> > threadEdges;
};



#endif
