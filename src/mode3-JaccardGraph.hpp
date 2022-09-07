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

#include "mode3-SegmentPairInformation.hpp"
#include "MemoryMappedVector.hpp"
#include "cstdint.hpp"

namespace shasta {
    namespace mode3 {
        class JaccardGraph;
        class JaccardGraphEdge;
    }
}



class shasta::mode3::JaccardGraphEdge {
public:

    // The source segment.
    uint64_t segmentId0;

    // The target segment.
    uint64_t segmentId1;

    // The SegmentPairInformation computed by
    // mode3::AssemblyGraph::analyzeSegmentPair
    // when called for (segmentId0, segmentId1), in this order.
    SegmentPairInformation segmentPairInformation;
};



class shasta::mode3::JaccardGraph {
public:
    MemoryMapped::Vector<JaccardGraphEdge> edges;
};



#endif
