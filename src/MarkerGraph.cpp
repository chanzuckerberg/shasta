#include "MarkerGraph.hpp"
using namespace ::shasta;
using namespace ChanZuckerberg::shasta;

const MarkerGraph::VertexId MarkerGraph::invalidVertexId = std::numeric_limits<VertexId>::max();
const MarkerGraph::EdgeId MarkerGraph::invalidEdgeId = std::numeric_limits<EdgeId>::max();
const MarkerGraph::CompressedVertexId
	MarkerGraph::invalidCompressedVertexId = std::numeric_limits<uint64_t>::max();



// Locate the edge given the vertices.
const MarkerGraph::Edge*
    MarkerGraph::findEdge(Uint40 source, Uint40 target) const
{
    const auto edgesWithThisSource = edgesBySource[source];
    for(const uint64_t i: edgesWithThisSource) {
        const Edge& edge = edges[i];
        if(edge.target == target) {
            return &edge;
        }
    }
    return 0;
}


MarkerGraph::EdgeId MarkerGraph::findEdgeId(Uint40 source, Uint40 target) const
{
	const Edge* edgePointer = findEdge(source, target);
	SHASTA_ASSERT(edgePointer);
	return edgePointer - edges.begin();
}
