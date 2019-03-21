#include "ReadGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



// Create a vertex for each of the two oriented reads
// corresponding to readCount reads.
// Used to create a DynamicReadGraph representing
// the entire global read graph.
DynamicReadGraph::DynamicReadGraph(ReadId readCount)
{
    // for()
}


// Count the triangles that have an edge as one of the sides.
size_t ReadGraph::countTriangles(uint32_t edgeId01) const
{
    const Edge& edge01 = edges[edgeId01];

    // Find neighbors of each of the vertices of this edge.
    array<vector<OrientedReadId>, 2> neighbors;
    for(size_t i=0; i<2; i++) {
        const OrientedReadId orientedReadId = edge01.orientedReadIds[i];
        const MemoryAsContainer<const uint32_t>edgeIds = connectivity[orientedReadId.getValue()];
        for(const uint32_t edgeId: edgeIds) {
            if(edgeId == edgeId01) {
                continue;
            }
            neighbors[i].push_back(edges[edgeId].getOther(orientedReadId));
        }
        sort(neighbors[i].begin(), neighbors[i].end());
    }

    // Find common neighbors.
    vector<OrientedReadId> commonNeighbors;
    std::set_intersection(
        neighbors[0].begin(), neighbors[0].end(),
        neighbors[1].begin(), neighbors[1].end(),
        back_inserter(commonNeighbors));
    return commonNeighbors.size();
}
