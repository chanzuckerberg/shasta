// Shasta.
#include "Assembler.hpp"
#include "AssemblyPathGraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>


void Assembler::detangle()
{
    // Check that we have what we need.
    SHASTA_ASSERT(markerGraph.edgeMarkerIntervals.isOpen());
    SHASTA_ASSERT(assemblyGraph.vertices.isOpen);
    SHASTA_ASSERT(assemblyGraph.reverseComplementVertex.isOpen);
    SHASTA_ASSERT(assemblyGraph.edges.isOpen);
    SHASTA_ASSERT(assemblyGraph.reverseComplementEdge.isOpen);
    SHASTA_ASSERT(assemblyGraph.edgesBySource.isOpen());
    SHASTA_ASSERT(assemblyGraph.edgesByTarget.isOpen());
    SHASTA_ASSERT(assemblyGraph.edgeLists.isOpen());
    SHASTA_ASSERT(assemblyGraph.reverseComplementEdge.isOpen);

    // Create the AssemblyGraphPath.
    // Initially, it is a faithful copy of the assembly graph.
    AssemblyPathGraph graph(assemblyGraph);



    // Fill in the oriented read ids of the edges.
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph) {
        AssemblyPathGraphEdge& edge = graph[e];

        // At this stage the path must be a single assembly graph edge.
        SHASTA_ASSERT(edge.path.size() == 1);
        const AssemblyGraph::EdgeId edgeId = edge.path.front();

        // Get the marker graph edges corresponding to this assembly graph edge.
        const auto markerGraphEdgeIds = assemblyGraph.edgeLists[edgeId];

        // Loop over these marker graph edges.
        std::set<OrientedReadId> orientedReadIds;
        for(const MarkerGraph::EdgeId markerGraphEdgeId: markerGraphEdgeIds) {

            // Loop over the marker intervals of this marker graph edge.
            const auto markerIntervals = markerGraph.edgeMarkerIntervals[markerGraphEdgeId];
            for(const MarkerInterval& markerInterval: markerIntervals) {
                orientedReadIds.insert(markerInterval.orientedReadId);
            }
        }

        // Copy the oriented read ids to the AssemblyPathGraphEdge.
        copy(orientedReadIds.begin(), orientedReadIds.end(),
            back_inserter(edge.orientedReadIds));

        // Also store the path length, measured on the marker graph.
        edge.pathLength = markerGraphEdgeIds.size();
    }
}

