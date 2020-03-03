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


    // Create the tangles.
    graph.createTangles();

    // Do the detangling.
    const double basesPerMarker =
        double(assemblerInfo->baseCount) /
        double(markers.totalSize()/2);
    graph.detangle(basesPerMarker);



    // Use the detangled AssemblyPathGraph to create a new AssemblyGraph.
    AssemblyGraph newAssemblyGraph;



    // Create the vertices of the new AssemblyGraph.
    vector< pair<MarkerGraph::VertexId, AssemblyPathGraph::vertex_descriptor> > newVertices;
    BGL_FORALL_VERTICES(v, graph, AssemblyPathGraph) {
        const AssemblyGraph::VertexId oldVertexId = graph[v].vertexId;
        const MarkerGraph::VertexId markerGraphVertexId = assemblyGraph.vertices[oldVertexId];
        newVertices.push_back(make_pair(markerGraphVertexId, v));
    }
    sort(newVertices.begin(), newVertices.end());

    // The new vertices are sorted by marker graph vertex id.
    // The position in the newVertices vector is the vertex id in the new assembly graph.
    newAssemblyGraph.vertices.createNew(
        largeDataName("New-AssemblyGraphVertices"),
        largeDataPageSize);
    newAssemblyGraph.vertices.resize(newVertices.size());
    for(AssemblyGraph::VertexId vertexId=0; vertexId<newVertices.size(); vertexId++) {
        newAssemblyGraph.vertices[vertexId] = newVertices[vertexId].first;
    }
    cout << "The detangled assembly graph has " <<
        newVertices.size() << " vertices." << endl;



    // Find the reverse complement of each vertex.
    newAssemblyGraph.reverseComplementVertex.createNew(
        largeDataName("New-AssemblyGraphReverseComplementVertex"), largeDataPageSize);
    newAssemblyGraph.reverseComplementVertex.resize(newAssemblyGraph.vertices.size());
    for(AssemblyGraph::VertexId vertexId=0; vertexId<newAssemblyGraph.vertices.size(); vertexId++) {
        const AssemblyPathGraph::vertex_descriptor v = newVertices[vertexId].second;
        const AssemblyPathGraphVertex& vertex = graph[v];

        // Use the Rc suffix for reverse complement ids.
        const AssemblyPathGraph::vertex_descriptor vRc = vertex.reverseComplementVertex;
        const AssemblyPathGraphVertex& vertexRc = graph[vRc];
        const AssemblyGraph::VertexId oldAssemblyGraphVertexIdRc = vertexRc.vertexId;
        const MarkerGraph::VertexId markerGraphVertexIdRc = assemblyGraph.vertices[oldAssemblyGraphVertexIdRc];

        // Look it up in our sorted newVertices vector.
        const auto it = std::lower_bound(newVertices.begin(), newVertices.end(),
            make_pair(markerGraphVertexIdRc, vRc));
        SHASTA_ASSERT(it != newVertices.end());
        SHASTA_ASSERT(*it == make_pair(markerGraphVertexIdRc, vRc));

        newAssemblyGraph.reverseComplementVertex[vertexId] = it - newVertices.begin();
    }

    // Sanity check on the reverse complement vertices.
    for(AssemblyGraph::VertexId vertexId=0; vertexId<newAssemblyGraph.vertices.size(); vertexId++) {
        const AssemblyGraph::VertexId vertexIdRc = newAssemblyGraph.reverseComplementVertex[vertexId];
        SHASTA_ASSERT(newAssemblyGraph.reverseComplementVertex[vertexIdRc] == vertexId);
    }



    // Create edges of the new assembly graph.
    vector<AssemblyPathGraph::edge_descriptor> newEdges;
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph) {
        newEdges.push_back(e);
    }
    cout << "The detangled assembly graph has " <<
        newEdges.size() << " edges." << endl;

    newAssemblyGraph.edges.createNew(
        largeDataName("New-AssemblyGraphEdges"),
        largeDataPageSize);
    newAssemblyGraph.edgeLists.createNew(
        largeDataName("New-AssemblyGraphEdgeLists"),
        largeDataPageSize);
    for(AssemblyGraph::EdgeId newEdgeId=0; newEdgeId<newEdges.size(); newEdgeId++) {
        const AssemblyPathGraph::edge_descriptor e = newEdges[newEdgeId];
        const AssemblyPathGraphEdge edge = graph[e];
        const AssemblyPathGraph::vertex_descriptor v0 = source(e, graph);
        const AssemblyPathGraph::vertex_descriptor v1 = target(e, graph);

        // Find the corresponding vertex ids in the new assembly graph.
        const auto p0 = make_pair(assemblyGraph.vertices[graph[v0].vertexId], v0);
        const auto it0 = std::lower_bound(newVertices.begin(), newVertices.end(), p0);
        SHASTA_ASSERT(it0 != newVertices.end());
        SHASTA_ASSERT(*it0 == p0);
        const AssemblyGraph::VertexId newVertexId0 = it0 - newVertices.begin();

        const auto p1 = make_pair(assemblyGraph.vertices[graph[v1].vertexId], v1);
        const auto it1 = std::lower_bound(newVertices.begin(), newVertices.end(), p1);
        SHASTA_ASSERT(it1 != newVertices.end());
        SHASTA_ASSERT(*it1 == p1);
        const AssemblyGraph::VertexId newVertexId1 = it1 - newVertices.begin();

        // Create ands store the new AssemblyGraph::Edge.
        AssemblyGraph::Edge newEdge;
        newEdge.source = newVertexId0;
        newEdge.target = newVertexId1;
        newAssemblyGraph.edges.push_back(newEdge);

        // Now store the marker graph path corresponding to this edge.
        newAssemblyGraph.edgeLists.appendVector();
        for(const AssemblyGraph::EdgeId oldAssemblyGraphEdgeId: edge.path) {
            const span<MarkerGraph::EdgeId> partialPath = assemblyGraph.edgeLists[oldAssemblyGraphEdgeId];
            for(const MarkerGraph::EdgeId markerGraphEdgeId: partialPath) {
                newAssemblyGraph.edgeLists.append(markerGraphEdgeId);
            }
        }

    }
    SHASTA_ASSERT(newAssemblyGraph.edges.size() == newAssemblyGraph.edgeLists.size());



    // Compute connectivity of the new assembly graph.
    newAssemblyGraph.edgesBySource.createNew(
        largeDataName("New-AssemblyGraphEdgesBySource"),
        largeDataPageSize);
    newAssemblyGraph.edgesByTarget.createNew(
        largeDataName("New-AssemblyGraphEdgesByTarget"),
        largeDataPageSize);
    newAssemblyGraph.computeConnectivity();
}

