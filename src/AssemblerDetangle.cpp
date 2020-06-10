// Shasta.
#include "Assembler.hpp"
#include "AssemblyPathGraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

#include <cstdlib>


void Assembler::detangle()
{
    AssemblyGraph& assemblyGraph = *assemblyGraphPointer;

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

    cout << timestamp << "Before detangling, the assembly graph has " <<
        assemblyGraph.vertices.size() << " vertices and " <<
        assemblyGraph.edges.size() << " edges." << endl;

    // Create the AssemblyGraphPath.
    // Initially, it is a faithful copy of the assembly graph.
    AssemblyPathGraph graph(assemblyGraph);



    // Fill in the oriented read ids of the edges.
    cout << timestamp << "Filling in oriented reads." << endl;
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
    cout << timestamp << "Creating the tangles." << endl;
    graph.createTangles();

    // Do the detangling.
    const double basesPerMarker =
        double(assemblerInfo->baseCount) /
        double(markers.totalSize()/2);
    graph.detangle(basesPerMarker, assemblyGraph);



    // Use the detangled AssemblyPathGraph to create a new AssemblyGraph.
    shared_ptr<AssemblyGraph> newAssemblyGraphPointer = make_shared<AssemblyGraph>();
    AssemblyGraph& newAssemblyGraph = *newAssemblyGraphPointer;



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
    std::map<AssemblyPathGraph::edge_descriptor, AssemblyGraph::EdgeId> newEdgesMap;
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph) {
        newEdgesMap.insert(make_pair(e, newEdges.size()));
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
    ofstream csv("DetangleMap.csv");
    csv << "Path before detangle,Edge after detangle\n";
    for(AssemblyGraph::EdgeId newEdgeId=0; newEdgeId<newEdges.size(); newEdgeId++) {
        const AssemblyPathGraph::edge_descriptor e = newEdges[newEdgeId];
        const AssemblyPathGraphEdge edge = graph[e];
        csv << edge << "," << newEdgeId << "\n";
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


    // Find reverse complement edges of the new assembly graph.
    newAssemblyGraph.reverseComplementEdge.createNew(
        largeDataName("New-AssemblyGraphReverseComplementEdge"), largeDataPageSize);
    newAssemblyGraph.reverseComplementEdge.resize(newAssemblyGraph.edges.size());
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<newAssemblyGraph.edges.size(); edgeId++) {
        const AssemblyPathGraph::edge_descriptor e = newEdges[edgeId];
        const AssemblyPathGraph::edge_descriptor eRc = graph[e].reverseComplementEdge;
        newAssemblyGraph.reverseComplementEdge[edgeId] = newEdgesMap[eRc];
    }



    // Fill in coverage metrics for the detangled assembly graph.
    for(AssemblyGraph::EdgeId assemblyGraphEdgeId=0;
        assemblyGraphEdgeId<newAssemblyGraph.edges.size(); assemblyGraphEdgeId++) {
        AssemblyGraph::Edge& assemblyGraphEdge = newAssemblyGraph.edges[assemblyGraphEdgeId];
        const span<MarkerGraph::EdgeId> path = newAssemblyGraph.edgeLists[assemblyGraphEdgeId];

        // Vertex coverage.
        // Only internal vertices contribute to this -
        // the first and last vertex don't contribute.
        // If there is only one marker graph edge, there are
        // no internal vertices, and in that case vertex coverage
        // metrics are set to zero.
        if(path.size() == 1) {
            assemblyGraphEdge.minVertexCoverage = 0;
            assemblyGraphEdge.averageVertexCoverage = 0;
            assemblyGraphEdge.maxVertexCoverage = 0;
        } else {
            assemblyGraphEdge.minVertexCoverage = std::numeric_limits<uint32_t>::max();
            assemblyGraphEdge.maxVertexCoverage = 0;
            uint64_t sum = 0;
            for(uint64_t i=1; i<path.size(); i++) {
                const MarkerGraph::EdgeId markerGraphEdgeId = path[i];
                const MarkerGraph::Edge& markerGraphEdge = markerGraph.edges[markerGraphEdgeId];
                const MarkerGraph::VertexId markerGraphVertexId = markerGraphEdge.source;
                const uint32_t coverage = uint32_t(markerGraph.vertexCoverage(markerGraphVertexId));
                assemblyGraphEdge.minVertexCoverage = min(assemblyGraphEdge.minVertexCoverage, coverage);
                assemblyGraphEdge.maxVertexCoverage = max(assemblyGraphEdge.maxVertexCoverage, coverage);
                sum += coverage;
            }
            assemblyGraphEdge.averageVertexCoverage = uint32_t(sum / (path.size() - 1));
        }


        // Edge coverage.
        assemblyGraphEdge.minEdgeCoverage = std::numeric_limits<uint32_t>::max();
        assemblyGraphEdge.maxEdgeCoverage = 0;
        uint64_t sum = 0;
        for(uint64_t i=0; i<path.size(); i++) {
            const MarkerGraph::EdgeId markerGraphEdgeId = path[i];
            const uint32_t coverage = uint32_t(markerGraph.edgeMarkerIntervals.size(markerGraphEdgeId));
            assemblyGraphEdge.minEdgeCoverage = min(assemblyGraphEdge.minEdgeCoverage, coverage);
            assemblyGraphEdge.maxEdgeCoverage = max(assemblyGraphEdge.maxEdgeCoverage, coverage);
            sum += coverage;
        }
        assemblyGraphEdge.averageEdgeCoverage = uint32_t(sum / path.size());
    }



    // Create the marker to assembly table for the detangled marker graph.
    newAssemblyGraph.markerToAssemblyTable.createNew(
        largeDataName("New-MarkerToAssemblyTable"),
        largeDataPageSize);
    newAssemblyGraph.createMarkerToAssemblyTable(markerGraph.edges.size());



    // Now replace the old tangled assembly graph with the detangled one.
    assemblyGraph.remove();
    assemblyGraphPointer = 0;
    newAssemblyGraph.vertices.rename(
        largeDataName("AssemblyGraphVertices"));
    newAssemblyGraph.reverseComplementVertex.rename(
        largeDataName("AssemblyGraphReverseComplementVertex"));
    newAssemblyGraph.edges.rename(
        largeDataName("AssemblyGraphEdges"));
    newAssemblyGraph.edgeLists.rename(
        largeDataName("AssemblyGraphEdgeLists"));
    newAssemblyGraph.edgesBySource.rename(
        largeDataName("AssemblyGraphEdgesBySource"));
    newAssemblyGraph.edgesByTarget.rename(
        largeDataName("AssemblyGraphEdgesByTarget"));
    newAssemblyGraph.reverseComplementEdge.rename(
        largeDataName("AssemblyGraphReverseComplementEdge"));
    newAssemblyGraph.markerToAssemblyTable.rename(
        largeDataName("MarkerToAssemblyTable"));
    assemblyGraphPointer = newAssemblyGraphPointer;


}

