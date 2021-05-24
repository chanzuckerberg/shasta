#include "Mode1-AssemblyGraph.hpp"
#include "approximateTopologicalSort.hpp"
#include "findMarkerId.hpp"
#include "Marker.hpp"
#include "orderPairs.hpp"
#include "SHASTA_ASSERT.hpp"
#include "shortestPathBoundedDistance.hpp"
using namespace shasta;
using namespace Mode1;

#include <map>


Mode1::AssemblyGraph::AssemblyGraph(
    uint64_t minEdgeCoverage,
    uint64_t minEdgeCoveragePerStrand,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph) :
    minEdgeCoverage(minEdgeCoverage),
    minEdgeCoveragePerStrand(minEdgeCoveragePerStrand),
    markers(markers),
    markerGraph(markerGraph)
{
    // Create the vertices.
    createVertices(minEdgeCoverage, minEdgeCoveragePerStrand);
    createMarkerGraphToAssemblyGraphTable();

    // Create the edges.
    computePseudoPaths();
    createEdges();
    cout << "The initial assembly graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;

    // Handle jumps.
    const uint64_t maxSkip = 60;   // EXPOSE WHEN CODE STABILIZES
    handleJumps(uint32_t(maxSkip));

    // Transitive reduction.
    uint64_t maxDistanceMarkers = 60;   // EXPOSE WHEN CODE STABILIZES
    transitiveReduction(maxDistanceMarkers);
    cout << "After transitive reduction, the assembly graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;

    // Graphviz output.
    approximateTopologicalSort();
    writeGraphviz("Mode1-AssemblyGraph.dot");
}



void Mode1::AssemblyGraph::createVertices(
    uint64_t minEdgeCoverage,
    uint64_t minEdgeCoveragePerStrand)
{
    using EdgeId = MarkerGraph::EdgeId;
    AssemblyGraph& assemblyGraph = *this;

    // Check that we have what we need.
    SHASTA_ASSERT(markers.isOpen());
    SHASTA_ASSERT(markerGraph.vertices().isOpen());
    SHASTA_ASSERT(markerGraph.edges.isOpen);
    SHASTA_ASSERT(markerGraph.edgeMarkerIntervals.isOpen());
    SHASTA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    SHASTA_ASSERT(markerGraph.reverseComplementEdge.isOpen);

    const bool debug = false;



    // Flags that will be set for marker graph that have already been used.
    const EdgeId edgeCount = markerGraph.edges.size();
    vector<bool> edgeWasUsed(edgeCount, false);



    // Work vectors used in the main loop below.
    vector<EdgeId> nextEdges;
    vector<EdgeId> previousEdges;
    vector<EdgeId> chain;
    vector<EdgeId> reverseComplementedChain;



    // Main loop over all edges of the marker graph.
    // At each iteration we find a new linear chain of edges which
    // generates a vertex (not an edge) of our assembly graph.
    uint64_t totalChainEdgeCount = 0;
    for(EdgeId startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {

        // If this edge has already been used, skip it.
        if(edgeWasUsed[startEdgeId]) {
            continue;
        }

        // If this edge has insufficient coverage, skip it.
        if(not markerGraphEdgeHasSufficientCoverage(startEdgeId)) {
            continue;
        }

        if(debug) {
            cout << "startEdgeId " << startEdgeId << endl;
        }

        // Follow the chain forward.
        edgeWasUsed[startEdgeId] = true;
        EdgeId edgeId = startEdgeId;
        bool isCircularChain = false;
        nextEdges.clear();
        while(true) {
            const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
            const MarkerGraph::VertexId v = edge.target;

            /// Check if we have reached the end of the chain.
            if( (markerGraphVertexInDegree(v)  != 1) or
                (markerGraphVertexOutDegree(v) != 1)) {
                break;
            }

            edgeId = getMarkerGraphUniqueNextEdge(edgeId);
            if(edgeId == MarkerGraph::invalidEdgeId) {
                break;
            }
            if(edgeId == startEdgeId) {
                if(debug) {
                    cout << "circular" << endl;
                }
                isCircularChain = true;
                break;
            }
            nextEdges.push_back(edgeId);
            if(debug) {
                cout << "Forward " << edgeId << " (" <<
                    markerGraph.reverseComplementEdge[edgeId] << ")\n";
            }

            // Mark it as used.
            SHASTA_ASSERT(not edgeWasUsed[edgeId]);
            edgeWasUsed[edgeId] = true;
        }

        // Follow the chain backward.
        previousEdges.clear();
        if(not isCircularChain) {
            edgeId = startEdgeId;
            while(true) {
                const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
                const MarkerGraph::VertexId v = edge.source;

                /// Check if we have reached the end of the chain.
                if( (markerGraphVertexInDegree(v)  != 1) or
                    (markerGraphVertexOutDegree(v) != 1)) {
                    break;
                }

                edgeId = getMarkerGraphUniquePreviousEdge(edgeId);
                if(edgeId == MarkerGraph::invalidEdgeId) {
                    break;
                }
                previousEdges.push_back(edgeId);
                if(debug) {
                    cout << "Backward " << edgeId << " (" <<
                        markerGraph.reverseComplementEdge[edgeId] << ")\n";
                }

                // Mark it as used.
                SHASTA_ASSERT(not edgeWasUsed[edgeId]);
                edgeWasUsed[edgeId] = true;
            }
        }


        // Gather the chain.
        chain.clear();
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter< vector<EdgeId> >(chain));
        chain.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter< vector<EdgeId> >(chain));
        totalChainEdgeCount += chain.size();

        if(debug) {
            for(const EdgeId edgeId: chain) {
                cout << "Chain " << edgeId <<  " (" <<
                    markerGraph.reverseComplementEdge[edgeId] << ")\n";
            }
        }


        // Also construct the reverse complemented chain.
        reverseComplementedChain.clear();
        for(const EdgeId edgeId: chain) {
            reverseComplementedChain.push_back(markerGraph.reverseComplementEdge[edgeId]);
        }
        std::reverse(reverseComplementedChain.begin(), reverseComplementedChain.end());

        if(debug) {
            for(const EdgeId edgeId: reverseComplementedChain) {
                cout << "Reverse complemented chain " << edgeId <<
                    " (" << markerGraph.reverseComplementEdge[edgeId] << ")\n";
            }
        }


        // Figure out if the reverse complemented chain is the same
        // as the original chain. This can happen in exceptional cases.
        bool isSelfComplementary = false;
        if(!isCircularChain) {
            isSelfComplementary = (chain == reverseComplementedChain);
        } else {

            // For a circular chain the test is more complex.
            // We check if the reverse complement of the first edge
            // is in the chain.
            isSelfComplementary =
                find(chain.begin(), chain.end(), reverseComplementedChain.front()) != chain.end();
        }
        if(isSelfComplementary) {
            cout << "Found a self-complementary chain." << endl;
        }

        // Mark as used the edges of the reverse complement chain.
        if(not isSelfComplementary) {
            for(const MarkerGraph::EdgeId edgeId: reverseComplementedChain) {
                SHASTA_ASSERT(not edgeWasUsed[edgeId]);
                edgeWasUsed[edgeId] = true;
            }
            totalChainEdgeCount += reverseComplementedChain.size();
        }


        // Add a vertex corresponding to this chain.
        const vertex_descriptor v = add_vertex(assemblyGraph);
        AssemblyGraphVertex& vertex = assemblyGraph[v];
        vertex.markerGraphEdgeIds = chain;

        // Also add a vertex for the reverse complemented chain.
        if(isSelfComplementary) {
            vertex.vRc = v;
        } else {
            const vertex_descriptor vRc = add_vertex(assemblyGraph);
            AssemblyGraphVertex& vertexRc = assemblyGraph[vRc];
            vertexRc.markerGraphEdgeIds = reverseComplementedChain;
            vertexRc.vRc = v;
            vertex.vRc = vRc;
        }
    }

    cout << "Out of " << edgeCount << " marker graph edges, " << totalChainEdgeCount <<
        " were used to generate the initial assembly graph." << endl;

}


// Return true if a given marker graph edge has sufficient coverage
// (total and for each strand).
bool Mode1::AssemblyGraph::markerGraphEdgeHasSufficientCoverage(MarkerGraph::EdgeId edgeId) const
{
    // Check coverage.
    if(markerGraph.edgeCoverage(edgeId) < minEdgeCoverage) {
        return false;
    }

    // Check coverage for each strand.
    const array<uint64_t, 2> strandCoverage = markerGraph.edgeStrandCoverage(edgeId);
    if( (strandCoverage[0] < minEdgeCoveragePerStrand) or
        (strandCoverage[1] < minEdgeCoveragePerStrand)) {
        return false;
    }

    return true;

}



// Return the out-degree of a marker graph vertex,
// counting only edges with sufficient coverage.
uint64_t Mode1::AssemblyGraph::markerGraphVertexOutDegree(MarkerGraph::VertexId v) const
{
    const auto edgeIds = markerGraph.edgesBySource[v];

    uint64_t outDegree = 0;
    for(const MarkerGraph::EdgeId edgeId: edgeIds) {
        if(markerGraphEdgeHasSufficientCoverage(edgeId)) {
            ++outDegree;
        }
    }
    return outDegree;
}



// Return the in-degree of a marker graph vertex,
// counting only edges with sufficient coverage.
uint64_t Mode1::AssemblyGraph::markerGraphVertexInDegree(MarkerGraph::VertexId v) const
{

    const auto edgeIds = markerGraph.edgesByTarget[v];

    uint64_t inDegree = 0;
    for(const MarkerGraph::EdgeId edgeId: edgeIds) {
        if(markerGraphEdgeHasSufficientCoverage(edgeId)) {
            ++inDegree;
        }
    }
    return inDegree;
}



// Return the unique next/previous marker graph edge for a given marker graph edge,
// or MarkerGraph::invalidEdgeId if there are none or more than one.
// The next/previous marker graph edge is chosen among the
// ones with sufficient coverage.
MarkerGraph::EdgeId Mode1::AssemblyGraph::getMarkerGraphUniqueNextEdge(
    MarkerGraph::EdgeId edgeId01) const
{
    const MarkerGraph::Edge& edge01 = markerGraph.edges[edgeId01];
    MarkerGraph::EdgeId nextEdgeId = MarkerGraph::invalidEdgeId;
    const MarkerGraph::VertexId v1 = edge01.target;

    // Loop over possible next edges.
    for(const MarkerGraph::EdgeId edgeId12: markerGraph.edgesBySource[v1]) {

        // Check coverage.
        if(not markerGraphEdgeHasSufficientCoverage(edgeId12)) {
            continue;
        }

        // This is a possibility.
        if(nextEdgeId == MarkerGraph::invalidEdgeId) {
            // This is the first we found.
            nextEdgeId = edgeId12;
        } else {
            // This is not the first we found.
            return MarkerGraph::invalidEdgeId;
        }
    }

    // If getting here, we only found one.
    return nextEdgeId;
}

MarkerGraph::EdgeId Mode1::AssemblyGraph::getMarkerGraphUniquePreviousEdge(
    MarkerGraph::EdgeId edgeId10) const
{
    const MarkerGraph::Edge& edge10 = markerGraph.edges[edgeId10];
    MarkerGraph::EdgeId previousEdgeId = MarkerGraph::invalidEdgeId;
    const MarkerGraph::VertexId v1 = edge10.source;

    // Loop over possible next edges.
    for(const MarkerGraph::EdgeId edgeId21: markerGraph.edgesByTarget[v1]) {

        // Check coverage.
        if(not markerGraphEdgeHasSufficientCoverage(edgeId21)) {
            continue;
        }

        // This is a possibility.
        if(previousEdgeId == MarkerGraph::invalidEdgeId) {
            // This is the first we found.
            previousEdgeId = edgeId21;
        } else {
            // This is not the first we found.
            return MarkerGraph::invalidEdgeId;
        }
    }

    // If getting here, we only found one.
    return previousEdgeId;
}




// For each marker graph edge, store the Mode1::AssemblyGraph vertex
// that it is on. Can be null_vertex() for marker graph edges not associated
// with a Mode1::AssemblyGraph vertex.
// Indexed by MarkerGraph::EdgeId.
void Mode1::AssemblyGraph::createMarkerGraphToAssemblyGraphTable()
{
    const AssemblyGraph& graph = *this;

    markerGraphToAssemblyGraphTable.clear();
    markerGraphToAssemblyGraphTable.resize(markerGraph.edges.size(), null_vertex());

    BGL_FORALL_VERTICES(v, graph, AssemblyGraph) {
        const AssemblyGraphVertex& vertex = graph[v];
        for(const MarkerGraph::EdgeId edgeId: vertex.markerGraphEdgeIds) {
            SHASTA_ASSERT(markerGraphToAssemblyGraphTable[edgeId] == null_vertex());
            markerGraphToAssemblyGraphTable[edgeId] = v;
        }
    }
}



void Mode1::AssemblyGraph::computePseudoPaths()
{
    // Get the number of reads and oriented reads.
    const ReadId orientedReadCount = ReadId(markers.size());
    SHASTA_ASSERT((orientedReadCount % 2) == 0);
    const ReadId readCount = orientedReadCount / 2;

    // Start with empty pseudo-paths.
    pseudoPaths.clear();
    pseudoPaths.resize(orientedReadCount);

    // Conmpute pseudo-paths of all oriented reads.
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            PseudoPath& pseudoPath = pseudoPaths[orientedReadId.getValue()];
            computePseudoPath(orientedReadId, pseudoPath);
        }
    }
}



void Mode1::AssemblyGraph::computePseudoPath(
    OrientedReadId orientedReadId,
    PseudoPath& pseudoPath)
{

    // Start with an empty pseudo-path.
    pseudoPath.clear();

    // Access the markers for this oriented read.
    const span<const CompressedMarker> orientedReadMarkers = markers[orientedReadId.getValue()];
    const uint32_t markerCount = uint32_t(orientedReadMarkers.size());

    // The first MarkerId for this oriented read.
    const MarkerId markerIdOffset = orientedReadMarkers.begin() - markers.begin();

    // Loop over markers of this oriented read.
    vertex_descriptor vPrevious = null_vertex();
    for(uint32_t ordinal0=0; ordinal0<markerCount; /*increment later*/) {
        const MarkerId markerId0 = markerIdOffset + ordinal0;
        const MarkerGraph::VertexId v0 = markerGraph.vertexTable[markerId0];

        // If there is no vertex for this marker, do nothing and try the next ordinal.
        if(v0 == MarkerGraph::invalidCompressedVertexId) {
            ++ordinal0;
            continue;
        }

        // We found the marker graph vertex corresponding to this marker.
        // Amount its outgoing edges, found the one followed by this oriented read.
        const auto edgeIds = markerGraph.edgesBySource[v0];
        MarkerGraph::EdgeId edgeId01 = MarkerGraph::invalidEdgeId;
        uint32_t ordinal1 = std::numeric_limits<uint32_t>::max();
        for(const MarkerGraph::EdgeId e01: edgeIds) {
            const span<const MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[e01];
            for(const MarkerInterval& markerInterval: markerIntervals) {
                if((markerInterval.orientedReadId == orientedReadId) and
                    markerInterval.ordinals[0] == ordinal0) {
                    // We found it.
                    edgeId01 = e01;
                    ordinal1 = markerInterval.ordinals[1];
                    break;
                }
            }
            if(edgeId01 != MarkerGraph::invalidEdgeId) {
                break;
            }
        }

        // If we did not find it, we must have reached the end of the oriented read.
        // Check that there are no more marker graph vertices on this oriented read.
        if(edgeId01 == MarkerGraph::invalidEdgeId) {
            for(uint32_t ordinal1=ordinal0+1; ordinal1<markerCount; ordinal1++) {
                const MarkerId markerId1 = markerIdOffset + ordinal1;
                const MarkerGraph::VertexId v1 = markerGraph.vertexTable[markerId1];
                SHASTA_ASSERT(v1 == MarkerGraph::invalidCompressedVertexId);
            }
            break;
        }

        // If getting here, we know that marker graph edge edgeId01 contains the
        // marker interval (ordinal0, ordinal1) for this oriented read.
        // Find the Mode1::AssemblyGraph vertecx that this marker graph edge is on
        // and add it to the pseudo-path ,if different from the previous one.
        const vertex_descriptor v = markerGraphToAssemblyGraphTable[edgeId01];
        if(v != null_vertex()) {
            if(v != vPrevious) {
                pseudoPath.push_back(v);
                vPrevious = v;
            }

        }

        // Continue the loop on ordinals.
        ordinal0 = ordinal1;

    }
}



void Mode1::AssemblyGraph::createEdges()
{
    // Get the number of reads and oriented reads.
    const ReadId orientedReadCount = ReadId(markers.size());
    SHASTA_ASSERT((orientedReadCount % 2) == 0);
    const ReadId readCount = orientedReadCount / 2;

    // Loop over oriented reads.
    // Look at pseudo-path transitions.
    std::map< pair<vertex_descriptor, vertex_descriptor>, vector<OrientedReadId> > pathMap;
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);

            // Get its pseudo-path.
            const PseudoPath& pseudoPath = pseudoPaths[orientedReadId.getValue()];

            for(uint64_t i=1; i<pseudoPath.size(); i++) {
                const vertex_descriptor v0 = pseudoPath[i-1];
                const vertex_descriptor v1 = pseudoPath[i];
                pathMap[make_pair(v0,v1)].push_back(orientedReadId);
            }
        }
    }

    // Each entry in the pathMap generates an edge.
    for(const auto& p: pathMap) {
        const vertex_descriptor v0 = p.first.first;
        const vertex_descriptor v1 = p.first.second;
        const vector<OrientedReadId>& orientedReadIds = p.second;
        add_edge(v0, v1, AssemblyGraphEdge(orientedReadIds), *this);
    }
}



// Given an edge e01 v0->v1, return true if the edge corresponds to a "jump"
// in the marker graph. This is the case if the last marker graph vertex
// of the v0 marker graph path is not the same as the first marker graph edge of
// the v1 marker graph path.
bool Mode1::AssemblyGraph::isMarkerGraphJump(edge_descriptor e01) const
{
    const AssemblyGraph& graph = *this;

    // Access the assembly graph vertices.
    const vertex_descriptor v0 = source(e01, graph);
    const vertex_descriptor v1 = target(e01, graph);
    const AssemblyGraphVertex vertex0 = graph[v0];
    const AssemblyGraphVertex vertex1 = graph[v1];

    // Access the marker graph edges we need to check.
    const MarkerGraph::EdgeId lastMarkerGraphEdgeId0  = vertex0.markerGraphEdgeIds.back();
    const MarkerGraph::EdgeId firstMarkerGraphEdgeId1 = vertex1.markerGraphEdgeIds.front();
    const MarkerGraph::Edge lastMarkerGraphEdge0  = markerGraph.edges[lastMarkerGraphEdgeId0];
    const MarkerGraph::Edge firstMarkerGraphEdge1 = markerGraph.edges[firstMarkerGraphEdgeId1];

    // Now it's easy to know if we have a jump.
    return lastMarkerGraphEdge0.target != firstMarkerGraphEdge1.source;
}



// Transitive reduction up to the specified distance, expressed in markers.
void Mode1::AssemblyGraph::transitiveReduction(uint64_t maxDistanceMarkers)
{
    AssemblyGraph& graph = *this;
    const bool debug = false;

    // Gather the edges in order of increasing coverage.
    vector< pair<edge_descriptor, uint64_t> > sortedEdges;
    BGL_FORALL_EDGES(e, graph, AssemblyGraph) {
        sortedEdges.push_back(make_pair(e, graph[e].orientedReadIds.size()));
    }
    sort(sortedEdges.begin(), sortedEdges.end(),
        OrderPairsBySecondOnly<edge_descriptor, uint64_t>());

    // Create an edgeLength map containing for each edge the number of marker graph edges
    // corresponding to the source vertex of that assembly graph edge.
    // This map will be used to find shortest paths below.
    std::map<edge_descriptor, uint64_t> edgeLength;
    BGL_FORALL_EDGES(e, graph, AssemblyGraph) {
        const vertex_descriptor v0 = source(e, graph);
        edgeLength.insert(make_pair(e, graph[v0].markerGraphEdgeIds.size()));
    }

    // Process the edges in order of increasing coverage.
    vector<edge_descriptor> path;
    std::set<edge_descriptor> removedEdges;
    for(const auto& p: sortedEdges) {
        const edge_descriptor e = p.first;
        if(removedEdges.find(e) != removedEdges.end()) {
            continue;
        }
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        if(debug) {
            cout << "Working on edge " << getVertexId(v0) << "->" << getVertexId(v1) << "\n";
        }

        // Temporary changes to the edgeLength map.
        BGL_FORALL_OUTEDGES(v0, e01, graph, AssemblyGraph) {
            edgeLength[e01] = 0;
        }
        edgeLength[e] = std::numeric_limits<uint64_t>::max();

        // Compute the shortest path from v0 to v1, excluding this edge.
        const uint64_t pathLength =
            shortestPathBoundedDistance(graph, v0, v1, maxDistanceMarkers, edgeLength, path);
        if(debug) {
            cout << "Found a path with " << path.size() << " edges and length " << pathLength << ":\n";
            for(const edge_descriptor e: path) {
                const vertex_descriptor v0 = source(e, graph);
                const vertex_descriptor v1 = target(e, graph);
                cout << getVertexId(v0) << "->" << getVertexId(v1) << "\n";
            }
        }

        // Restore the edgeLength.
        BGL_FORALL_OUTEDGES(v0, e01, graph, AssemblyGraph) {
            edgeLength[e01] = graph[v0].markerGraphEdgeIds.size();
        }

        // If we found a path, remove this edge and its reverse complement.
        if(not path.empty()) {
            boost::remove_edge(e, graph);
            removedEdges.insert(e);
            const edge_descriptor eRc = getReverseComplementEdge(e);
            if(eRc != e) {
                boost::remove_edge(eRc, graph);
                removedEdges.insert(eRc);
            }
        }

    }
    cout << "Transitive reduction removed " << removedEdges.size() << " edges." << endl;
}

// Return the reverse complement of an edge.
Mode1::AssemblyGraph::edge_descriptor Mode1::AssemblyGraph::getReverseComplementEdge(
    edge_descriptor e) const
{
    const AssemblyGraph& graph = *this;
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);
    const vertex_descriptor v0Rc = graph[v0].vRc;
    const vertex_descriptor v1Rc = graph[v1].vRc;

    edge_descriptor eRc;
    bool wasFound = false;
    tie(eRc, wasFound) = edge(v1Rc, v0Rc, graph);
    SHASTA_ASSERT(wasFound);
    return eRc;
}



void Mode1::AssemblyGraph::approximateTopologicalSort()
{
    Mode1::AssemblyGraph& graph = *this;

    vector<pair<uint64_t, edge_descriptor> > edgeTable;
    BGL_FORALL_EDGES(e, graph, AssemblyGraph) {
        edgeTable.push_back(make_pair(graph[e].orientedReadIds.size(), e));
    }
    sort(edgeTable.begin(), edgeTable.end(),
        std::greater< pair<uint64_t, edge_descriptor> >());

    vector<edge_descriptor> sortedEdges;
    for(const auto& p: edgeTable) {
        sortedEdges.push_back(p.second);
    }

    shasta::approximateTopologicalSort(graph, sortedEdges);

}



// For edges that correspond to jumps in the marker graph,
// compute marker intervals for a possible secondary marker graph edge.
void Mode1::AssemblyGraph::handleJumps(uint32_t maxSkip)
{
    AssemblyGraph& graph = *this;
    const bool debug = false;

    // Use "a" suffixes for assembly graph vertices and edges.
    // Use "m" suffixes for marker graph vertices and edges.

    // Loop over all edges.
    BGL_FORALL_EDGES(e01a, graph, AssemblyGraph) {

        // If not a marker graph jump, skip it.
        if(not isMarkerGraphJump(e01a)) {
            continue;
        }
        AssemblyGraphEdge& edge01a = graph[e01a];
        edge01a.isMarkerGraphJump = true;

        // Compute marker intervals for a possible secondary marker graph edge.

        // Access the assembly graph vertices.
        const vertex_descriptor v0a = source(e01a, graph);
        const vertex_descriptor v1a = target(e01a, graph);
        const AssemblyGraphVertex& vertex0a = graph[v0a];
        const AssemblyGraphVertex& vertex1a = graph[v1a];

        // Access the marker graph vertices.
        const MarkerGraph::EdgeId edgeId0m = vertex0a.markerGraphEdgeIds.back();
        const MarkerGraph::EdgeId edgeId1m = vertex1a.markerGraphEdgeIds.front();
        const MarkerGraph::Edge& edge0m = markerGraph.edges[edgeId0m];
        const MarkerGraph::Edge& edge1m = markerGraph.edges[edgeId1m];
        const MarkerGraph::VertexId vertexId0m = edge0m.target;
        const MarkerGraph::VertexId vertexId1m = edge1m.source;

        // debug = (vertexId0m == 239798) and (vertexId1m == 235325);

        // Access the markers for the two marker graph vertices.
        const span<const MarkerId> markers0 = markerGraph.getVertexMarkerIds(vertexId0m);
        const span<const MarkerId> markers1 = markerGraph.getVertexMarkerIds(vertexId1m);
        const auto begin0 = markers0.begin();
        const auto begin1 = markers1.begin();
        const auto end0 = markers0.end();
        const auto end1 = markers1.end();

        if(debug) {
            cout << "Working on assembly graph edge " <<
                getVertexId(v0a) << "->" << getVertexId(v1a) << "\n";
            cout << "Marker graph vertices: " << vertexId0m << " " << vertexId1m << "\n";
        }



        // Create marker intervals.
        // For each vertex, the markers are sorted by OrientedReadId and ordinal.
        auto it0 = begin0;
        auto it1 = begin1;
        while((it0!=end0) and (it1!=end1)) {

            // Find the oriented read ids.
            OrientedReadId orientedReadId0;
            OrientedReadId orientedReadId1;
            tie(orientedReadId0, ignore) = findMarkerId(*it0, markers);
            tie(orientedReadId1, ignore) = findMarkerId(*it1, markers);

            // If not the same oriented read id, increment the iterator
            // corresponding to the lower oriented read id.
            if(orientedReadId0 < orientedReadId1) {
                ++it0;
                continue;
            }
            if(orientedReadId1 < orientedReadId0) {
                ++it1;
                continue;
            }
            SHASTA_ASSERT(orientedReadId0 == orientedReadId1);
            const OrientedReadId orientedReadId = orientedReadId0;

            // Each of the two vertices can have more than one marker on the same oriented read.
            // Find the streaks with this oriented read id in markers0 and marker1.
            auto streakBegin0 = it0;
            auto streakEnd0 = streakBegin0 + 1;
            while(streakEnd0 != end0) {
                OrientedReadId o;
                tie(o, ignore) = findMarkerId(*streakEnd0, markers);
                if(o != orientedReadId) {
                    break;
                }
                ++streakEnd0;
            }
            auto streakBegin1 = it1;
            auto streakEnd1 = streakBegin1 + 1;
            while(streakEnd1 != end1) {
                OrientedReadId o;
                tie(o, ignore) = findMarkerId(*streakEnd1, markers);
                if(o != orientedReadId) {
                    break;
                }
                ++streakEnd1;
            }


            // In most cases the streaks are going to be a single marker,
            // but we have to be prepared to handle the general case.
            for(auto jt0=streakBegin0; jt0!=streakEnd0; ++jt0) {
                uint32_t ordinal0;
                tie(ignore, ordinal0) = findMarkerId(*jt0, markers);

                // Find the best ordinal1 for this ordinal0.
                uint32_t bestOrdinal1 = std::numeric_limits<uint32_t>::max();
                for(auto jt1=streakBegin1; jt1!=streakEnd1; ++jt1) {
                    uint32_t ordinal1;
                    tie(ignore, ordinal1) = findMarkerId(*jt1, markers);
                    if(ordinal1 <= ordinal0) {
                        continue;
                    }
                    if(ordinal1 > ordinal0 + maxSkip) {
                        continue;
                    }
                    if(ordinal1 < bestOrdinal1) {
                        bestOrdinal1 = ordinal1;
                    }
                }
                if(bestOrdinal1 == std::numeric_limits<uint32_t>::max()) {
                    continue;
                }
                const uint32_t ordinal1 = bestOrdinal1;
                edge01a.markerIntervals.push_back(MarkerInterval(orientedReadId, ordinal0, ordinal1));
                if(debug) {
                    cout << orientedReadId << " " << ordinal0 << " " << ordinal1 << "\n";
                }
            }

            // Point the iterators to the end of the streaks.
            it0 = streakEnd0;
            it1 = streakEnd1;
        }

        if(debug and edge01a.markerIntervals.empty()) {
            cout << "Unable to add secondary edge at assembly graph vertices " <<
                getVertexId(v0a) << " " << getVertexId(v1a) <<
                " marker graph vertices " <<
                vertexId0m << " " << vertexId1m << endl;
        }

    }

}


/*
Output in Graphviz format.
To display the neighborhood of a vertex:
CreateLocalSubgraph.py Mode1-AssemblyGraph.dot 4024 10
dot -O -T svg LocalSubgraph-Mode1-AssemblyGraph.dot

The first parameter is the Mode1::AssemblyGraph vertex id
(same as the first marker graph edge id).
The second parameter is the distance.
*/

void Mode1::AssemblyGraph::writeGraphviz(const string& fileName) const
{
    ofstream graphOut(fileName);
    writeGraphviz(graphOut);
}
void Mode1::AssemblyGraph::writeGraphviz(ostream& s) const
{
    using Graph = AssemblyGraph;
    const Graph& graph = *this;

    s << "digraph Mode1AssemblyGraph{\n"
        "node [shape=rectangle]";


    // Vertices.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const AssemblyGraphVertex& vertex = graph[v];
        const MarkerGraph::EdgeId firstMarkerGraphEdgeId = vertex.markerGraphEdgeIds.front();
        const MarkerGraph::EdgeId lastMarkerGraphEdgeId  = vertex.markerGraphEdgeIds.back();
        const MarkerGraph::Edge& firstMarkerGraphEdge = markerGraph.edges[firstMarkerGraphEdgeId];
        const MarkerGraph::Edge& lastMarkerGraphEdge  = markerGraph.edges[lastMarkerGraphEdgeId];
        const MarkerGraph::VertexId firstMarkerGraphVertex = firstMarkerGraphEdge.source;
        const MarkerGraph::VertexId lastMarkerGraphVertex  = lastMarkerGraphEdge.target;

        s << getVertexId(v) << " [width=\"" <<
            0.1*sqrt(double(vertex.markerGraphEdgeIds.size())) <<
            "\" label=\"" <<
            "v0 " << firstMarkerGraphVertex << "\\n" <<
            "v1 " << lastMarkerGraphVertex << "\\n" <<
            "e0 " << firstMarkerGraphEdgeId << "\\n" <<
            "e1 " << lastMarkerGraphEdgeId << "\\n" <<
            "n " << vertex.markerGraphEdgeIds.size() <<
            "\"];\n";
    }



    // Edges.
    BGL_FORALL_EDGES(e, graph, Graph) {
        const AssemblyGraphEdge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        s << getVertexId(v0) << "->" <<
            getVertexId(v1) << " [penwidth=\"" <<
            0.3*double(graph[e].orientedReadIds.size()) <<
            "\" label=\"" << graph[e].orientedReadIds.size();
        for(uint64_t i=0; i<edge.orientedReadIds.size(); i++) {
            const OrientedReadId orientedReadId = edge.orientedReadIds[i];
            if((i % 4) == 0) {
                s << "\\n";
            } else {
                s << " ";
            }
            s << orientedReadId;

        }
        s<< "\"";

        if(not graph[e].isDagEdge) {
            s << " constraint=false";
        }

        if(edge.isMarkerGraphJump) {
            if(edge.markerIntervals.empty()) {
                s << " color=red";
            } else {
                s << " color=cyan";
            }
        }

        s << "];\n";
    }

    s << "}\n";
}
