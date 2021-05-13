#include "Mode1-AssemblyGraph.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;
using namespace Mode1;


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
    createVertices(minEdgeCoverage, minEdgeCoveragePerStrand);
}



void Mode1::AssemblyGraph::createVertices(
    uint64_t minEdgeCoverage,
    uint64_t minEdgeCoveragePerStrand)
{
    using EdgeId = MarkerGraph::EdgeId;

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
                cout << "Forward " << edgeId << endl;
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
                    cout << "Backward " << edgeId << endl;
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


        // Also construct the reverse complemented chain.
        reverseComplementedChain.clear();
        for(const EdgeId edgeId: chain) {
            reverseComplementedChain.push_back(markerGraph.reverseComplementEdge[edgeId]);
        }
        std::reverse(reverseComplementedChain.begin(), reverseComplementedChain.end());



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

