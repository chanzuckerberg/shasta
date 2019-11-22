// Shasta.
#include "DirectedReadGraph.hpp"
#include "Alignment.hpp"
#include "LocalDirectedReadGraph.hpp"
#include "orderPairs.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "chrono.hpp"
#include "memory.hpp"
#include <queue>


void DirectedReadGraph::createVertices(ReadId readCount)
{
    vertices.resize(2*readCount);
    for(ReadId readId=0; readId<readCount; readId++) {
        const VertexId vertexId0 = OrientedReadId(readId, 0).getValue();
        const VertexId vertexId1 = OrientedReadId(readId, 1).getValue();
        vertices[vertexId0].reverseComplementedVertexId = vertexId1;
        vertices[vertexId1].reverseComplementedVertexId = vertexId0;
    }
}


// Add a pair of edges corresponding to an alignment.
void DirectedReadGraph::addEdgePair(const AlignmentData& alignment)
{
    const bool debug = false;

    const ReadId readId0 = alignment.readIds[0];
    const ReadId readId1 = alignment.readIds[1];
    const bool isSameStrand = alignment.isSameStrand;

    if(debug) {
        cout << "Working on alignment " << readId0 << " " << readId1 << " " <<
            int(isSameStrand) << endl;
        alignment.info.write(cout);
    }

    // Add the first edge.
    OrientedReadId orientedReadId0(readId0, 0);
    OrientedReadId orientedReadId1(readId1, isSameStrand ? 0 : 1);
    AlignmentInfo alignmentInfo = alignment.info;
    const EdgeId edgeId0 = addEdge(orientedReadId0, orientedReadId1, alignmentInfo);

    // Add the second edge.
    orientedReadId0.flipStrand();
    orientedReadId1.flipStrand();
    swap(orientedReadId0, orientedReadId1);
    alignmentInfo.reverseComplement();
    alignmentInfo.swap();
    const EdgeId edgeId1 = addEdge(orientedReadId0, orientedReadId1, alignmentInfo);

    // Store reverse complemented edge ids.
    getEdge(edgeId0).reverseComplementedEdgeId = edgeId1;
    getEdge(edgeId1).reverseComplementedEdgeId = edgeId0;

}



// Add an edge 0->1, reversing the direction if necessary
DirectedReadGraph::EdgeId DirectedReadGraph::addEdge(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    AlignmentInfo alignmentInfo)
{
    const bool debug = false;

    // Get the read ids and strands.
   const ReadId readId0 = orientedReadId0.getReadId();
   const ReadId readId1 = orientedReadId1.getReadId();
   const ReadId strand0 = orientedReadId0.getStrand();
   const ReadId strand1 = orientedReadId1.getStrand();

   // Sanity check: don't allow alignment with self on
   // same strand or opposste strands.
   SHASTA_ASSERT(readId0 != readId1);

   // Compute the offset at center with the current orientation.
   int twiceOffsetAtCenter = alignmentInfo.twiceOffsetAtCenter();



   // Figure out if we need to swap the reads.
   bool swapNeeded = false;
   if(twiceOffsetAtCenter < 0) {

       // The offset is negative. We need to swap.
       swapNeeded = true;

   } else if(twiceOffsetAtCenter == 0) {

       // The offset is zero.
       // We need to break ties in a way that leaves the read graph
       // invariant under reverse complementing.
       // See comments at the beginning of this function.

       if(strand0==0 and strand1==0) {
           // Both are on strand 0. We need readId0 < readId1.
           swapNeeded = readId1 < readId0;
       } else if(strand0==1 and strand1==1) {
           // Both are on strand 1. We need readId1 < readId0.
           swapNeeded = readId0 < readId1;
       } else {
           // The two reads are on opposite strands. We need strand0=0.
           swapNeeded = strand0 == 1;
       }

   } else {

       // The offset is positive. We don't need to swap.
       SHASTA_ASSERT(twiceOffsetAtCenter > 0);
       swapNeeded = false;
   }

   // Do the swap, if necessary.
   if(swapNeeded) {
       swap(orientedReadId0, orientedReadId1);
       alignmentInfo.swap();
   }

   // Sanity check.
   twiceOffsetAtCenter = alignmentInfo.twiceOffsetAtCenter();
   SHASTA_ASSERT(twiceOffsetAtCenter >= 0);

   // Add the edge.
   const EdgeId edgeId = BaseClass::addEdge(
       orientedReadId0.getValue(),
       orientedReadId1.getValue(),
       DirectedReadGraphEdge(alignmentInfo));

   if(debug) {
       cout << "Adding edge " << orientedReadId0 << " -> " << orientedReadId1 << endl;
       alignmentInfo.write(cout);

   }
   return edgeId;
}



// Make sure the graph is invariant under reverse complementing.
void DirectedReadGraph::check()
{
    // Check that if we reverse complement a vertex twice
    // we get the same vertex.
    for(VertexId v0=0; v0<vertices.size(); v0++) {
        const Vertex& vertex0 = getVertex(v0);
        const VertexId v1 = vertex0.reverseComplementedVertexId;
        SHASTA_ASSERT(v1  != v0);
        const Vertex& vertex1 = getVertex(v1);
        SHASTA_ASSERT(vertex1.reverseComplementedVertexId == v0);
    }

    // Check that if we reverse complement an edge twice
    // we get the same edge.
    // Also check that the offset of the reverse complemented edge
    // is the same.
    for(EdgeId e0=0; e0<edges.size(); e0++) {
        const Edge& edge0 = getEdge(e0);
        const EdgeId e1 = edge0.reverseComplementedEdgeId;
        SHASTA_ASSERT(e1 != e0);
        const Edge& edge1 = getEdge(e1);
        SHASTA_ASSERT(edge0.reverseComplementedEdgeId == e1);
        SHASTA_ASSERT(edge0.alignmentInfo.twiceOffsetAtCenter() ==
            edge1.alignmentInfo.twiceOffsetAtCenter());

        // Also check the vertices.
        SHASTA_ASSERT(source(e0) == getVertex(target(e1)).reverseComplementedVertexId);
        SHASTA_ASSERT(target(e0) == getVertex(source(e1)).reverseComplementedVertexId);
    }

}



// Create a LocalDirectedReadGraph.
bool DirectedReadGraph::extractLocalSubgraph(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    uint64_t minAlignedMarkerCount,
    uint64_t maxOffsetAtCenter,
    double minAlignedFraction,
    float minTransitiveCoverage,
    bool allowTransitiveReductionEdges,
    double timeout,
    LocalDirectedReadGraph& graph)
{
    const auto startTime = steady_clock::now();

    // The local graph must start empty.
    SHASTA_ASSERT(boost::num_vertices(graph) == 0);

    // Construct our edge filter.
    EdgeFilter edgeFilter(minAlignedMarkerCount,
        2*maxOffsetAtCenter,
        minAlignedFraction, minTransitiveCoverage,
        allowTransitiveReductionEdges);

    // Get the vertices in this neighborhood.
    std::map<VertexId, uint64_t> distanceMap;
    if(not findNeighborhood(orientedReadId.getValue(), maxDistance,
        edgeFilter,
        true, true, timeout,
        distanceMap)) {
        graph.clear();
        return false;
    }

    // Add them to the local subgraph.
    for(const pair<VertexId, uint64_t>& p: distanceMap) {
        const VertexId vertexId = p.first;
        const uint64_t distance = p.second;
        const DirectedReadGraphVertex& vertex = getVertex(vertexId);
        graph.addVertex(OrientedReadId(ReadId(vertexId)),
            vertex.baseCount, vertex.markerCount, distance);
    }



    // Add the edges to the local subgraph.
    using boost::vertices;  // Hide DirectedReadGraph::vertices for BGL_FORALL_VERTICES.

    // Loop over vertices of the local subgraph.
    BGL_FORALL_VERTICES(v0, graph, LocalDirectedReadGraph) {

        // See if we exceeded the timeout.
        if(timeout>0. && (seconds(steady_clock::now() - startTime) > timeout)) {
            graph.clear();
            return false;
        }

        // Find the corresponding vertex in the global graph.
        const VertexId vertexId0 = graph[v0].orientedReadId.getValue();
        const OrientedReadId orientedReadId0 = OrientedReadId(OrientedReadId::Int(vertexId0));

        // Find its out-edges.
        const MemoryAsContainer<EdgeId> outEdges0 = outEdges(vertexId0);

        // Loop over the out-edges.
        for(const EdgeId edgeId: outEdges0) {
            if(not edgeFilter.allowEdge(edgeId, getEdge(edgeId))) {
                continue;
            }

            // Find the target vertex.
            const VertexId vertexId1 = target(edgeId);

            // Find the corresponding local vertex, if any.
            const OrientedReadId orientedReadId1 = OrientedReadId(OrientedReadId::Int(vertexId1));
            const LocalDirectedReadGraph::vertex_descriptor v1 =
                graph.getVertex(orientedReadId1);

            // If no such local vertex, skip.
            if(v1 == LocalDirectedReadGraph::null_vertex()) {
                continue;
            }

            // Add the edge to the local subgraph.
            const AlignmentInfo& alignmentInfo = getEdge(edgeId).alignmentInfo;
            graph.addEdge(orientedReadId0, orientedReadId1,
                alignmentInfo,
                getEdge(edgeId).wasRemovedByTransitiveReduction == 1,
                getEdge(edgeId).transitiveCoverage);
        }
    }

    return true;
}



void DirectedReadGraph::transitiveReduction(
    double offsetTolerance0,
    double offsetTolerance1)
{
    const bool debug = false;
    using Edge = DirectedReadGraphEdge;

    // Mark all edges as not removed by transitive reduction.
    for(EdgeId edgeId=0; edgeId<edges.size(); edgeId++) {
        Edge& edge = getEdge(edgeId);
        edge.wasRemovedByTransitiveReduction = 0;
        edge.transitiveCoverage = 1.;
    }

    // Sort edges by decreasing offset.
    vector< pair<EdgeId, uint64_t> > edgeTable;
    for(EdgeId edgeId=0; edgeId<edges.size(); edgeId++) {
        Edge& edge = getEdge(edgeId);
        edgeTable.push_back(make_pair(edgeId, edge.alignmentInfo.twiceOffsetAtCenter()));
    }
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnlyGreater<EdgeId, uint64_t>());



    // Priority queue used by the BFS below.
    // It stores pairs (VertexId, distance),
    // where distance is sum of twice offsets at centers.
    std::priority_queue<
        pair<VertexId, uint64_t>,
        vector< pair<VertexId, uint64_t> >,
        OrderPairsBySecondOnlyGreater<VertexId, uint64_t>
        > q;

    // The vertices visited by the current BFS.
    vector<VertexId> visitedVertices;

    // Flags set for vertices that were visited in the current BFS.
    vector<bool> wasVisited(vertices.size(), false);

    // The edge that led us to each visited vertex.
    vector<EdgeId> predecessorEdge(vertices.size(), invalidEdgeId);

    // The path v0->...->v1.
    vector<EdgeId> path;



    // Process all edges in order of decreasing offset.
    for(const auto& p: edgeTable) {

        // Get some information about this edge.
        const EdgeId edgeId = p.first;
        const uint64_t twiceOffsetAtCenter = p.second;
        const Edge& edge = getEdge(edgeId);
        if(edge.wasRemovedByTransitiveReduction) {
            // This edge was already removed. This can happen
            // because every time we remove an edge we also remove
            // the reverse complemented edge.
            continue;
        }
        const VertexId u0 = source(edgeId);
        const VertexId u1 = target(edgeId);
        const OrientedReadId orientedReadId0 = OrientedReadId(OrientedReadId::Int(u0));
        const OrientedReadId orientedReadId1 = OrientedReadId(OrientedReadId::Int(u1));

        if(debug) {
            cout << "Working on edge " << edgeId << " " << orientedReadId0 << "->" << orientedReadId1 <<
                " with offset " << double(twiceOffsetAtCenter)/2. << endl;
        }

        // Look for the shortest path between vertex0 and vertex1 that:
        // - Does not use edge vertex0->vertex1.
        // - Has total offset sufficiently similar to the offset of edge vertex0->vertex1.

        // Compute the allowable range for twice the total offset of the path.
        const uint64_t offsetToleranceOnTwiceOffset = uint64_t(
            2. * offsetTolerance0 + offsetTolerance1 * double(twiceOffsetAtCenter));
        const uint64_t minPathTwiceOffset =
            (offsetToleranceOnTwiceOffset < twiceOffsetAtCenter) ?
            (twiceOffsetAtCenter - offsetToleranceOnTwiceOffset) :
            0;
        const uint64_t maxPathTwiceOffset =  twiceOffsetAtCenter + offsetToleranceOnTwiceOffset;
        if(debug) {
            cout << "Allowable offset range for path: " << minPathTwiceOffset/2 << " " << maxPathTwiceOffset/2 << endl;
        }



        // Do a forward BFS starting at u0.
        // Use twice offset at center as edge length.
        SHASTA_ASSERT(visitedVertices.empty());
        SHASTA_ASSERT(q.empty());
        q.push(make_pair(u0, 0));
        visitedVertices.push_back(u0);
        wasVisited[u0] = true;
        bool done = false;
        while(not q.empty()) {

            // Dequeue a vertex.
            const auto& p = q.top();
            const VertexId v0 = p.first;
            const uint64_t distance0 = p.second;
            q.pop();
            if(debug) {
                cout << "Dequeued " << OrientedReadId(OrientedReadId::Int(v0)) << " at distance " << distance0/2 << endl;
            }

            // Loop over its out-edges, skipping the ones
            // that were already removed.
            for(const EdgeId edgeId01: outEdges(v0)) {
                if(debug) {
                    cout << "Encountered edge " << edgeId01 << endl;
                }

                // If this is the edge we are working on, skip it.
                if(edgeId01 == edgeId) {
                    continue;
                }
                const Edge& edge = getEdge(edgeId01);

                // If this edge was already removed during transitive reduction, skip it.
                if(edge.wasRemovedByTransitiveReduction) {
                    continue;
                }

                // Get the target vertex.
                const VertexId v1 = target(edgeId01);

                if(debug) {
                    cout << "Found " << OrientedReadId(OrientedReadId::Int(v1)) << endl;
                }

                // If already visited, skip.
                if(wasVisited[v1]) {
                    if(debug) {
                        cout << "Already visited." << endl;
                    }
                    continue;
                }

                // Compute the updated distance.
                const uint64_t distance1 = distance0 + edge.alignmentInfo.twiceOffsetAtCenter();
                if(debug) {
                    cout << "distance1 " << distance1 << endl;
                }

                // If we got too far, skip.
                if(distance1 > maxPathTwiceOffset) {
                    if(debug) {
                        cout << "Distance is too large " << distance1 << endl;
                    }
                    continue;
                }


                // Did we find u1?
                if(v1==u1) {
                    if(distance1 >= minPathTwiceOffset) {
                        Edge& edge = getEdge(edgeId);
                        Edge& reverseComplementedEdge = getEdge(edge.reverseComplementedEdgeId);
                        edge.wasRemovedByTransitiveReduction = true;
                        reverseComplementedEdge.wasRemovedByTransitiveReduction = true;
                        done = true;
                        if(debug) {
                            cout << "Transitive reduction removed edge " << edgeId << " " <<
                                orientedReadId0 << "->" << orientedReadId1 <<
                                " with offset " << double(twiceOffsetAtCenter)/2 << endl;
                        }

                        // To update transitive coverage, we need to reconstruct the path
                        // we found between v0 and v1.
                        SHASTA_ASSERT(path.empty());
                        EdgeId pathEdgeId = edgeId01;
                        while(true) {
                            if(debug) {
                                cout << "Found path edge " << pathEdgeId << endl;
                            }
                            SHASTA_ASSERT(pathEdgeId != invalidEdgeId);
                            SHASTA_ASSERT(not getEdge(pathEdgeId).wasRemovedByTransitiveReduction);
                            path.push_back(pathEdgeId);
                            const VertexId v = source(pathEdgeId);
                            if(v == u0) {
                                break;
                            }
                            pathEdgeId = predecessorEdge[v];
                        }

                        // Donate transitive coverage to the path edges.
                        const float coverageIncrement = edge.transitiveCoverage / float(path.size());
                        // const float coverageIncrement = edge.transitiveCoverage; // EXPERIMENT
                        edge.transitiveCoverage = 0.;
                        reverseComplementedEdge.transitiveCoverage = 0.;
                        for(const EdgeId pathEdgeId: path) {
                            SHASTA_ASSERT(pathEdgeId != invalidEdgeId);
                            Edge& pathEdge = getEdge(pathEdgeId);
                            SHASTA_ASSERT(not pathEdge.wasRemovedByTransitiveReduction);
                            pathEdge.transitiveCoverage += coverageIncrement;
                            getEdge(pathEdge.reverseComplementedEdgeId).transitiveCoverage += coverageIncrement;
                        }
                        path.clear();


                        break;
                    } else {
                        // We found u1, but the distance is too small. Keep going.
                        if(debug) {
                            cout << "Found u1, but the distance is too small." << endl;
                        }
                        continue;
                    }
                } else {

                    SHASTA_ASSERT(v1 != u1);
                    SHASTA_ASSERT(!wasVisited[v1]);
                    SHASTA_ASSERT(predecessorEdge[v1] == invalidEdgeId);
                    wasVisited[v1] = true;
                    predecessorEdge[v1] = edgeId01;
                    visitedVertices.push_back(v1);
                    q.push(make_pair(v1, distance1));
                    if(debug) {
                        cout << "Enqueued " << OrientedReadId(OrientedReadId::Int(v1)) << " at distance " << distance1/2 << endl;
                    }
                }
            }

            if(done) {
                break;
            }
        }


        // Clean up after this BFS.
        for(const VertexId v: visitedVertices) {
            wasVisited[v] = false;
            predecessorEdge[v] = invalidEdgeId;
        }
        visitedVertices.clear();
        while(not q.empty()) {
            q.pop();
        }

    }

    if(debug) {
        cout << "Edges removed during transitive reduction:" << endl;
        for(EdgeId edgeId=0; edgeId<edges.size(); edgeId++) {
            if(getEdge(edgeId).wasRemovedByTransitiveReduction) {
                cout << edgeId << " ";
            }
        }
        cout << endl;
    }
}
