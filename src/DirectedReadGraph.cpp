// Shasta.
#include "DirectedReadGraph.hpp"
#include "LocalDirectedReadGraph.hpp"
#include "Alignment.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



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
void DirectedReadGraph::extractLocalSubgraph(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    LocalDirectedReadGraph& graph)
{
    // The local graph must start empty.
    SHASTA_ASSERT(boost::num_vertices(graph) == 0);

    // Get the vertices in this neighborhood.
    std::map<VertexId, uint64_t> distanceMap;
    findNeighborhood(orientedReadId.getValue(), maxDistance, true, true, distanceMap);

    // Add them to the local subgraph.
    for(const pair<VertexId, uint64_t>& p: distanceMap) {
        const VertexId vertexId = p.first;
        const uint64_t distance = p.second;
        graph.addVertex(OrientedReadId(ReadId(vertexId)),
            getVertex(vertexId).baseCount, distance);
    }



    // Add the edges to the local subgraph.
    using boost::vertices;  // Hide DirectedReadGraph::vertices for BGL_FORALL_VERTICES.

    // Loop over vertices of the local subgraph.
    BGL_FORALL_VERTICES(v0, graph, LocalDirectedReadGraph) {

        // Find the corresponding vertex in the global graph.
        const VertexId vertexId0 = graph[v0].orientedReadId.getValue();
        const OrientedReadId orientedReadId0 = OrientedReadId(OrientedReadId::Int(vertexId0));

        // Find its out-edges.
        const MemoryAsContainer<EdgeId> outEdges0 = outEdges(vertexId0);

        // Loop over the out-edges.
        for(const EdgeId edgeId: outEdges0) {

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
                alignmentInfo.twiceOffsetAtCenter(),
                alignmentInfo.markerCount);
        }
    }

}

