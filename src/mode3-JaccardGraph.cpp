#include "mode3-JaccardGraph.hpp"
#include "deduplicate.hpp"
#include "mode3.hpp"
#include "orderVectors.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

#include "fstream.hpp"



// Create a JaccardGraph with the given number of vertices
// (one for each segment) and no edges.
JaccardGraph::JaccardGraph(uint64_t segmentCount)
{
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        vertexTable.push_back(add_vertex(JaccardGraphVertex(segmentId), *this));
    }
}


void AssemblyGraph::createJaccardGraph(
    size_t threadCount
    )
{
    cout << timestamp << "createJaccardGraph begins." << endl;

    // Create the JaccardGraph and its vertices.
    const uint64_t segmentCount = markerGraphPaths.size();
    jaccardGraphPointer = make_shared<JaccardGraph>(segmentCount);
    JaccardGraph& jaccardGraph = *jaccardGraphPointer;

    // Compute edges, in parallel.
    jaccardGraph.threadEdges.resize(threadCount);
    const uint64_t batchSize = 100;
    setupLoadBalancing(segmentCount, batchSize);
    runThreads(&AssemblyGraph::createJaccardGraphThreadFunction, threadCount);



    // Store the edges in the JaccardGraph.
    for(const auto& threadEdges: jaccardGraph.threadEdges) {
        for(const JaccardGraphEdgeInfo& info: threadEdges) {
            const uint64_t segmentId0 = info.segmentId0;
            const uint64_t segmentId1 = info.segmentId1;
            const JaccardGraph::vertex_descriptor v0 = jaccardGraph.vertexTable[segmentId0];
            const JaccardGraph::vertex_descriptor v1 = jaccardGraph.vertexTable[segmentId1];

            JaccardGraph::edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(v0, v1, jaccardGraph);
            if(not edgeExists) {
                boost::add_edge(v0, v1,
                    JaccardGraphEdge(info.segmentPairInformation, info.direction),
                    jaccardGraph);
            } else {
                jaccardGraph[e].wasFoundInDirection[info.direction] = true;
            }
        }
    }
    jaccardGraph.threadEdges.clear();



    // A "strong edge" is an edge that was found in both directions.
    // An "almost isolated" vertex is a vertex that has no
    // incoming/outgoing strong edges.
    // Remove all edges to "almost isolated" vertices.
    {
        vector<JaccardGraph::edge_descriptor> edgesToBeRemoved;
        BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {

            // Figure out if this vertex has strong edges.
            bool hasStrongEdges = false;
            BGL_FORALL_OUTEDGES(v, e, jaccardGraph, JaccardGraph) {
                if(jaccardGraph[e].wasFoundInBothDirections()) {
                    hasStrongEdges = true;
                    break;
                }
            }
            if(not hasStrongEdges) {
                BGL_FORALL_INEDGES(v, e, jaccardGraph, JaccardGraph) {
                    if(jaccardGraph[e].wasFoundInBothDirections()) {
                        hasStrongEdges = true;
                        break;
                    }
                }
            }

            // If not, flag its edges for removal.
            if(not hasStrongEdges) {
                BGL_FORALL_OUTEDGES(v, e, jaccardGraph, JaccardGraph) {
                    edgesToBeRemoved.push_back(e);
                }
                BGL_FORALL_INEDGES(v, e, jaccardGraph, JaccardGraph) {
                    edgesToBeRemoved.push_back(e);
                }
            }
        }

        // Remove the edges we flagged.
        deduplicate(edgesToBeRemoved);
        for(const JaccardGraph::edge_descriptor e: edgesToBeRemoved) {
            remove_edge(e, jaccardGraph);
        }
    }




    cout << "The Jaccard graph has " << num_vertices(jaccardGraph) <<
        " vertices (segments) and " << num_edges(jaccardGraph) << " edges." << endl;

    // Write out the Jaccard graph.
    // This only writes the edges, so isolated vertices are not included.
    {
        ofstream dot("JaccardGraph.dot");
        dot << "digraph JaccardGraph {" << endl;
        BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
            const JaccardGraphEdge& edge = jaccardGraph[e];
            const JaccardGraph::vertex_descriptor v0 = source(e, jaccardGraph);
            const JaccardGraph::vertex_descriptor v1 = target(e, jaccardGraph);
            const uint64_t segmentId0 = jaccardGraph[v0].segmentId;
            const uint64_t segmentId1 = jaccardGraph[v1].segmentId;
            dot << segmentId0 << "->" << segmentId1;
            if(edge.wasFoundInDirection[0]) {
                if(edge.wasFoundInDirection[1]) {
                    // Found in both directions, leave black.
                } else {
                    // Only found in the forward direction.
                    dot << "[color=red]";
                }
            } else {
                if(edge.wasFoundInDirection[1]) {
                    // Only found in the backward direction.
                    dot << "[color=green]";
                } else {
                    SHASTA_ASSERT(0);
                }
            }
            dot << ";\n";
        }
        dot << "}" << endl;
    }

    cout << timestamp << "createJaccardGraph ends." << endl;
}



void AssemblyGraph::createJaccardGraphThreadFunction(size_t threadId)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all segments assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; ++segmentId) {
            createJaccardGraphEdges(segmentId, jaccardGraphPointer->threadEdges[threadId]);
        }
    }
}



void AssemblyGraph::createJaccardGraphEdges(
    uint64_t segmentId,
    vector<JaccardGraphEdgeInfo>& edges)
{
    for(uint64_t direction=0; direction<2; direction++) {
        createJaccardGraphEdges(segmentId, direction, edges);
    }
}



// This follows an algorithm similar to the one used by createAssemblyPath3.
void AssemblyGraph::createJaccardGraphEdges(
    uint64_t primarySegmentId,
    uint64_t direction,
    vector<JaccardGraphEdgeInfo>& edges)
{
    // EXPOSE WHEN CODE STABILIZES.
    // FOR NOW THESE SHOULD BE THE SAME AS IN AssemblyGraph::createAssemblyPath3.
    const uint64_t minCommonForLink = 3;
    const uint64_t minCommonForPrimary = 3;
    const double minJaccard = 0.7;
    const int32_t minLinkSeparation = -20;

    // We start from primarySegmentId
    // and move in the specified direction until we find segmentId1 with
    // sufficiently high Jaccard similarity and number of
    // common oriented reads with primarySegmentId.
    // At each step, we choose the link that has the most common oriented
    // reads with the primarySegmentId.
    SegmentOrientedReadInformation infoPrimary;
    getOrientedReadsOnSegment(primarySegmentId, infoPrimary);
    JaccardGraphEdgeInfo edge;
    edge.direction = direction;
    uint64_t segmentId0 = primarySegmentId;
    std::set<uint64_t> previousSegments;
    while(true) {

        // Loop over outgoing or incoming links of segmentId0.
        // Find the link with the most common reads with the primarySegmentId.
        const auto linkIds = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
        if(linkIds.empty()) {
            return;
        }
        uint64_t linkIdBest = invalid<uint64_t>;
        uint64_t commonOrientedReadCountBest = 0;
        for(const uint64_t linkId: linkIds) {

            // If link separation is too negative, skip it.
            // The goal here is to avoid cycles in paths.
            const Link& link = links[linkId];
            if(link.separation < minLinkSeparation) {
                continue;
            }

            // Count the number of common oriented reads between the reference segment and this link.
            uint64_t commonOrientedReadCount;
            analyzeSegmentLinkPair(primarySegmentId, linkId, commonOrientedReadCount);

            // If better than the one we have it, record it.
            if(commonOrientedReadCount > commonOrientedReadCountBest) {
                linkIdBest = linkId;
                commonOrientedReadCountBest = commonOrientedReadCount;
            }
        }
        if(commonOrientedReadCountBest < minCommonForLink) {
            return;
        }
        const uint64_t linkId = linkIdBest;

        // Get the segment at the other side of this link.
        const Link& link = links[linkId];
        const uint64_t segmentId1 = (direction==0) ? link.segmentId1 : link.segmentId0;

        // Check that we haven't been here before.
        if(previousSegments.contains(segmentId1)) {
            break;
        }
        previousSegments.insert(segmentId1);

        // Check segmentId1 against the primary segment.
        SegmentOrientedReadInformation info1;
        getOrientedReadsOnSegment(segmentId1, info1);
        if(direction == 0) {
            analyzeSegmentPair(
                    primarySegmentId, segmentId1,
                    infoPrimary, info1,
                    markers, edge.segmentPairInformation);
        } else {
            analyzeSegmentPair(
                segmentId1, primarySegmentId,
                info1, infoPrimary,
                markers, edge.segmentPairInformation);
        }

        // If the Jaccard similarity is high, we found the Jaccard graph edge
        // we were looking for.
        if( edge.segmentPairInformation.commonCount >= minCommonForPrimary and
            edge.segmentPairInformation.jaccard() >= minJaccard) {
            if(direction == 0) {
                edge.segmentId0 = primarySegmentId;
                edge.segmentId1 = segmentId1;
            } else {
                edge.segmentId0 = segmentId1;
                edge.segmentId1 = primarySegmentId;
            }
            edges.push_back(edge);
            return;
        }

        segmentId0 = segmentId1;
    }
}

