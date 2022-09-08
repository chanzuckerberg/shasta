#include "mode3-JaccardGraph.hpp"
#include "deduplicate.hpp"
#include "mode3.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

#include <boost/graph/iteration_macros.hpp>

#include "fstream.hpp"



void AssemblyGraph::createJaccardGraph(
    size_t threadCount
    )
{
    cout << timestamp << "createJaccardGraph begins." << endl;

    // Compute edges, in parallel.
    jaccardGraph.threadEdges.clear();
    jaccardGraph.threadEdges.resize(threadCount);
    const uint64_t batchSize = 100;
    setupLoadBalancing(markerGraphPaths.size(), batchSize);
    runThreads(&AssemblyGraph::createJaccardGraphThreadFunction, threadCount);



    // Because we do this in both directions, two things can happen:
    // - An edge can be found twice.
    // - Vertex in-degree and out-degree are often equal to 1 but
    // are in principle unlimited.

    // Represent the edges we found using a boost::Graph.
    using Graph = boost::adjacency_list<
        boost::listS, boost::vecS, boost::bidirectionalS,
        boost::no_property, SegmentPairInformation>;
    Graph graph(markerGraphPaths.size());
    for(const auto& threadEdges: jaccardGraph.threadEdges) {
        for(const auto& jaccardGraphEdge: threadEdges) {
            const uint64_t segmentId0 = jaccardGraphEdge.segmentId0;
            const uint64_t segmentId1 = jaccardGraphEdge.segmentId1;
            bool edgeExists = false;
            tie(ignore, edgeExists) = boost::edge(segmentId0, segmentId1, graph);
            if(not edgeExists) {
                add_edge(segmentId0, segmentId1, jaccardGraphEdge.segmentPairInformation, graph);
            }
        }
    }



    // When a vertex as in-degree/out-degree > 1,
    // of all the outgoing/incoming edges we only keep the one
    // with the smallest estimated gap.
    vector<Graph::edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(out_degree(v, graph) > 1) {
            vector<pair<int64_t, Graph::edge_descriptor> > outgoingEdges; // Gap abs value
            BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
                const int64_t offset = graph[e].offset;
                const uint64_t segmentId0 = source(e, graph);
                const int64_t length0 = int64_t(markerGraphPaths.size(segmentId0));
                const int64_t gap = offset - length0;
                outgoingEdges.push_back(make_pair(labs(gap), e));
            }
            sort(outgoingEdges.begin(), outgoingEdges.end());
            for(uint64_t i=1; i<outgoingEdges.size(); i++) {
                edgesToBeRemoved.push_back(outgoingEdges[i].second);
            }
        }
        if(in_degree(v, graph) > 1) {
            vector<pair<int64_t, Graph::edge_descriptor> > incomingEdges; // Gap abs value
            BGL_FORALL_INEDGES(v, e, graph, Graph) {
                const int64_t offset = graph[e].offset;
                const uint64_t segmentId0 = source(e, graph);
                const int64_t length0 = int64_t(markerGraphPaths.size(segmentId0));
                const int64_t gap = offset - length0;
                incomingEdges.push_back(make_pair(labs(gap), e));
            }
            sort(incomingEdges.begin(), incomingEdges.end());
            for(uint64_t i=1; i<incomingEdges.size(); i++) {
                edgesToBeRemoved.push_back(incomingEdges[i].second);
            }
        }
    }
    deduplicate(edgesToBeRemoved);
    for(const auto e: edgesToBeRemoved) {
        remove_edge(e, graph);
    }

    cout << "The Jaccard graph has " << num_vertices(graph) <<
        " vertices (segments) and " << num_edges(graph) << " edges." << endl;

    // Write them out.
    ofstream dot("JaccardGraph.dot");
    dot << "digraph JaccardGraph {" << endl;
    BGL_FORALL_EDGES(e, graph, Graph) {
        dot << source(e, graph) << "->" << target(e, graph) << "\n";
    }
    dot << "}" << endl;
    cout << timestamp << "createJaccardGraph ends." << endl;
}



void AssemblyGraph::createJaccardGraphThreadFunction(size_t threadId)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all segments assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; ++segmentId) {
            createJaccardGraphEdges(segmentId, jaccardGraph.threadEdges[threadId]);
        }
    }
}



void AssemblyGraph::createJaccardGraphEdges(
    uint64_t segmentId,
    vector<JaccardGraphEdge>& edges)
{
    for(uint64_t direction=0; direction<2; direction++) {
        createJaccardGraphEdges(segmentId, direction, edges);
    }
}



// This follows an algorithm similar to the one used by createAssemblyPath3.
void AssemblyGraph::createJaccardGraphEdges(
    uint64_t primarySegmentId,
    uint64_t direction,
    vector<JaccardGraphEdge>& edges)
{
    // EXPOSE WHEN CODE STABILIZES.
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
    JaccardGraphEdge edge;
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

