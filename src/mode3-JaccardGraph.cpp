#include "mode3-JaccardGraph.hpp"
#include "deduplicate.hpp"
#include "mode3.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

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

    // Combine the edges found by all threads.
    vector<JaccardGraphEdge> edges;
    for(const auto& v: jaccardGraph.threadEdges) {
        copy(v.begin(), v.end(), back_inserter(edges));
    }

    // Deduplicate them.
    deduplicate(edges);
    cout << "The Jaccard graph has " << markerGraphPaths.size() <<
        " vertices (segments) and " << edges.size() << " edges." << endl;

    // Write them out.
    ofstream dot("JaccardGraph.dot");
    dot << "digraph JaccardGraph {" << endl;
    for(const auto& edge: edges) {
        dot << edge.segmentId0 << "->" << edge.segmentId1 << "\n";
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

