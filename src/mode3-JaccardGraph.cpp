#include "mode3-JaccardGraph.hpp"
#include "deduplicate.hpp"
#include "mode3.hpp"
#include "orderVectors.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
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
    const uint64_t n = markerGraphPaths.size();
    Graph graph(n);
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



    // When a vertex has multiple outgoing/incoming edges,
    // remove the ones with very large gap.
    vector<Graph::edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {

        // Outgoing.
        if(out_degree(v, graph) > 1) {
            int64_t minGap = std::numeric_limits<int64_t>::max();
            BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
                const int64_t offset = graph[e].offset;
                const uint64_t segmentId0 = source(e, graph);
                const int64_t length0 = int64_t(markerGraphPaths.size(segmentId0));
                const int64_t gap = offset - length0;
                minGap = min(minGap, labs(gap));
            }
            const int64_t gapThreshold = 20 + 3 * minGap;  // *** EXPOSE WHEN CODE STABILIZES
            BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
                const int64_t offset = graph[e].offset;
                const uint64_t segmentId0 = source(e, graph);
                const int64_t length0 = int64_t(markerGraphPaths.size(segmentId0));
                const int64_t gap = offset - length0;
                if(gap > gapThreshold) {
                    edgesToBeRemoved.push_back(e);
                }
            }
        }

        // Incoming.
        if(in_degree(v, graph) > 1) {
            int64_t minGap = std::numeric_limits<int64_t>::max();
            BGL_FORALL_INEDGES(v, e, graph, Graph) {
                const int64_t offset = graph[e].offset;
                const uint64_t segmentId0 = source(e, graph);
                const int64_t length0 = int64_t(markerGraphPaths.size(segmentId0));
                const int64_t gap = offset - length0;
                minGap = min(minGap, labs(gap));
            }
            const int64_t gapThreshold = 20 + 3 * minGap;  // *** EXPOSE WHEN CODE STABILIZES
            BGL_FORALL_INEDGES(v, e, graph, Graph) {
                const int64_t offset = graph[e].offset;
                const uint64_t segmentId0 = source(e, graph);
                const int64_t length0 = int64_t(markerGraphPaths.size(segmentId0));
                const int64_t gap = offset - length0;
                if(gap > gapThreshold) {
                    edgesToBeRemoved.push_back(e);
                }
            }
        }
    }
    deduplicate(edgesToBeRemoved);
    for(const auto e: edgesToBeRemoved) {
        remove_edge(e, graph);
    }

    cout << "The Jaccard graph has " << num_vertices(graph) <<
        " vertices (segments) and " << num_edges(graph) << " edges." << endl;

    // Write out the Jaccard graph
    {
        ofstream dot("JaccardGraph.dot");
        dot << "digraph JaccardGraph {" << endl;
        BGL_FORALL_EDGES(e, graph, Graph) {
            dot << source(e, graph) << "->" << target(e, graph) << "\n";
        }
        dot << "}" << endl;
    }

#if 0
    // At this point no vertex has in-degree or out-degree greater than 1.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        SHASTA_ASSERT(out_degree(v, graph) <= 1);
        SHASTA_ASSERT(in_degree(v, graph) <= 1);
    }

    // Compute connected components.
    // Since no vertex has in-degree or out-degree greater than 1,
    // all connected component are guaranteed to have one of two structures:
    // 1. A linear structure beginning at a vertex with in-degree 0,
    //    out-degree 1, ending at a vertex with in-degree 1, out-degree 0,
    //    and other wise consisting of vertices with in-degree and out-degree 1.
    //    This includes isolated vertices.
    // 2. A circular structure consisting entirely of vertices with in-degree and out-degree 1.
    //    This is not common.
    vector<bool> wasFound(n, false);
    vector< vector<uint64_t> > components;
    vector<bool> isCircularComponent;
    for(uint64_t startSegmentId=0; startSegmentId<n; startSegmentId++) {
        if(wasFound[startSegmentId]) {
            continue;
        }
        // cout << "Start a new component at " << startSegmentId << endl;

        // Start a new connected component here,
        wasFound[startSegmentId] = true;
        uint64_t segmentId0 = startSegmentId;
        vector<uint64_t> forwardVertices;
        vector<uint64_t> backwardVertices;
        bool isCircular = false;

        // Look forward.
        while(true) {
            Graph::out_edge_iterator begin, end;
            tie(begin, end) = out_edges(segmentId0, graph);
            if(end == begin) {
                // segmentId0 has out-degree 0
                break;
            }
            const Graph::edge_descriptor e = *begin;
            const uint64_t segmentId1 = target(e, graph);
            // cout << "Forward found " << segmentId1 << endl;
            if(segmentId1 == startSegmentId) {
                isCircular = true;
                break;
            }
            SHASTA_ASSERT(not wasFound[segmentId1]);
            wasFound[segmentId1] = true;
            forwardVertices.push_back(segmentId1);
            segmentId0 = segmentId1;
        }

        // Look backward.
        if(not isCircular) {
            segmentId0 = startSegmentId;
            while(true) {
                Graph::in_edge_iterator begin, end;
                tie(begin, end) = in_edges(segmentId0, graph);
                if(end == begin) {
                    // segmentId0 has in-degree 0
                    break;
                }
                const Graph::edge_descriptor e = *begin;
                const uint64_t segmentId1 = source(e, graph);
                // cout << "Backward found " << segmentId1 << endl;
                SHASTA_ASSERT(segmentId1 != startSegmentId);
                SHASTA_ASSERT(not wasFound[segmentId1]);
                wasFound[segmentId1] = true;
                backwardVertices.push_back(segmentId1);
                segmentId0 = segmentId1;
            }
        }

        // Store it.
        components.resize(components.size() + 1);
        isCircularComponent.push_back(isCircular);
        vector<uint64_t>& component = components.back();
        copy(backwardVertices.rbegin(), backwardVertices.rend(), back_inserter(component));
        component.push_back(startSegmentId);
        copy(forwardVertices.begin(), forwardVertices.end(), back_inserter(component));
    }
    sort(components.begin(), components.end(), OrderVectorsByDecreasingSize<uint64_t>());

    // Write out the components.
    ofstream csv("JaccardGraphComponents.csv");
    csv << "Id,Size,Circular,SegmentIds\n";
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        csv << componentId << ",";
        const vector<uint64_t>& component = components[componentId];
        csv << component.size() << ",";
        if(isCircularComponent[componentId]) {
            csv << "Yes,";
        } else {
            csv << "No,";
        }
        for(const uint64_t segmentId: component) {
            csv << segmentId << ",";
        }
        csv << "\n";
    }

    // Write a histogram of component sizes.
    {
        vector<uint64_t> histogram;
        for(const vector<uint64_t>& component: components) {
            const uint64_t componentSize = component.size();
            if(componentSize >= histogram.size()) {
                histogram.resize(componentSize + 1, 0);
            }
            ++histogram[componentSize];
        }
        ofstream csv("JaccardGraphComponentSizeHistogram.csv");
        csv << "Size,Frequency,Segments,Cumulative segments\n";
        uint64_t cumulativeSegments = 0;
        for(uint64_t componentSize=histogram.size()-1; componentSize>0; componentSize--) {
            const uint64_t frequency = histogram[componentSize];
            const uint64_t segments = componentSize * frequency;
            cumulativeSegments += segments;
            if(frequency > 0) {
                csv << componentSize << ",";
                csv << frequency << ",";
                csv << segments << ",";
                csv << cumulativeSegments << "\n";
            }
        }
    }
#endif

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

