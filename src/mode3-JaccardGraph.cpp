#include "mode3-JaccardGraph.hpp"
#include "deduplicate.hpp"
#include "mode3.hpp"
#include "orderPairs.hpp"
#include "orderVectors.hpp"
#include "timestamp.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
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
    cout << "The total number of segments in the assembly graph is " << segmentCount << endl;
    jaccardGraphPointer = make_shared<JaccardGraph>(segmentCount);
    JaccardGraph& jaccardGraph = *jaccardGraphPointer;

    // Compute edges, in parallel.
    jaccardGraph.threadEdges.resize(threadCount);
    const uint64_t batchSize = 100;
    setupLoadBalancing(segmentCount, batchSize);
    runThreads(&AssemblyGraph::createJaccardGraphThreadFunction, threadCount);
    jaccardGraph.storeEdges();
    jaccardGraph.writeGraphviz("JaccardGraph.dot", false, false);
    jaccardGraph.writeGraphviz("JaccardGraph-Labeled.dot", false, true);
    jaccardGraph.writeEdgesCsv("JaccardGraphEdges.csv");
    cout << "The Jaccard graph has " << num_vertices(jaccardGraph) <<
        " vertices (segments) and " << num_edges(jaccardGraph) << " edges." << endl;

#if 0
    // Remove all weak vertices.
    jaccardGraph.removeWeakVertices();
    cout << "After removing weak vertices, the Jaccard graph has " << num_vertices(jaccardGraph) <<
        " vertices (segments) and " << num_edges(jaccardGraph) << " edges." << endl;
    jaccardGraph.writeGraphviz("JaccardGraph1.dot", false, false);
    jaccardGraph.writeGraphviz("JaccardGraph1-Labeled.dot", false, true);
#endif

    // Store the cluster id of each segment (uninitialized for now).
    createNew(clusterIds, "Mode3-ClusterIds");
    jaccardGraph.findClusters(segmentCount, clusterIds);

    // Create the ExpandedJaccardGraph.
    ExpandedJaccardGraph expandedJaccardGraph(jaccardGraph);
    expandedJaccardGraph.writeGraphviz("ExpandedJaccardGraph.dot");

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
    const double minJaccard = 0.75;
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
            edge.segmentPairInformation.rawJaccard() >= minJaccard) {   // ****** USING RAWJACCARD INSTEAD OF JACCARD
            if(direction == 0) {
                edge.segmentId0 = primarySegmentId;
                edge.segmentId1 = segmentId1;
            } else {
                edge.segmentId0 = segmentId1;
                edge.segmentId1 = primarySegmentId;
                reverse(edge.segmentIds.begin(), edge.segmentIds.end());
            }
            edges.push_back(edge);
            return;
        }

        edge.segmentIds.push_back(segmentId1);
        segmentId0 = segmentId1;
    }
}



// This storesin the Jaccard graph the edges found by all threads.
void JaccardGraph::storeEdges()
{
    JaccardGraph& jaccardGraph = *this;

    for(const auto& threadEdges: threadEdges) {
        for(const JaccardGraphEdgeInfo& info: threadEdges) {

            const uint64_t segmentId0 = info.segmentId0;
            const uint64_t segmentId1 = info.segmentId1;
            const JaccardGraph::vertex_descriptor v0 = vertexTable[segmentId0];
            const JaccardGraph::vertex_descriptor v1 = vertexTable[segmentId1];

            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(v0, v1, jaccardGraph);
            if(not edgeExists) {
                boost::add_edge(v0, v1,
                    JaccardGraphEdge(info.segmentPairInformation, info.direction, info.segmentIds),
                    jaccardGraph);
            } else {
                jaccardGraph[e].wasFoundInDirection[info.direction] = true;
            }
        }
    }
    threadEdges.clear();
}



// A strong vertex is one that is incident to at least one strong edge.
bool JaccardGraph::isStrongVertex(vertex_descriptor v) const
{
    const JaccardGraph& jaccardGraph = *this;

    // Check the out-edges.
    BGL_FORALL_OUTEDGES(v, e, jaccardGraph, JaccardGraph) {
        if(jaccardGraph[e].isStrong()) {
            return true;
        }
    }

    // Check the in-edges.
    BGL_FORALL_INEDGES(v, e, jaccardGraph, JaccardGraph) {
        if(jaccardGraph[e].isStrong()) {
            return true;
        }
    }

    // We did not find any strong edges.
    return false;
}




// Remove all weak vertices.
void JaccardGraph::removeWeakVertices()
{
    JaccardGraph& jaccardGraph = *this;

    // Find the vertices we are going to remove.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        if(not isStrongVertex(v)) {
            verticesToBeRemoved.push_back(v);
        }
    }

    // Remove the vertices we flagged.
    for(const vertex_descriptor v: verticesToBeRemoved) {
        removeVertex(v);
    }

}



// Remove a vertex, making sure to update the vertexTable.
void JaccardGraph::removeVertex(vertex_descriptor v)
{
    JaccardGraph& jaccardGraph = *this;
    const uint64_t segmentId = jaccardGraph[v].segmentId;
    vertexTable[segmentId] = null_vertex();
    clear_vertex(v, jaccardGraph);
    remove_vertex(v, jaccardGraph);
}



void JaccardGraph::writeGraphviz(
    const string& fileName,
    bool includeIsolatedVertices,
    bool writeLabels) const
{
    ofstream file(fileName);
    writeGraphviz(file, includeIsolatedVertices, writeLabels);
}



void JaccardGraph::writeGraphviz(
    ostream& graphOut,
    bool includeIsolatedVertices,
    bool writeLabels) const
{
    const JaccardGraph& jaccardGraph = *this;

    graphOut << "digraph JaccardGraph {" << endl;

    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        if( includeIsolatedVertices or
            in_degree(v, jaccardGraph) or
            out_degree(v, jaccardGraph)) {
            graphOut << jaccardGraph[v].segmentId;
            if(writeLabels) {
                graphOut << " [label=" << jaccardGraph[v].segmentId << "]";
            }
            graphOut << ";\n";
        }
    }

    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        const JaccardGraphEdge& edge = jaccardGraph[e];
        const JaccardGraph::vertex_descriptor v0 = source(e, jaccardGraph);
        const JaccardGraph::vertex_descriptor v1 = target(e, jaccardGraph);
        const uint64_t segmentId0 = jaccardGraph[v0].segmentId;
        const uint64_t segmentId1 = jaccardGraph[v1].segmentId;

        graphOut << segmentId0 << "->" << segmentId1 << "[";

        // Color the edge based on the direction flags.
        if(edge.wasFoundInDirection[0]) {
            if(edge.wasFoundInDirection[1]) {
                // Found in both directions.
                graphOut << " color=black";
            } else {
                // Only found in the forward direction.
                graphOut << " color=red";
            }
        } else {
            if(edge.wasFoundInDirection[1]) {
                // Only found in the backward direction.
                graphOut << " color=green";
            } else {
                SHASTA_ASSERT(0);
            }
        }

        if(writeLabels) {
            graphOut << " label=\"";
            for(const uint64_t segmentId: edge.segmentIds) {
                graphOut << segmentId << "\\n";
            }
            graphOut << "\"";
        }
        graphOut << "];\n";
    }

    graphOut << "}" << endl;

}



// Write edges in csv format.
void JaccardGraph::writeEdgesCsv(const string& fileName) const
{
    ofstream file(fileName);
    writeEdgesCsv(file);
}
void JaccardGraph::writeEdgesCsv(ostream& csv) const
{
    const JaccardGraph& jaccardGraph = *this;

    csv << "SegmentId0,SegmentId1,FoundForward,FoundBackward,SegmentId\n";
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        const JaccardGraphEdge& edge = jaccardGraph[e];
        const JaccardGraph::vertex_descriptor v0 = source(e, jaccardGraph);
        const JaccardGraph::vertex_descriptor v1 = target(e, jaccardGraph);
        const uint64_t segmentId0 = jaccardGraph[v0].segmentId;
        const uint64_t segmentId1 = jaccardGraph[v1].segmentId;

        for(const uint64_t segmentId: edge.segmentIds) {
            csv << segmentId0 << ",";
            csv << segmentId1 << ",";
            csv << int(edge.wasFoundInDirection[0]) << ",";
            csv << int(edge.wasFoundInDirection[1]) << ",";
            csv << segmentId << "\n";
        }
    }
}



// Compute connected component and store the component
// (define as a cluster) that each segment belongs to.
void JaccardGraph::findClusters(
    uint64_t segmentCount,
    MemoryMapped::Vector<uint64_t>& clusterIds)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minClusterSize = 5;

    const JaccardGraph& jaccardGraph = *this;

    // This must be called without removing any vertices.
    SHASTA_ASSERT(num_vertices(jaccardGraph) == segmentCount);

    // Compute connected components.
    vector<uint64_t> rank(segmentCount);
    vector<uint64_t> parent(segmentCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        disjointSets.make_set(segmentId);
    }
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        const JaccardGraph::vertex_descriptor v0 = source(e, jaccardGraph);
        const JaccardGraph::vertex_descriptor v1 = target(e, jaccardGraph);
        const uint64_t segmentId0 = jaccardGraph[v0].segmentId;
        const uint64_t segmentId1 = jaccardGraph[v1].segmentId;
        disjointSets.union_set(segmentId0, segmentId1);
    }

    // Gather the segments in each connected component.
    vector< vector<uint64_t> > components(segmentCount);
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        const uint64_t componentId = disjointSets.find_set(segmentId);
        components[componentId].push_back(segmentId);
    }

    // Sort the components by decreasing size.
    vector< pair<uint64_t, uint64_t> > componentTable; // pair(componentId, componentSize)
    for(uint64_t componentId=0; componentId<segmentCount; componentId++) {
        const uint64_t componentSize = components[componentId].size();
        if(componentSize >= minClusterSize) {
            componentTable.push_back(make_pair(componentId, componentSize));
        }
    }
    sort(componentTable.begin(), componentTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    // Store the cluster ids.
    clusterIds.resize(segmentCount);
    fill(clusterIds.begin(), clusterIds.end(), invalid<uint64_t>);
    for(uint64_t newComponentId=0; newComponentId<componentTable.size(); newComponentId++) {
        const auto& p = componentTable[newComponentId];
        const uint64_t oldComponentId = p.first;
        const uint64_t componentSize = p.second;
        const vector<uint64_t>& component = components[oldComponentId];
        SHASTA_ASSERT(component.size() == componentSize);
        for(const uint64_t segmentId: component) {
            clusterIds[segmentId] = newComponentId;
        }
    }

}



// Construction of the ExpandedJaccardGraph.
// Each vertex of the JaccardGraph generates a vertex in the ExpandedJaccardGraph.
// Each edge of the JaccardGraph generates a linear chain of vertices
// in the ExpandedJaccardGraph.
ExpandedJaccardGraph::ExpandedJaccardGraph(const JaccardGraph& jaccardGraph)
{
    using Graph = ExpandedJaccardGraph;
    Graph& graph = *this;

    // Generate the vertices.
    std::map<JaccardGraph::vertex_descriptor, Graph::vertex_descriptor> vertexMap;
    BGL_FORALL_VERTICES(v, jaccardGraph, JaccardGraph) {
        const Graph::vertex_descriptor u = add_vertex(
            ExpandedJaccardGraphVertex(jaccardGraph[v].segmentId), graph);
        vertexMap.insert(make_pair(v, u));
    }



    // Each edge of the JaccardGraph generates a linear chain of vertices
    // in the ExpandedJaccardGraph.
    BGL_FORALL_EDGES(e, jaccardGraph, JaccardGraph) {
        const JaccardGraph::vertex_descriptor v0 = source(e, jaccardGraph);
        const JaccardGraph::vertex_descriptor v1 = target(e, jaccardGraph);
        const Graph::vertex_descriptor u0 = vertexMap[v0];
        const Graph::vertex_descriptor u1 = vertexMap[v1];
        const vector<uint64_t>& segmentIds = jaccardGraph[e].segmentIds;

        Graph::vertex_descriptor u = u0;
        for(const uint64_t segmentId: segmentIds) {
            const Graph::vertex_descriptor w = add_vertex(
                ExpandedJaccardGraphVertex(segmentId), graph);
            add_edge(u, w, graph);
            u = w;
        }
        add_edge(u, u1, graph);
    }
}



void ExpandedJaccardGraph::writeGraphviz(const string& fileName) const
{
    ofstream s(fileName);
    writeGraphviz(s);
}
void ExpandedJaccardGraph::writeGraphviz(ostream& s) const
{
    using Graph = ExpandedJaccardGraph;
    const Graph& graph = *this;

    s << "digraph ExpandedJaccardGraph {" << endl;

    // We can't use the segment ids to identify vertices
    // because each segment id can appear multiple times.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        s << "\"" << v << "\" [label=" << graph[v].segmentId << "];\n";
    }

    BGL_FORALL_EDGES(e, graph, Graph) {
        const Graph::vertex_descriptor v0 = source(e, graph);
        const Graph::vertex_descriptor v1 = target(e, graph);

        s << "\"" << v0 << "\"->\"" << v1 << "\";\n";
    }

    s << "}" << endl;

}

