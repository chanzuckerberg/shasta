
// Shasta
#include "mode3.hpp"
#include "findMarkerId.hpp"
#include "MarkerGraph.hpp"
#include "ReadFlags.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <map>
#include <queue>



DynamicAssemblyGraph::DynamicAssemblyGraph(
    const MemoryMapped::Vector<ReadFlags>& readFlags,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize,
    size_t threadCount) :
    MultithreadedObject<DynamicAssemblyGraph>(*this),
    readFlags(readFlags),
    markerGraph(markerGraph),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize),
    threadCount(threadCount)
{
    createVertices(markers);
    computeMarkerGraphEdgeTable();

    computePseudoPaths();
    writePseudoPaths("PseudoPaths.csv");

    createEdges();
}



DynamicAssemblyGraph::~DynamicAssemblyGraph()
{
    if(markerGraphEdgeTable.isOpen) {
        markerGraphEdgeTable.remove();
    }
    if(pseudoPaths.isOpen()) {
        pseudoPaths.remove();
    }
}



// Each  linear chain of marker graph edges generates a vertex
// of the DynamicAssemblyGraph.
void DynamicAssemblyGraph::createVertices(
     const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers)
{
    const MarkerGraph::EdgeId edgeCount = markerGraph.edges.size();
    vector<bool> wasFound(edgeCount, false);

    using MarkerGraphPath = vector<MarkerGraph::EdgeId>;
    MarkerGraphPath nextEdges;
    MarkerGraphPath previousEdges;
    MarkerGraphPath path;
    MarkerGraphPath reverseComplementedPath;

    // Main loop over all edges of the marker graph.
    // At each iteration we find a new linear path of edges.
    for(MarkerGraph::EdgeId startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {

        // If we already found this edge, skip it.
        // It is part of a path we already found.
        if(wasFound[startEdgeId]) {
            continue;
        }

        // Follow the path forward.
        nextEdges.clear();
        MarkerGraph::EdgeId edgeId = startEdgeId;
        bool isCircular = false;
        while(true) {
            const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
            const MarkerGraph::VertexId v1 = edge.target;
            const auto outEdges = markerGraph.edgesBySource[v1];
            if(outEdges.size() != 1) {
                break;
            }
            const auto inEdges = markerGraph.edgesByTarget[v1];
            if(inEdges.size() != 1) {
                break;
            }
            edgeId = outEdges[0];
            if(edgeId == startEdgeId) {
                isCircular = true;
                break;
            }
            nextEdges.push_back(edgeId);
            SHASTA_ASSERT(not wasFound[edgeId]);
        }

        // Follow the path backward.
        previousEdges.clear();
        if(!isCircular) {
            edgeId = startEdgeId;
            while(true) {
                const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
                const MarkerGraph::VertexId v0 = edge.source;
                const auto outEdges = markerGraph.edgesBySource[v0];
                if(outEdges.size() != 1) {
                    break;
                }
                const auto inEdges = markerGraph.edgesByTarget[v0];
                if(inEdges.size() != 1) {
                    break;
                }
                edgeId = inEdges[0];
                previousEdges.push_back(edgeId);
                SHASTA_ASSERT(not wasFound[edgeId]);
            }
        }

        // Gather the path.
        path.clear();
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter(path));
        path.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter(path));

        // Mark all the edges in the path as found.
        for(const MarkerGraph::EdgeId edgeId: path) {
            if(wasFound[edgeId]) {
                cout << "Assertion failed at " << edgeId << endl;
                SHASTA_ASSERT(0);
            }
            wasFound[edgeId] = true;
        }



        // Check if this path originated from one of the read graph
        // components that we want to assemble
        // (one of each reverse complemented pair).
        // If not, don't store it.
        // This way we do a single-stranded assembly.
        {
            const MarkerGraph::Edge& edge = markerGraph.edges[path.front()];
            const MarkerGraphVertexId vertexId = edge.source;
            const span<const MarkerId> markerIds = markerGraph.getVertexMarkerIds(vertexId);
            const MarkerId markerId = markerIds[0];
            const OrientedReadId orientedReadId = findMarkerId(markerId, markers).first;
            const ReadId readId = orientedReadId.getReadId();
            const Strand strand = orientedReadId.getStrand();
            if(readFlags[readId].strand != strand) {
                continue;
            }
        }

        // This path generates a new vertex of the assembly graph.
        boost::add_vertex(DynamicAssemblyGraphVertex(path), *this);
    }



    // Check that all edges of the marker graph were found.
    SHASTA_ASSERT(find(wasFound.begin(), wasFound.end(), false) == wasFound.end());

}



MarkerGraphEdgeInfo::MarkerGraphEdgeInfo(
    MarkerGraph::EdgeId edgeIdArgument, bool isVirtualArgument)
{
    isVirtual = uint64_t(isVirtualArgument & 1);
    edgeId = edgeIdArgument & 0x7fffffffffffffffULL;
}



DynamicAssemblyGraphVertex::DynamicAssemblyGraphVertex(
    const vector<MarkerGraph::EdgeId>& pathArgument)
{
    for(const MarkerGraph::EdgeId edgeId: pathArgument) {
        path.push_back(MarkerGraphEdgeInfo(edgeId, false));
    }
}



// For each marker graph edge, store in the marker graph edge table
// the DynamicAssemblyGraph vertex (segment)
// and position, if any.
// This is needed when computing pseudopaths.
void DynamicAssemblyGraph::computeMarkerGraphEdgeTable()
{
    G& g = *this;

    // Initialize the marker graph edge table.
    markerGraphEdgeTable.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "tmp-mode3-MarkerGraphEdgeTable"),
        largeDataPageSize);
    markerGraphEdgeTable.resize(markerGraph.edges.size());
    fill(markerGraphEdgeTable.begin(), markerGraphEdgeTable.end(), make_pair(null_vertex(), 0));

    // Store a vector of all vertices.
    markerGraphEdgeTableData.allVertices.clear();
    BGL_FORALL_VERTICES(v, g, DynamicAssemblyGraph) {
        markerGraphEdgeTableData.allVertices.push_back(v);
    }

    // Fill in the marker graph edge table.
    const uint64_t batchSize = 100;
    setupLoadBalancing(markerGraphEdgeTableData.allVertices.size(), batchSize);
    runThreads(&G::computeMarkerGraphEdgeTableThreadFunction, threadCount);

    // Clean up.
    markerGraphEdgeTableData.allVertices.clear();
    markerGraphEdgeTableData.allVertices.shrink_to_fit();
}



void DynamicAssemblyGraph::computeMarkerGraphEdgeTableThreadFunction(size_t threadId)
{
    G& g = *this;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all vertices assigned to this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            const vertex_descriptor v = markerGraphEdgeTableData.allVertices[i];
            const vector<MarkerGraphEdgeInfo>& path = g[v].path;

            // Loop over the path of this vertex (segment).
            for(uint64_t position=0; position<path.size(); position++) {
                const MarkerGraphEdgeInfo& info = path[position];

                // Skip virtual edges.
                if(info.isVirtual) {
                    continue;
                }

                // Store the marker graph edge table entry for this edge.
                const MarkerGraph::EdgeId edgeId = info.edgeId;
                SHASTA_ASSERT(edgeId < markerGraphEdgeTable.size());
                markerGraphEdgeTable[edgeId] = make_pair(v, position);
            }
        }

    }
}




void DynamicAssemblyGraph::computePseudoPaths()
{
    pseudoPaths.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "tmp-mode3-PseudoPaths"),
        largeDataPageSize);

    uint64_t batchSize = 1000;
    pseudoPaths.beginPass1(2 * readFlags.size());
    setupLoadBalancing(markerGraphEdgeTable.size(), batchSize);
    runThreads(&G::computePseudoPathsPass1, threadCount);
    pseudoPaths.beginPass2();
    setupLoadBalancing(markerGraphEdgeTable.size(), batchSize);
    runThreads(&G::computePseudoPathsPass2, threadCount);
    pseudoPaths.endPass2();

    batchSize = 100;
    setupLoadBalancing(pseudoPaths.size(), batchSize);
    runThreads(&G::sortPseudoPaths, threadCount);
}



void DynamicAssemblyGraph::computePseudoPathsPass1(size_t threadId)
{
    computePseudoPathsPass12(1);
}



void DynamicAssemblyGraph::computePseudoPathsPass2(size_t threadId)
{
    computePseudoPathsPass12(2);
}



void DynamicAssemblyGraph::computePseudoPathsPass12(uint64_t pass)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker graph edges assigned to this batch.
        for(MarkerGraph::EdgeId edgeId=begin; edgeId!=end; ++edgeId) {
            const auto& p = markerGraphEdgeTable[edgeId];

            // Get the DynamicAssemblyGraph vertex and position, if any.
            const vertex_descriptor v = p.first;
            if(v == null_vertex()) {
                continue;
            }
            const uint32_t position = p.second;

            // Loop over the marker intervals of this read.
            const auto markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
            for(const MarkerInterval& markerInterval: markerIntervals) {
                const OrientedReadId orientedReadId = markerInterval.orientedReadId;

                if(pass == 1) {
                    pseudoPaths.incrementCountMultithreaded(orientedReadId.getValue());
                } else {
                    PseudoPathEntry pseudoPathEntry;
                    pseudoPathEntry.v = v;
                    pseudoPathEntry.position = position;
                    pseudoPathEntry.ordinals = markerInterval.ordinals;
                    pseudoPaths.storeMultithreaded(orientedReadId.getValue(), pseudoPathEntry);
                }
            }
        }
    }
}



void DynamicAssemblyGraph::sortPseudoPaths(size_t threadId)
{
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker graph edges assigned to this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            auto pseudoPath = pseudoPaths[i];
            sort(pseudoPath.begin(), pseudoPath.end());
        }
    }
}




void DynamicAssemblyGraph::writePseudoPaths(const string& fileName) const
{
    ofstream csv(fileName);
    writePseudoPaths(csv);
}



void DynamicAssemblyGraph::writePseudoPaths(ostream& csv) const
{
    csv << "OrientedReadId,Vertex,Position,Ordinal0,Ordinal1\n";

    for(ReadId readId=0; readId<readFlags.size(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto pseudoPath = pseudoPaths[orientedReadId.getValue()];

            for(const auto& pseudoPathEntry: pseudoPath) {
                csv << orientedReadId << ",";
                csv << pseudoPathEntry.v << ",";
                csv << pseudoPathEntry.position << ",";
                csv << pseudoPathEntry.ordinals[0] << ",";
                csv << pseudoPathEntry.ordinals[1] << "\n";
            }
        }
    }
}



void DynamicAssemblyGraph::createEdges()
{
    G& g = *this;

    // Gather transitions, and store them keyed by the
    // pair of vertices involved.
    using Key = pair<vertex_descriptor, vertex_descriptor>;
    using Value = vector< pair<OrientedReadId, Transition> >;
    std::map< Key, Value> m;


    for(ReadId readId=0; readId<readFlags.size(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto pseudoPath = pseudoPaths[orientedReadId.getValue()];

            if(pseudoPath.size() < 2) {
                continue;
            }

            for(uint64_t i=1; i<pseudoPath.size(); i++) {
                const auto& previous = pseudoPath[i-1];
                const auto& current = pseudoPath[i];
                if(previous.v == current.v) {
                    continue;
                }

                const Key key = make_pair(previous.v, current.v);
                m[key].push_back(
                    make_pair(orientedReadId, Transition({previous, current})));

            }
        }
    }


    ofstream csv("Transitions.csv");
    for(const auto& p: m) {
        const vertex_descriptor v0 = p.first.first;
        const vertex_descriptor v1 = p.first.second;
        for(const auto& q: p.second) {
            const OrientedReadId orientedReadId = q.first;
            // const Transition& transition = q.second.second;
            csv << v0 << ",";
            csv << v1 << ",";
            csv << orientedReadId << "\n";
        }
    }



    ofstream dot("Transitions.dot");
    const uint64_t minCoverage = 1;
    dot << "digraph G {\n";
    BGL_FORALL_VERTICES(v, g, G) {
        dot << "\"" << v << "\";\n";
    }
    for(const auto& p: m) {
        const vertex_descriptor v0 = p.first.first;
        const vertex_descriptor v1 = p.first.second;
        const uint64_t coverage = p.second.size();
        if(coverage >= minCoverage) {
            dot << "\"" << v0 << "\"->\"";
            dot << v1 << "\" [penwidth=" << 0.2*double(coverage) << "];\n";
        }
    }
    dot << "}\n";

    for(const auto& p: m) {
        const vertex_descriptor v0 = p.first.first;
        const vertex_descriptor v1 = p.first.second;
        const uint64_t coverage = p.second.size();
        if(coverage >= minCoverage) {
            add_edge(v0, v1, DynamicAssemblyGraphEdge(coverage), g);
        }
    }

}



// The mode3::AssemblyGraph stores a permanent and immutable copy
// of the DynamicAssemblyGraph in memory mapped data structures,
// so it can be accessed inthe http server and in the Python API.
mode3::AssemblyGraph::AssemblyGraph(
    const DynamicAssemblyGraph& dynamicAssemblyGraph,
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize) :
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize)
{

    // Create a gfa segment for each vertex of the DynamicAssemblyGraph.
    std::map<DynamicAssemblyGraph::vertex_descriptor, uint64_t> m;
    uint64_t segmentId = 0;
    paths.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "Mode3-Paths"),
        largeDataPageSize);
    BGL_FORALL_VERTICES(v, dynamicAssemblyGraph, DynamicAssemblyGraph) {
        paths.appendVector(dynamicAssemblyGraph[v].path);
        m.insert(make_pair(v, segmentId++));
    }

    // Create a gfa Link for each edge of the DynamicAssemblyGraph.
    links.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "Mode3-Links"),
        largeDataPageSize);
    BGL_FORALL_EDGES(e, dynamicAssemblyGraph, DynamicAssemblyGraph) {
        const DynamicAssemblyGraph::vertex_descriptor v0 = source(e, dynamicAssemblyGraph);
        const DynamicAssemblyGraph::vertex_descriptor v1 = target(e, dynamicAssemblyGraph);
        const uint64_t segmentId0 = m[v0];
        const uint64_t segmentId1 = m[v1];
        links.push_back(Link(segmentId0, segmentId1, dynamicAssemblyGraph[e].coverage));
    }
    createConnectivity();

    cout << "The mode 3 assembly graph has " << paths.size() << " segments and " <<
        links.size() << " links." << endl;
}



// Constructor from binary data.
mode3::AssemblyGraph::AssemblyGraph(const string& largeDataFileNamePrefix) :
    largeDataFileNamePrefix(largeDataFileNamePrefix)
{
    paths.accessExistingReadOnly(largeDataFileNamePrefix + "Mode3-Paths");
    links.accessExistingReadOnly(largeDataFileNamePrefix + "Mode3-Links");
    linksBySource.accessExistingReadOnly(largeDataFileNamePrefix + "Mode3-LinksBySource");
    linksByTarget.accessExistingReadOnly(largeDataFileNamePrefix + "Mode3-LinksByTarget");
}



void mode3::AssemblyGraph::createConnectivity()
{
    linksBySource.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "Mode3-LinksBySource"),
        largeDataPageSize);
    linksByTarget.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "Mode3-LinksByTarget"),
        largeDataPageSize);

    linksBySource.beginPass1(links.size());
    linksByTarget.beginPass1(links.size());
    for(uint64_t linkId=0; linkId<links.size(); linkId++) {
        const Link& link = links[linkId];
        linksBySource.incrementCount(link.segmentId0);
        linksByTarget.incrementCount(link.segmentId1);
    }
    linksBySource.beginPass2();
    linksByTarget.beginPass2();
    for(uint64_t linkId=0; linkId<links.size(); linkId++) {
        const Link& link = links[linkId];
        linksBySource.store(link.segmentId0, linkId);
        linksByTarget.store(link.segmentId1, linkId);
    }
    linksBySource.endPass2();
    linksByTarget.endPass2();
}



// Create the LocalAssemblyGraph using a BFS
// that starts at the specified vertex and moves away
// (in both directions) up to the specified distance
mode3::LocalAssemblyGraph::LocalAssemblyGraph(
    const AssemblyGraph& assemblyGraph,
    uint64_t startSegmentId,
    uint64_t maxDistance)
{
    LocalAssemblyGraph& localAssemblyGraph= *this;

    // The BFS queue.
    std::queue<uint64_t> q;

    // Map segments in the AssemblyGraph to vertices in
    // the LocalAssemblyGraph.
    std::map<uint64_t, vertex_descriptor> segmentMap;

    // Initialize the BFS.
    q.push(startSegmentId);
    const vertex_descriptor vStart = addVertex(startSegmentId, 0, assemblyGraph.paths[startSegmentId]);
    segmentMap.insert(make_pair(startSegmentId, vStart));



    // BFS.
    while(not q.empty()) {

        // Dequeue a segment.
        const uint64_t segmentId0 = q.front();
        q.pop();
        const vertex_descriptor v0 = segmentMap[segmentId0];
        const uint64_t distance0 = localAssemblyGraph[v0].distance;
        const uint64_t distance1 = distance0 + 1;

        // Loop over children.
        for(const uint64_t linkId: assemblyGraph.linksBySource[segmentId0]) {
            const Link& link = assemblyGraph.links[linkId];
            const uint64_t segmentId1 = link.segmentId1;
            if(segmentMap.find(segmentId1) != segmentMap.end()) {
                // We already encountered this segment.
                continue;
            }
            const vertex_descriptor v1 = addVertex(segmentId1, distance1, assemblyGraph.paths[segmentId1]);
            segmentMap.insert(make_pair(segmentId1, v1));
            if(distance1 < maxDistance) {
                q.push(segmentId1);
            }
        }

        // Loop over parents.
        for(const uint64_t linkId: assemblyGraph.linksByTarget[segmentId0]) {
            const Link& link = assemblyGraph.links[linkId];
            const uint64_t segmentId1 = link.segmentId0;
            if(segmentMap.find(segmentId1) != segmentMap.end()) {
                // We already encountered this segment.
                continue;
            }
            const vertex_descriptor v1 = addVertex(segmentId1, distance1, assemblyGraph.paths[segmentId1]);
            segmentMap.insert(make_pair(segmentId1, v1));
            if(distance1 < maxDistance) {
                q.push(segmentId1);
            }
        }
    }



    // Add the edges.
    for(const auto& p: segmentMap) {
        const uint64_t segmentId0 = p.first;
        const vertex_descriptor v0 = p.second;

        for(const uint64_t linkId: assemblyGraph.linksBySource[segmentId0]) {
            const Link& link = assemblyGraph.links[linkId];
            const uint64_t segmentId1 = link.segmentId1;
            const auto it1 = segmentMap.find(segmentId1);
            if(it1 == segmentMap.end()) {
                continue;
            }
            const vertex_descriptor v1 = it1->second;
            boost::add_edge(v0, v1, LocalAssemblyGraphEdge(link.coverage), localAssemblyGraph);
        }
    }
}



mode3::LocalAssemblyGraphVertex::LocalAssemblyGraphVertex(
    uint64_t segmentId,
    uint64_t distance,
    const span<const MarkerGraphEdgeInfo> pathArgument) :
    segmentId(segmentId),
    distance(distance)
{
    copy(pathArgument.begin(), pathArgument.end(), back_inserter(path));
}



mode3::LocalAssemblyGraphVertex::LocalAssemblyGraphVertex() :
    segmentId(0),
    distance(0)
{
}



mode3::LocalAssemblyGraph::vertex_descriptor mode3::LocalAssemblyGraph::addVertex(
    uint64_t segmentId,
    uint64_t distance,
    const span<const MarkerGraphEdgeInfo> path)
{
    return add_vertex(LocalAssemblyGraphVertex(segmentId, distance, path), *this);
}




void mode3::LocalAssemblyGraph::writeGraphviz(const string& fileName) const
{
    ofstream s(fileName);
    writeGraphviz(s);
}



void mode3::LocalAssemblyGraph::writeGraphviz(ostream& s) const
{
    const LocalAssemblyGraph localAssemblyGraph = *this;

    s << "digraph LocalAssemblyGraph {\n";

    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        s << localAssemblyGraph[v].segmentId << ";\n";
    }

    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const vertex_descriptor v0 = source(e, localAssemblyGraph);
        const vertex_descriptor v1 = target(e, localAssemblyGraph);
        s << localAssemblyGraph[v0].segmentId << "->";
        s << localAssemblyGraph[v1].segmentId <<
            "[penwidth=" << 0.2 * double(localAssemblyGraph[e].coverage) << "];\n";
    }

    s << "}\n";
}
