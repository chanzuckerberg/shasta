
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
        boost::add_vertex(DynamicAssemblyGraphVertex(path, nextVertexId++), *this);
    }



    // Check that all edges of the marker graph were found.
    SHASTA_ASSERT(find(wasFound.begin(), wasFound.end(), false) == wasFound.end());

    cout << "The initial assembly graph has " << num_vertices(*this) <<
        " vertices." << endl;

}



MarkerGraphEdgeInfo::MarkerGraphEdgeInfo(
    MarkerGraph::EdgeId edgeIdArgument, bool isVirtualArgument)
{
    isVirtual = uint64_t(isVirtualArgument & 1);
    edgeId = edgeIdArgument & 0x7fffffffffffffffULL;
}



DynamicAssemblyGraphVertex::DynamicAssemblyGraphVertex(
    const vector<MarkerGraph::EdgeId>& pathArgument,
    uint64_t vertexId) :
    vertexId(vertexId)
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
    const G& g = *this;

    csv << "OrientedReadId,VertexId,Position,Ordinal0,Ordinal1\n";

    for(ReadId readId=0; readId<readFlags.size(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto pseudoPath = pseudoPaths[orientedReadId.getValue()];

            for(const auto& pseudoPathEntry: pseudoPath) {
                csv << orientedReadId << ",";
                csv << g[pseudoPathEntry.v].vertexId << ",";
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
            csv << g[v0].vertexId << ",";
            csv << g[v1].vertexId << ",";
            csv << orientedReadId << "\n";
        }
    }



    ofstream dot("Transitions.dot");
    const uint64_t minCoverage = 1;
    dot << "digraph G {\n";
    for(const auto& p: m) {
        const vertex_descriptor v0 = p.first.first;
        const vertex_descriptor v1 = p.first.second;
        const uint64_t coverage = p.second.size();
        if(coverage >= minCoverage) {
            dot << g[v0].vertexId << "->";
            dot << g[v1].vertexId << " [penwidth=" << 0.2*double(coverage) << "];\n";
        }
    }
    dot << "}\n";
}

