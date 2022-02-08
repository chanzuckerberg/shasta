
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
    size_t threadCount) :
    MultithreadedObject<DynamicAssemblyGraph>(*this),
    markerGraph(markerGraph)
{
    createVertices(readFlags, markers);
    computeOrdinalRanges(threadCount);
}



// Each  linear chain of marker graph edges generates a vertex
// of the DynamicAssemblyGraph.
void DynamicAssemblyGraph::createVertices(
    const MemoryMapped::Vector<ReadFlags>& readFlags,
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

    cout << "The initial assembly graph has " << num_vertices(*this) <<
        " vertices." << endl;

}



DynamicAssemblyGraphVertex::MarkerGraphEdgeInfo::MarkerGraphEdgeInfo(
    MarkerGraph::EdgeId edgeIdArgument, bool isVirtualArgument)
{
    isVirtual = uint64_t(isVirtualArgument & 1);
    edgeId = edgeIdArgument & 0x7fffffffffffffffULL;

}



DynamicAssemblyGraphVertex::DynamicAssemblyGraphVertex(const vector<MarkerGraph::EdgeId>& path)
{
    for(const MarkerGraphEdgeId edgeId: path) {
        markerGraphEdges.push_back(MarkerGraphEdgeInfo(edgeId, false));
    }
}



// Compute ordinal ranges for all vertices in the graph.
void DynamicAssemblyGraph::computeOrdinalRanges(size_t threadCount)
{
    computeOrdinalRangesData.allVertices.clear();
    BGL_FORALL_VERTICES(v, *this, DynamicAssemblyGraph) {
        computeOrdinalRangesData.allVertices.push_back(v);
    }

    const uint64_t batchSize = 100;
    setupLoadBalancing(computeOrdinalRangesData.allVertices.size(), batchSize);
    runThreads(&DynamicAssemblyGraph::computeOrdinalRangesThreadFunction, threadCount);
}



void DynamicAssemblyGraph::computeOrdinalRangesThreadFunction(size_t threadId)
{
    auto& g = *this;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices assigned to this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            const vertex_descriptor v = computeOrdinalRangesData.allVertices[i];

            g[v].computeOrdinalRanges(markerGraph);
        }
    }
}



// Compute ordinal ranges for a vertex.
void DynamicAssemblyGraphVertex::computeOrdinalRanges(const MarkerGraph& markerGraph)
{
    std::map<OrientedReadId, pair<uint32_t, uint32_t> > m;

    // Loop over non-virtual marker graph edges.
    for(const MarkerGraphEdgeInfo& info: markerGraphEdges) {
        if(info.isVirtual) {
            continue;
        }

        // Loop over marker intervals of this edge.
        const auto& markerIntervals = markerGraph.edgeMarkerIntervals[info.edgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            const uint32_t ordinal0 = markerInterval.ordinals[0];
            const uint32_t ordinal1 = markerInterval.ordinals[1];
            SHASTA_ASSERT(ordinal0 < ordinal1);

            // Update the map for this read.
            auto it = m.find(orientedReadId);
            if(it == m.end()) {
                m.insert(make_pair(orientedReadId, make_pair(ordinal0, ordinal1)));
            } else {
                auto& p = it->second;
                p.first = min(p.first, ordinal0);
                p.second = max(p.second, ordinal1);
            }

        }
    }

    // Store in the vertex what we found.
    ordinalRanges.clear();
    for(const auto& p: m) {
        OrdinalRange ordinalRange;
        ordinalRange.orientedReadId = p.first;
        ordinalRange.minOrdinal = p.second.first;
        ordinalRange.maxOrdinal = p.second.second;
        ordinalRanges.push_back(ordinalRange);
    }
}
