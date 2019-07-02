#ifndef SHASTA_STATIC_EXECUTABLE

// Shasta.
#include "LocalMarkerGraph.hpp"
#include "approximateTopologicalSort.hpp"
#include "findMarkerId.hpp"
#include "orderPairs.hpp"
using namespace ::shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



LocalMarkerGraph::LocalMarkerGraph(
    uint32_t k,
    LongBaseSequences& reads,
    const MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& readRepeatCounts,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MemoryMapped::Vector<MarkerGraph::CompressedVertexId>& globalMarkerGraphVertex,
    const ConsensusCaller& consensusCaller
    ) :
    k(k),
    reads(reads),
    readRepeatCounts(readRepeatCounts),
    markers(markers),
    globalMarkerGraphVertex(globalMarkerGraphVertex),
    consensusCaller(consensusCaller)
{

}


// Find out if a vertex with the given MarkerGraph::VertexId exists.
// If it exists, return make_pair(true, v).
// Otherwise, return make_pair(false, null_vertex());
std::pair<bool, LocalMarkerGraph::vertex_descriptor>
    LocalMarkerGraph::findVertex(MarkerGraph::VertexId vertexId) const
{
    const auto it = vertexMap.find(vertexId);
    if(it == vertexMap.end()) {
        return make_pair(false, null_vertex());
    } else {
        const vertex_descriptor v = it->second;
        return make_pair(true, v);
    }
}


// Add a vertex with the given MarkerGraph::VertexId
// and return its vertex descriptor.
// A vertex with this MarkerGraph::VertexId must not exist.
LocalMarkerGraph::vertex_descriptor
    LocalMarkerGraph::addVertex(
    MarkerGraph::VertexId vertexId,
    int distance,
    MemoryAsContainer<MarkerId> vertexMarkers)
{
    // Check that the vertex does not already exist.
    SHASTA_ASSERT(vertexMap.find(vertexId) == vertexMap.end());

    // Add the vertex and store it in the vertex map.
    const vertex_descriptor v = add_vertex(LocalMarkerGraphVertex(vertexId, distance), *this);
    vertexMap.insert(make_pair(vertexId, v));

    // Fill in the marker information for this vertex.
    LocalMarkerGraphVertex& vertex = (*this)[v];
    vertex.markerInfos.reserve(vertexMarkers.size());
    for(const MarkerId markerId: vertexMarkers) {
        LocalMarkerGraphVertex::MarkerInfo markerInfo;
        markerInfo.markerId = markerId;
        tie(markerInfo.orientedReadId, markerInfo.ordinal) =
            findMarkerId(markerId, markers);
        vertex.markerInfos.push_back(markerInfo);
    }

    return v;
}



// Get the KmerId for a vertex.
KmerId LocalMarkerGraph::getKmerId(vertex_descriptor v) const
{
    const LocalMarkerGraphVertex& vertex = (*this)[v];
    SHASTA_ASSERT(!vertex.markerInfos.empty());
    const MarkerId firstMarkerId = vertex.markerInfos.front().markerId;
    const CompressedMarker& firstMarker = markers.begin()[firstMarkerId];
    const KmerId kmerId = firstMarker.kmerId;

    // Sanity check that all markers have the same kmerId.
    // At some point this can be removed.
    for(const auto& markerInfo: vertex.markerInfos){
        const CompressedMarker& marker = markers.begin()[markerInfo.markerId];
        SHASTA_ASSERT(marker.kmerId == kmerId);
    }

    return kmerId;
}



// Get the repeat counts for a MarkerInfo of a vertex.
vector<uint8_t> LocalMarkerGraph::getRepeatCounts(
    const LocalMarkerGraphVertex::MarkerInfo& markerInfo) const
{
    const OrientedReadId orientedReadId = markerInfo.orientedReadId;
    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();
    const CompressedMarker& marker = markers.begin()[markerInfo.markerId];

    const auto& counts = readRepeatCounts[readId];

    vector<uint8_t> v(k);
    for(size_t i=0; i<k; i++) {
        if(strand == 0) {
            v[i] = counts[marker.position + i];
        } else {
            v[i] = counts[counts.size() - 1 - marker.position - i];
        }
    }

    return v;
}



// Fill in the ConsensusInfo's for each vertex.
void LocalMarkerGraph::computeVertexConsensusInfo()
{

    LocalMarkerGraph& graph = *this;
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        computeVertexConsensusInfo(v);
    }
}
void LocalMarkerGraph::computeVertexConsensusInfo( vertex_descriptor v)
{

    // Short-hands for the graph and the vertex.
    LocalMarkerGraph& graph = *this;
    LocalMarkerGraphVertex& vertex = graph[v];

    // Get the marker k-mer of this vertex.
    const KmerId kmerId = graph.getKmerId(v);
    const Kmer kmer(kmerId, k);
    // Resize the consensus info's for the vertex.
    vertex.coverages.resize(k);

    // Loop over all markers of this vertex.
    for(const auto& markerInfo: vertex.markerInfos) {

        // Get the repeat counts for this marker.
        const vector<uint8_t> counts = graph.getRepeatCounts(markerInfo);
        SHASTA_ASSERT(counts.size() == k);

        // Increment coverage.
        for(size_t position=0; position<k; position++) {
            vertex.coverages[position].addRead(
                AlignedBase(kmer[position]),
                markerInfo.orientedReadId.getStrand(),
                counts[position]);
        }
    }
}



// Store sequence information in the edge.
// This version takes as input a vector of the
// LocalMarkerGraphEdge::Info that caused the edge to be created.
void LocalMarkerGraph::storeEdgeInfo(
    edge_descriptor e,
    const vector<MarkerInterval>& intervals)
{
    LocalMarkerGraph& graph = *this;
    LocalMarkerGraphEdge& edge = graph[e];

    // Map to store the oriented read ids and ordinals, grouped by sequence.
    std::map<LocalMarkerGraphEdge::Sequence, vector<MarkerIntervalWithRepeatCounts> > sequenceTable;
    for(const MarkerInterval& interval: intervals) {
        const CompressedMarker& marker0 = markers.begin(interval.orientedReadId.getValue())[interval.ordinals[0]];
        const CompressedMarker& marker1 = markers.begin(interval.orientedReadId.getValue())[interval.ordinals[1]];

        // Fill in the sequence information and, if necessary, the base repeat counts.
        LocalMarkerGraphEdge::Sequence sequence;
        MarkerIntervalWithRepeatCounts intervalWithRepeatCounts(interval);
        if(marker1.position <= marker0.position + k) {
            sequence.overlappingBaseCount = uint8_t(marker0.position + k - marker1.position);
            const auto& repeatCounts = readRepeatCounts[interval.orientedReadId.getReadId()];
            for(uint32_t i=0; i<sequence.overlappingBaseCount; i++) {
                uint32_t position = marker1.position + i;
                uint8_t repeatCount = 0;
                if(interval.orientedReadId.getStrand() == 0) {
                    repeatCount = repeatCounts[position];
                } else {
                    repeatCount = repeatCounts[repeatCounts.size() - 1 - position];
                }
                intervalWithRepeatCounts.repeatCounts.push_back(repeatCount);
            }
        } else {
            sequence.overlappingBaseCount = 0;
            const auto read = reads[interval.orientedReadId.getReadId()];
            const uint32_t readLength = uint32_t(read.baseCount);
            for(uint32_t position=marker0.position+k;  position!=marker1.position; position++) {
                Base base;
                if(interval.orientedReadId.getStrand() == 0) {
                    base = read.get(position);
                } else {
                    base = read.get(readLength - 1 - position);
                    base.complementInPlace();
                }
                sequence.sequence.push_back(base);
            }
            const auto repeatCounts = readRepeatCounts[interval.orientedReadId.getReadId()];
            for(uint32_t position=marker0.position+k;  position!=marker1.position; position++) {
                uint8_t repeatCount;
                if(interval.orientedReadId.getStrand() == 0) {
                    repeatCount = repeatCounts[position];
                } else {
                    repeatCount = repeatCounts[readLength - 1 - position];
                }
                intervalWithRepeatCounts.repeatCounts.push_back(repeatCount);
            }

        }

        // Store it.
        sequenceTable[sequence].push_back(intervalWithRepeatCounts);

    }

    // Copy to the edge infos data structure.
    edge.infos.clear();
    copy(sequenceTable.begin(), sequenceTable.end(), back_inserter(edge.infos));

    // Sort by decreasing size of the infos vector.
    sort(edge.infos.begin(), edge.infos.end(),
        OrderPairsBySizeOfSecondGreater<
        LocalMarkerGraphEdge::Sequence,
        vector<MarkerIntervalWithRepeatCounts> >());

}



// Look for the ordinal for a given oriented read id.
// If found, returns pair(true, ordinal).
// Otherwise, returns pair(false, don't care).
// If more than an ordinal is found, the first one is returned.
pair<bool, uint32_t> LocalMarkerGraphVertex::getOrdinal(
    OrientedReadId orientedReadId) const
{
    for(const MarkerInfo& markerInfo: markerInfos) {
        if(markerInfo.orientedReadId == orientedReadId) {
            return make_pair(true, markerInfo.ordinal);
        }
    }
    return make_pair(false, std::numeric_limits<uint32_t>::max());
}



// Look for the ordinals for a given oriented read id.
// If found, returns true.
// If more than an ordinal pairs is found, the first one is returned.
bool LocalMarkerGraphEdge::getOrdinals(
    OrientedReadId orientedReadId,
    array<uint32_t, 2>& ordinals) const
{
    for(const pair<Sequence, vector<MarkerIntervalWithRepeatCounts> >& p: infos) {
        for(const MarkerIntervalWithRepeatCounts& interval: p.second) {
            if(interval.orientedReadId == orientedReadId) {
                ordinals = interval.ordinals;
                return true;
            }
        }
    }

    // If getting here, we did not find it.
    return false;
}



// Approximate topological sort, adding edges
// in order of decreasing coverage. The topological sort
// stored in LocalMarkerGrapg2Vertex::rank.
void LocalMarkerGraph::approximateTopologicalSort()
{
    LocalMarkerGraph& graph = *this;

    vector<pair<uint32_t, edge_descriptor> > edgeTable;
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph) {
        edgeTable.push_back(make_pair(graph[e].coverage(), e));
    }
    sort(edgeTable.begin(), edgeTable.end(),
        std::greater< pair<uint32_t, edge_descriptor> >());

    vector<edge_descriptor> sortedEdges;
    for(const auto& p: edgeTable) {
        sortedEdges.push_back(p.second);
    }

    shasta::approximateTopologicalSort(graph, sortedEdges);


    // Also store the vertices in topological sort order.
    vector< pair<size_t, vertex_descriptor> > vertexTable;
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        vertexTable.push_back(make_pair(graph[v].rank, v));
    }
    sort(vertexTable.begin(), vertexTable.end());
    topologicallySortedVertices.clear();
    for(const auto& p: vertexTable) {
        topologicallySortedVertices.push_back(p.second);
    }

}

#endif
