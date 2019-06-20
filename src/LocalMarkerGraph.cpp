#ifndef SHASTA_STATIC_EXECUTABLE

// Shasta
#include "LocalMarkerGraph.hpp"
#include "approximateTopologicalSort.hpp"
#include "findMarkerId.hpp"
#include "LongBaseSequence.hpp"
// #include "Histogram.hpp"
// #include "Marker.hpp"
#include "orderPairs.hpp"
// #include "ReadId.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// SeqAn.
#include <seqan/graph_msa.h>
#include <seqan/version.h>

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/topological_sort.hpp>

// Standard libraries.
#include "iterator.hpp"



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
    CZI_ASSERT(vertexMap.find(vertexId) == vertexMap.end());

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
    CZI_ASSERT(!vertex.markerInfos.empty());
    const MarkerId firstMarkerId = vertex.markerInfos.front().markerId;
    const CompressedMarker& firstMarker = markers.begin()[firstMarkerId];
    const KmerId kmerId = firstMarker.kmerId;

    // Sanity check that all markers have the same kmerId.
    // At some point this can be removed.
    for(const auto& markerInfo: vertex.markerInfos){
        const CompressedMarker& marker = markers.begin()[markerInfo.markerId];
        CZI_ASSERT(marker.kmerId == kmerId);
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
        CZI_ASSERT(counts.size() == k);

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



void LocalMarkerGraphEdge::computeCoverage()
{
    // The SeqAn alignment must have been computed.
    CZI_ASSERT(seqanAlignmentWasComputed);

    // The length of the alignment.
    // This includes gaps.
    const size_t n = seqan::length(seqan::row(seqanAlignment, 0));

    // The number of reads in the alignment.
    const size_t m = alignmentInfos.size();

    // Loop over all positions of the alignment.
    vector<size_t> positions(m, 0);
    coverages.resize(n);
    for(size_t i=0; i<n; i++) {
        Coverage& coverage =coverages[i];

        // Loop over all reads in the alignment to compute coverage
        // for each base and repeat count.
        for(size_t j=0; j<m; j++) {
            if(seqan::isGap(seqan::row(seqanAlignment, j), i)) {
                coverage.addRead(
                    AlignedBase::gap(),
                    alignmentInfos[j].orientedReadId.getStrand(),
                    0);
            } else {

                // Extract the read base and repeat count at this position
                // in the alignment.
                const Base base = alignmentInfos[j].sequence[positions[j]];
                const size_t repeatCount = alignmentInfos[j].repeatCounts[positions[j]];
                ++positions[j];

                // Increment coverage for this base and repeat count.
                coverage.addRead(
                    AlignedBase(base),
                    alignmentInfos[j].orientedReadId.getStrand(),
                    repeatCount);
            }
        }
    }
}



// Compute the set of vertices that corresponds to a given oriented read.
// Vertices are returned in a pair with the corresponding ordinal,
// sorted by the ordinal.
void LocalMarkerGraph::getOrientedReadVertices(
    OrientedReadId orientedReadId,
    vector< pair<uint32_t, vertex_descriptor> >& orientedReadVertices) const
{
    const LocalMarkerGraph& graph = *this;

    orientedReadVertices.clear();
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        bool found = true;
        uint32_t ordinal;
        tie(found, ordinal) = graph[v].getOrdinal(orientedReadId);
        if(found) {
            orientedReadVertices.push_back(make_pair(ordinal, v));
        }
    }
    sort(orientedReadVertices.begin(), orientedReadVertices.end());

}



// Compute the set of edges that corresponds to a given oriented read.
// Each edge is returned in a tuple containing the two ordinals
// for the given oriented read.
// The edges are computed sorted by the ordinals.
void LocalMarkerGraph::getOrientedReadEdges(
    OrientedReadId orientedReadId,
    vector< pair< array<uint32_t, 2>, edge_descriptor> >& orientedReadEdges) const
{
    const LocalMarkerGraph& graph = *this;

    orientedReadEdges.clear();
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph) {
        array<uint32_t, 2> ordinals;
        if(graph[e].getOrdinals(orientedReadId, ordinals)) {
            // const vertex_descriptor v0 = source(e, graph);
            // const vertex_descriptor v1 = target(e, graph);
            // const LocalMarkerGraphVertex& vertex0 = graph[v0];
            // const LocalMarkerGraphVertex& vertex1 = graph[v1];
            orientedReadEdges.push_back(make_pair(ordinals, e));
        }
    }
    sort(orientedReadEdges.begin(), orientedReadEdges.end());

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

#endif
