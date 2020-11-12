#ifndef SHASTA_ALIGNMENT_GRAPH_HPP
#define SHASTA_ALIGNMENT_GRAPH_HPP

/*******************************************************************************

Class AlignmentGraph is used to compute an alignment of the markers
of two oriented reads.

In the alignment graph, each vertex corresponds to a pair of markers
in the two oriented reads that have the same k-mer.

Edges in the alignment  graph are created between vertices
that are sufficiently close on both oriented reads.
Each edge is assigned a weight based on the distances
between the two pairs of markers on each of the
oriented reads, measured using ordinals (not position).

To find a good alignment, we find a shortest path in the graph.

*******************************************************************************/

// Shasta
#include "CompactUndirectedGraph.hpp"
#include "Marker.hpp"
#include "shortestPath.hpp"

// Standard library.
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {

    class AlignmentGraphVertex;
    class AlignmentGraphEdge;
    class AlignmentGraph;
    using AlignmentGraphBaseClass = CompactUndirectedGraph<
        AlignmentGraphVertex,
        AlignmentGraphEdge>;
    class Alignment;
    class AlignmentInfo;

    // Top level function to compute the marker alignment.
    void align(

        // Markers of the two oriented reads to be aligned, sorted by KmerId.
        const array<vector<MarkerWithOrdinal>, 2>& markers,

        // The maximum ordinal skip to be tolerated between successive markers
        // in the alignment.
        size_t maxSkip,

        // The maximum ordinal drift to be tolerated between successive markers
        // in the alignment. This is the drift of the two oriented reads
        // relative to each other.
        size_t maxDrift,

        // Marker frequency threshold.
        // When computing an alignment between two oriented reads,
        // marker kmers that appear more than this number of times
        // in either of the two oriented reads are discarded
        // (in both oriented reads).
        // Change to size_t when conversion completed.
        uint32_t maxMarkerFrequency,

        // Flag to control various types of debug output.
        bool debug,

        // The AlignmentGraph can be reused.
        // For performance, it should be reused when doing many alignments.
        AlignmentGraph&,

        // The computed alignment.
        // This should also be reused when performance is important.
        Alignment&,

        // Also create alignment summary information.
        AlignmentInfo&
        );
}



// Each vertex corresponds a pair of markers in the
// two oriented reads that have the same kmer.
class shasta::AlignmentGraphVertex {
public:

    // The KmerId of this marker.
    KmerId kmerId;

    // The index of this k-mer in each of the sorted markers vectors.
    array<size_t, 2> indexes;

    // The ordinals of this k-mer in each of the oriented reads.
    // This equals the index when the markers are sorted by position.
    array<size_t, 2> ordinals;

    // The position of this marker in each of the two sequences.
    array<int, 2> positions = array<int, 2>{
        std::numeric_limits<int>::min(), std::numeric_limits<int>::min()};

    // Data members used to find the shortest path.
    AlignmentGraphBaseClass::vertex_descriptor predecessor;
    uint64_t distance;
    uint8_t color;

    // Order by ordinal in the first sequence.
    bool operator<(const AlignmentGraphVertex& that) const
    {
        return ordinals[0] < that.ordinals[0];
    }
};



class shasta::AlignmentGraphEdge {
public:
    uint64_t weight;

    AlignmentGraphEdge(uint64_t weight) :
        weight(weight)
        {}
};



class shasta::AlignmentGraph : public AlignmentGraphBaseClass {
public:

    void create(
        const array<vector<MarkerWithOrdinal>, 2>&,
        uint32_t maxMarkerFrequency,
        size_t maxSkip,
        size_t maxDrift,
        bool debug,
        Alignment&,
        AlignmentInfo&);

private:

    // There is a vertex for each pair of markers with the same k-mer.
    // In addition, there is a start vertex and a finish vertex.
    vertex_descriptor vStart;
    vertex_descriptor vFinish;

    static void writeMarkers(
        const vector<MarkerWithOrdinal>&,
        const string& fileName
        );
    void createVertices(
        const array<vector<MarkerWithOrdinal>, 2>&,
        uint32_t maxMarkerFrequency);
    void writeVertices(const string& fileName) const;
    void createEdges(
        uint32_t markerCount0,
        uint32_t MarkerCount1,
        size_t maxSkip,
        size_t maxDrift);
    void writeEdges(const string& fileName) const;

    // Write in graphviz format, without the start and finish vertices.
    void writeGraphviz(const string& fileName) const;

    // Write an image representing the markers and the computed alignment
    // in 2-D ordinal space.
public:
    static void writeImage(
        const vector<MarkerWithOrdinal>&,
        const vector<MarkerWithOrdinal>&,
        const Alignment&,
        uint64_t markersPerPixel,
        const string& fileName);
private:

    // Data members used to find the shortest path.
    vector<vertex_descriptor> shortestPath;
    FindShortestPathQueue<AlignmentGraph> queue;
    void writeShortestPath(const string& fileName) const;

    // Flags that are set for markers whose k-mers
    // have frequency maxMarkerFrequency or less in
    // both oriented reads being aligned.
    // Indexed by [0 or 1][ordinal]
    // (this is the standard marker ordinal that counts all markers).
    array<vector<bool>, 2> isLowFrequencyMarker;

    // The corrected ordinals, keeping into account only low frequency markers.
    // Index by [01][ordinal].
    array<vector<uint32_t>, 2> correctedOrdinals;

};

#endif
