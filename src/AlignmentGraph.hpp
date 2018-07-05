#ifndef CZI_SHASTA_ALIGNMENT_GRAPH_HPP
#define CZI_SHASTA_ALIGNMENT_GRAPH_HPP

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

// shasta
#include "CompactUndirectedGraph.hpp"
#include "Marker.hpp"
#include "shortestPath.hpp"

// Standard library.
#include "utility.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {

        class Alignment;
        class AlignmentInfo;

        class AlignmentGraphVertex;
        class AlignmentGraphEdge;
        class AlignmentGraph;
        using AlignmentGraphBaseClass = CompactUndirectedGraph<
            AlignmentGraphVertex,
            AlignmentGraphEdge>;


        // Top level function to compute the marker alignment.
        void align(

            // Markers of the two oriented reads to be aligned, sorted by KmerId.
            const vector<MarkerWithOrdinal>& markers0,
            const vector<MarkerWithOrdinal>& markers1,

            // The maximum ordinal skip to be tolerated between successive markers
            // in the alignment.
            size_t maxSkip,

            // The  maximum number of vertices in the alignment graph
            // that we allow a single k-mer to generate.
            size_t maxVertexCountPerKmer,

            // The AlignmentGraph can be reused.
            // For performance, it should be reused when doing many alignments.
            AlignmentGraph& graph,

            // Flag to control various types of debug output.
            bool debug,

            // The computed alignment.
            // This should also be reused when performance is important.
            Alignment&
            );


    }
}



// Each vertex corresponds a pair of markers in the
// two oriented reads that have the same kmer.
class ChanZuckerberg::shasta::AlignmentGraphVertex {
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



class ChanZuckerberg::shasta::AlignmentGraphEdge {
public:
    uint64_t weight;

    AlignmentGraphEdge(uint64_t weight) :
        weight(weight)
        {}
};



class ChanZuckerberg::shasta::AlignmentGraph : public AlignmentGraphBaseClass {
public:

    void create(
        const vector<MarkerWithOrdinal>& kmers0,
        const vector<MarkerWithOrdinal>& kmers1,
        size_t maxSkip,
        size_t maxVertexCountPerKmer,
        bool debug,
        Alignment&);

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
        const vector<MarkerWithOrdinal>&,
        const vector<MarkerWithOrdinal>&,
        size_t maxVertexCountPerKmer);
    void writeVertices(const string& fileName) const;
    void createEdges(
        uint32_t markerCount0,
        uint32_t MarkerCount1,
        size_t maxSkip);
    void writeEdges(const string& fileName) const;

    // Write in graphviz format, without the start and finish vertices.
    void writeGraphviz(const string& fileName) const;

    // Write an image representing the markers and the computed alignment
    // in 2-D ordinal space.
    void writeImage(
        const vector<MarkerWithOrdinal>&,
        const vector<MarkerWithOrdinal>&,
        const Alignment&,
        const string& fileName) const;

    // Data members used to find the shortest path.
    vector<vertex_descriptor> shortestPath;
    FindShortestPathQueue<AlignmentGraph> queue;
    void writeShortestPath(const string& fileName) const;

};

#endif
