#ifndef CZI_SHASTA_LOCAL_MARKER_GRAPH2_HPP
#define CZI_SHASTA_LOCAL_MARKER_GRAPH2_HPP


/*******************************************************************************

The local marker graph created by class LocalMarkerGraph2 is a subgraph
of the global marker graph, created by starting at a given vertex,
and extending out to a specified distance in both directions.
Distance is number of edges on the global marker graph.

Like in the global marker graph, each vertex corresponds to
a group of aligned markers.



CLASS LocalMarkerGraph VERSUS CLASS LocalMarkerGraph2.

Assembler::createLocalMarkerGraph uses class LocalMarkerGraph
to create a local marker graph from a set of oriented reads,
without using the global marker graph.
This is computationally expensive as it requires
computations of alignments and merging of vertices.

Assembler::extractLocalMarkerGraph uses class LocalMarkerGraph2
to create a subgraph of the global marker graph.
This is computationally inexpensive as most of the computing
was done when the global marker graph was created.

*******************************************************************************/

// shasta.
#include "Kmer.hpp"
#include "MarkerId.hpp"
#include "MemoryAsContainer.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "iostream.hpp"
#include <map>
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {

        // Forward declaration of types declared in this file.
        class LocalMarkerGraph2Vertex;
        class LocalMarkerGraph2Edge;
        class LocalMarkerGraph2;
        using LocalMarkerGraph2BaseClass = boost::adjacency_list<
            boost::setS,
            boost::listS,
            boost::bidirectionalS,
            LocalMarkerGraph2Vertex,
            LocalMarkerGraph2Edge
            >;

        // Forward declarations of classes defined elsewhere.
        class CompressedMarker;
        class LongBaseSequences;
        namespace MemoryMapped {
            template<class T> class Vector;
            template<class Int, class T> class VectorOfVectors;
        }
    }
}



class ChanZuckerberg::shasta::LocalMarkerGraph2Vertex {
public:

    // The global vertex id of the vertex of the global marker
    // graph that corresponds to this vertex.
    GlobalMarkerGraphVertexId vertexId;

    // The distance from the start vertex.
    int distance;

    // The markers of this vertex.
    class MarkerInfo {
    public:
        MarkerId markerId;
        OrientedReadId orientedReadId;
        uint32_t ordinal;
    };
    vector<MarkerInfo> markerInfos;

    LocalMarkerGraph2Vertex(
        GlobalMarkerGraphVertexId vertexId,
        int distance) :
        vertexId(vertexId),
        distance(distance)
        {}
};



class ChanZuckerberg::shasta::LocalMarkerGraph2Edge {
public:

    // Class to describe the intervening sequence between
    // the two markers that define the edge.
    class Sequence {
    public:
        // The number of overlapping bases between the
        // marker of the source vertex and the
        // marker of the target vertex.
        uint8_t overlappingBaseCount;

        // The intervening sequence between the two markers.
        vector<Base> sequence;

        // There are three possible cases:
        // overlappingBaseCount>0 && sequence.empty():
        //     The sequence of the two markers overlap by
        //     overlappingBaseCount bases.
        // overlappingBaseCount==0 && !sequence.empty():
        //     Sequence stores the base sequence between the two markers.
        // overlappingBaseCount==0 && sequence.empty():
        //     The two markers immediately follow each other,
        //     without any intervening sequence.

        bool operator<(const Sequence& that) const
        {
            return tie(overlappingBaseCount, sequence) <
                tie(that.overlappingBaseCount, that.sequence);
        }
    };


    // Each Info object corresponds to an oriented read
    // that appears at startOrdinal in the source vertex
    // of this edge and at startOrdinal+1 in the target vertex.
    class Info {
    public:
        OrientedReadId orientedReadId;

        // The ordinals in the source and target vertex.
        array<uint32_t, 2> ordinals;

        Info(OrientedReadId orientedReadId, uint32_t ordinal0, uint32_t ordinal1) :
            orientedReadId(orientedReadId)
        {
            ordinals[0] = ordinal0;
            ordinals[1] = ordinal1;
        }
    };

    // The oriented vertices of this edge, grouped by sequence.
    // Sorted by decreasing number of supporting reads.
    vector< pair<Sequence, vector<Info> > > infos;

    // Consensus is the number of reads supporting the
    // strongest sequence.
    size_t consensus() const
    {
        if(infos.empty()) {
            return 0;
        } else {
            return infos.front().second.size();
        }
    }

    // Coverage is the total number of reads supporting this edge,
    // with any sequence.
    size_t coverage() const
    {
        size_t c = 0;
        for(const auto& p: infos) {
            c += p.second.size();
        }
        return c;
    }

};



class ChanZuckerberg::shasta::LocalMarkerGraph2 :
    public LocalMarkerGraph2BaseClass {
public:

    LocalMarkerGraph2(
        uint32_t k,
        LongBaseSequences& reads,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MemoryMapped::Vector<CompressedGlobalMarkerGraphVertexId>& globalMarkerGraphVertex
        );

    // Find out if a vertex with the given GlobalMarkerGraphVertexId exists.
    // If it exists, return make_pair(true, v).
    // Otherwise, return make_pair(false, null_vertex());
    pair<bool, vertex_descriptor> findVertex(GlobalMarkerGraphVertexId) const;

    // Add a vertex with the given GlobalMarkerGraphVertexId
    // and return its vertex descriptor.
    // A vertex with this GlobalMarkerGraphVertexId must not exist.
    vertex_descriptor addVertex(
        GlobalMarkerGraphVertexId,
        int distance,
        MemoryAsContainer<MarkerId> markers);

    // Get the KmerId for a vertex.
    KmerId getKmerId(vertex_descriptor) const;

    // Store sequence information in the edge.
    void storeEdgeInfo(edge_descriptor);

    // Write in Graphviz format.
    // There are two types of Graphviz output:
    // - Detailed: includes vertex and edge labels, uses dot layout.
    //   Dot layout becomes too slow when the graph has more than
    //   a few thousand vertices.
    // - Compact: no labels, vertices rendered as points, uses sfdp layout.
    //   Sfdp layout becomes slow when the graph has more than
    //   several tens of thousand vertices.
    // Typical graphviz command for rendering for both types:
    // dot -O -T svg fileName -Gsize=200
    // The graph contains a layout attribute, so this will invoke
    // the correct layout algorithm automatically.
    void write(
        ostream&,
        size_t minCoverage,
        size_t minConsensus,
        int maxDistance,
        bool detailed) const;
    void write(
        const string& fileName,
        size_t minCoverage,
        size_t minConsensus,
        int maxDistance,
        bool detailed) const;

private:

    // Map a global vertex id to a vertex descriptor for the local graph.
    std::map<GlobalMarkerGraphVertexId, vertex_descriptor> vertexMap;

    // The length of k-mers used as markers.
    uint32_t k;

    // Reference to the global data structure containing all reads and markers
    // (not just those in this local marker graph).
    LongBaseSequences& reads;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;

    // A reference to the vector containing the global marker graph vertex id
    // corresponding to each marker.
    // Indexed by MarkerId.
    const MemoryMapped::Vector<CompressedGlobalMarkerGraphVertexId>& globalMarkerGraphVertex;

    class Writer {
    public:
        Writer(
            const LocalMarkerGraph2&,
            size_t minCoverage,
            size_t minConsensus,
            int maxDistance,
            bool detailed);
        void operator()(ostream&) const;
        void operator()(ostream&, vertex_descriptor) const;
        void operator()(ostream&, edge_descriptor) const;
        const LocalMarkerGraph2& graph;
        size_t minCoverage;
        size_t minConsensus;
        int maxDistance;
        bool detailed;
    };
    friend class Writer;
};

#endif
