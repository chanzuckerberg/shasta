#ifndef SHASTA_FORKS_HPP
#define SHASTA_FORKS_HPP

/*******************************************************************************

A forward Fork is a marker graph vertex with out-degree greater than 1.
Each of the outgoing edges is a Branch of the Fork.

A backward Fork is a marker graph vertex with in-degree greater than 1.
Each of the incoming edges is a Branch of the Fork.

Forks can occur because of mixing in the marker graph between haplotypes
and/or between copies of similar/repeated sequence.

An edge can belong to at most one forward Fork and one backward Fork.

We assume that the marker graph contains a minimal amount of
erroneous edges. This can be achieved by using sufficiently large values of
of minCoverage and minEdgeCoverage.

Forks are used to separate haplotypes (phasing) and to separate copies of
similar/repeated sequence.

*******************************************************************************/

#include "MarkerGraph.hpp"

namespace shasta {
    class Forks;
    class CompressedMarker;
}



class shasta::Forks {
public:
    using VertexId = MarkerGraph::VertexId;
    using EdgeId = MarkerGraph::EdgeId;

    enum class ForkDirection {
        Forward = 0,
        Backward = 1
    };
    static string directionString(ForkDirection direction)
    {
        switch(direction) {
        case ForkDirection::Forward: return "Forward";
        case ForkDirection::Backward: return "Backward";
        default: SHASTA_ASSERT(0);
        }
    }

    class Branch {
    public:

        // The edge that this branch consists of.
        EdgeId edgeId;

        // The MarkerIntervals of this edge.
        vector<MarkerInterval> markerIntervals;

        // Constructor.
        Branch(EdgeId, const MarkerGraph&);

        void write(ostream&) const;
    };



    class Fork {
    public:

        // The vertex that defines this fork.
        VertexId vertexId;

        // Flag that tells us whether this is a forward of backward Fork.
        ForkDirection direction;

        // The Branches of this Fork.
        // If this is a forward branch, all the branch edges have vertexId as their source.
        // If this is a backward branch, all the branch edges have vertexId as their target.
        vector<Branch> branches;

        // Constructor.
        Fork(VertexId, ForkDirection, const MarkerGraph&);

        void write(ostream&) const;
    };

    // All of the forks in the marker graph.
    vector<Fork> forks;



    // A data structure that tells us which branches each marker graph edge
    // belongs to. Sorted by EdgeId.
    class EdgeInfo {
    public:
        EdgeId edgeId;
        const Fork* fork;
        uint64_t branch;    // The index of the branch in the fork.
        bool operator<(const EdgeInfo& that) const
        {
            return edgeId < that.edgeId;
        }
        EdgeInfo(EdgeId edgeId, const Fork* fork, uint64_t branch) :
            edgeId(edgeId),
            fork(fork),
            branch(branch)
        {}
    };
    vector<EdgeInfo> edgeInfos;
    void constructEdgeInfos();



    // We also need a data structure that tells us the sequence of
    // branch edges encountered by each OrientedReadId.
    // Indexed by OrientedReadId::getValue().
    // For each OrientedReadId, sorted by ordinals.
    class OrientedReadEdgeInfo {
    public:
        const EdgeInfo* edgeInfo;
        array<uint32_t, 2> ordinals;
        bool operator<(const OrientedReadEdgeInfo& that) const
        {
            return ordinals[0]+ordinals[1] < that.ordinals[0]+that.ordinals[1];
        }
        OrientedReadEdgeInfo(
            const EdgeInfo* edgeInfo,
            array<uint32_t, 2> ordinals) :
            edgeInfo(edgeInfo), ordinals(ordinals)
        {}
    };
    vector< vector<OrientedReadEdgeInfo> > orientedReadEdgeInfos;
    void constructOrientedReadEdgeInfos();



    // Const references to objects we need.
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MarkerGraph& markerGraph;

    // Constructor.
    Forks(
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&);

    // Analyze a single Fork.
    void analyze(VertexId, ForkDirection, uint32_t maxDistance) const;

    // Find a Fork with a given VertexId and direction.
    const Fork* findFork(VertexId, ForkDirection) const;

    // Find Forks that are close to a given Fork,
    // and return them with their approximate marker offsets.
    void findNearbyForks(
        const Fork&,                    // Our starting Fork
        uint32_t maxDistance,           // The maximum distance in markers
        vector< pair<const Fork*, int32_t> >& // pairs(Fork, marker offset)
        ) const;
};



#endif

