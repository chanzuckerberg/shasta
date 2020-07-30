#ifndef SHASTA_MINI_ASSEMBLY_MARKER_GRAPH_HPP
#define SHASTA_MINI_ASSEMBLY_MARKER_GRAPH_HPP

#include "Kmer.hpp"
#include "ReadId.hpp"
#include "MarkerGraph2.hpp"

#include "utility.hpp"

namespace shasta {
    class MiniAssemblyMarkerGraph;
}



// A graph class used for mini-assemblies.
class shasta::MiniAssemblyMarkerGraph : public MarkerGraph2<KmerId, uint64_t> {
public:
    using Graph = MiniAssemblyMarkerGraph;
    MiniAssemblyMarkerGraph(const vector<OrientedReadId>& orientedReadIds) :
        orientedReadIds(orientedReadIds) {}
    void removeSelfEdges();
    void removeLowCoverageEdges(
        uint64_t minCoverage,
        uint64_t minPerStrandEdgeCoverage);
    void removeIsolatedVertices();

    vector<OrientedReadId> orientedReadIds;

    // A path is a sequence of consecutive edges.
    using Path = vector<edge_descriptor>;



    // A bubble is a set of linear paths that have the same
    // start and end vertex.
    class Bubble {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        vector<Path> branches;

        // Which oriented reads are represented in each branch:
        // Indexed by [sequenceId][branch], where seuquenceId
        // identifies an oriented read by its index into orientedReadIds.
        // contains[sequenceId][branch] is present internally to that branch.
        vector< vector<bool> > contains;

        // Vector that tells us whic branch each oriented read appears in, or:
        // -1 if the oriented read appears in no branch.
        // -2 if the oriented read appears in more than one branch.
        // Indexes by seuenceid (index into orientedReadIds vector).
        vector<int64_t> branchTable;
     };
    vector<Bubble> bubbles;
    void findBubbles();
    void fillBubbleContainment(Bubble&);
};

#endif
