#ifndef SHASTA_BUBBLES_HPP
#define SHASTA_BUBBLES_HPP

/*******************************************************************************

Class to describe an analyze bubbles in the assembly graph.

*******************************************************************************/

#include "AssemblyGraph.hpp"

#include "vector.hpp"

namespace shasta {
    class Bubbles;
    class Assembler;
}

class shasta::Bubbles {
public:

    Bubbles(
        const Assembler&
    );

private:



    // For now we only consider diploid bubbles, defined using the
    // following strict criteria:
    // - Source vertex v0 has out-degree 2.
    // - Target vertex v1 has in-degree 2.
    // - There are two parallel edges eA and eB, both v0->v1.
    class Bubble {
    public:

        // The bubble as seen in the AssemblyGraph.
        AssemblyGraph::VertexId av0;
        AssemblyGraph::VertexId av1;
        array<AssemblyGraph::EdgeId, 2> aEdgeIds;

        // For convenience also store the MarkerGraph vertex ids.
        MarkerGraph::VertexId mv0;
        MarkerGraph::VertexId mv1;

        // The OrientedReadIds on each of the two sides.
        // Sorted by OrientedReadId, without duplicates.
        // Values that are present on both sides are removed.
        array<vector<OrientedReadId>, 2> orientedReadIds;

        // Constructor.
        Bubble(
            AssemblyGraph::VertexId av0,
            AssemblyGraph::VertexId av1,
            AssemblyGraph::EdgeId eA,
            AssemblyGraph::EdgeId eB,
            const MarkerGraph&,
            const AssemblyGraph&);

    private:
        void fillInOrientedReadIds(
            const MarkerGraph&,
            const AssemblyGraph&);
    };
    vector<Bubble> bubbles;
    void findBubbles();


    // A data structure that tells us, for each OrientedReadId,
    // which sides of which bubbles that OrientedReadId appears in.
    // Indexed by OrientedReadId::getValue().
    // Contains pairs (bubbleId, side)
    // where bubbleId is an index into the bubbles vector
    // and side is 0 or 1.
    vector< vector <pair <uint64_t, uint64_t> > > orientedReadsTable;
    void fillOrientedReadsTable();
    void writeOrientedReadsTable();


    // Figure out if two sequences differ only by copy numbers in
    // a 2- or 3-base repeat.
    static bool isShortRepeatCopyNumberDifference(
        const vector<Base>&,
        const vector<Base>&);

    const Assembler& assembler;

};

#endif

