#ifndef SHASTA_BUBBLES_HPP
#define SHASTA_BUBBLES_HPP

/*******************************************************************************

Class to describe an analyze bubbles in the assembly graph.

*******************************************************************************/

#include "AssemblyGraph.hpp"

#include <boost/graph/adjacency_list.hpp>

#include <limits>
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

        // Concordant/discordant sums computed by flagBadBubbles.
        uint64_t concordantSum = 0;
        uint64_t discordantSum = 0;
        double discordantRatio() const
        {
            return double(discordantSum) / double(concordantSum + discordantSum);
        }

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
    void writeBubbles();


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



    // Bubble graph.
    // An undirected graph in which each vertex represents a Bubble.
    // Two vertices are joined by an undirected edge if the corresponding Bubbles
    // have one or more common OrientedReadIds.
    class BubbleGraphVertex {
    public:
        uint64_t bubbleId;
        BubbleGraphVertex(uint64_t bubbleId) :
            bubbleId(bubbleId) {}
        BubbleGraphVertex() : bubbleId(std::numeric_limits<uint64_t>::max()) {}
    };
    class BubbleGraphEdge {
    public:

        // Store the number of common oriented reads for each pair of
        // branches in the two bubbles.
        // matrix[sideA][sideB] stores the number of OrientedReadIds
        // that appear on sideA of the "first" bubble and on
        // sideB of the "second" bubble of this edge.
        // The "first" bubble of the edge is the lowered numbered.
        array<array<uint64_t, 2>, 2> matrix;

        BubbleGraphEdge()
        {
            for(uint64_t sideA=0; sideA<2; sideA++) {
                for(uint64_t sideB=0; sideB<2; sideB++) {
                    matrix[sideA][sideB] = 0;
                }
            }
        }
    };
    // We use boost::vecS for thew vertices, so vertex_descriptors are the same
    // as bubble ids (indices into the bubbles vector).
    using BubbleGraphBaseClass =
        boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS,
        BubbleGraphVertex, BubbleGraphEdge>;
    class BubbleGraph: public BubbleGraphBaseClass {
    public:
        // The vertex descriptor corresponding to each bubbleId.
        vector<vertex_descriptor> vertexTable;
        vertex_descriptor vertexDescriptor(uint64_t bubbleId) const
        {
            return vertexTable[bubbleId];
        }
    };
    BubbleGraph bubbleGraph;
    void createBubbleGraph();
    void writeBubbleGraphGraphviz() const;

    // Use the BubbleGraph to flag bad bubbles.
    void flagBadBubbles();
    void removeBadBubbles(double discordantRatioThreshold);



    const Assembler& assembler;

};

#endif

