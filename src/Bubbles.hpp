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
    class AlignmentData;
    class OrientedReadPair;
}

class shasta::Bubbles {
public:

    Bubbles(
        const Assembler&,
        bool debug = false
    );

private:



    // Information about an oriented read in a bubble.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // The side of the bubble that this oriented read appears in.
        // This is an index into aEdgeIds.
        // OrientedReadIds that appear on more than one side are not stored.
        uint32_t side;

        // The minimum and maximum ordinal at which this OrientedRead
        // appears on this side of the bubble.
        uint32_t minOrdinal;
        uint32_t maxOrdinal;

        OrientedReadInfo(
            OrientedReadId orientedReadId,
            uint32_t side,
            uint32_t minOrdinal,
            uint32_t maxOrdinal) :
            orientedReadId(orientedReadId),
            side(side),
            minOrdinal(minOrdinal),
            maxOrdinal(maxOrdinal)
            {}
    };



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
        // Values that are present on more than one side are removed.
        vector<OrientedReadInfo> orientedReadInfos;

        // Return the number of oriented reads on a given side.
        uint64_t countOrientedReadsOnSide(uint64_t side) const;

        // Concordant/discordant sums computed by flagBadBubbles.
        uint64_t concordantSum = 0;
        uint64_t discordantSum = 0;
        uint64_t sum() const
        {
            return concordantSum + discordantSum;
        }
        double discordantRatio() const
        {
            return double(discordantSum) / double(sum());
        }

        // This flag is set for bubbles flagged as bad and removed from the BubbleGraph.
        bool isBad = false;

        // Terminal bubbles flag.
        bool isTerminalBackward = false; // No bubble "precedes" this bubble.
        bool isTerminalForward  = false; // No bubble "follows"  this bubble.

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
    void writeBubblesDetails();

    // Given two bubbles, return the number of oriented reads
    // common between the two and that:
    // - Reach bubble0 before bubble1.
    // - Reach bubble1 before bubble0.
    void evaluateRelativePosition(
        const Bubble& bubble0,
        const Bubble& bubble1,
        uint64_t& order01Count, // Number of of oriented reads that encounter bubble0 before bubble1
        uint64_t& order10Count  // Number of of oriented reads that encounter bubble1 before bubble0
    ) const;


    // A data structure that tells us, for each OrientedReadId,
    // which sides of which bubbles that OrientedReadId appears in.
    // Indexed by OrientedReadId::getValue().
    // Contains pairs (bubbleId, side)
    // where bubbleId is an index into the bubbles vector
    // and side is 0 or 1.
    // For each OrientedReadId, sorted by bubbleId and side.
    // Note an oriented read cannot appear on both sides of a bubble,
    // by construction.
    vector< vector <pair <uint64_t, uint64_t> > > orientedReadsTable;
    void fillOrientedReadsTable();
    void writeOrientedReadsTable();

    // Given two OrientedReadIds, use the orientedReadsTable
    // to count the number of times they appear on the same
    // side or opposite sides of the same bubble.
    void findOrientedReadsRelativePhase(
        OrientedReadId,
        OrientedReadId,
        uint64_t& sameSideCount,
        uint64_t& oppositeSideCount
    ) const;



    // Find OrientedReadIds that appear in at least one bubble
    // together with a given OrientedReadId.
    // They are returned sorted.
    void findNeighborOrientedReadIds(
        OrientedReadId,
        vector<OrientedReadId>&
    ) const;


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

        uint64_t componentId = std::numeric_limits<uint64_t>::max();
        uint64_t color; // Used to compute connected components.
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
        uint64_t diagonalCount() const
        {
            return matrix[0][0] + matrix[1][1];
        }
        uint64_t offDiagonalCount() const
        {
            return matrix[0][1] + matrix[1][0];
        }
        uint64_t totalCount() const
        {
            return matrix[0][0] + matrix[1][1]+ matrix[0][1] + matrix[1][0];
        }
        uint64_t concordantCount() const
        {
            return max(diagonalCount(), offDiagonalCount());
        }
        uint64_t discordantCount() const
        {
            return min(diagonalCount(), offDiagonalCount());
        }

        // Ambiguity of the edge is 0 if discordantCount() is 0
        // and 1 if discordantCount() = totalCount()/2, in which case
        // discordantCount() = concordantCount();
        double ambiguity() const
        {
            return double(2 * discordantCount()) / double(totalCount());
        }

        // Return the relative phase implied by this edge, which is
        // +1 if offdiagonalCount() is 0 and
        // -1 if diagonalCount() is 0.
        double relativePhase() const
        {
            const double diagonalRatio = double(diagonalCount()) / double(totalCount());
            return 2. * diagonalRatio - 1.;
        }

        // Number of oriented reads that reach each of the two bubbles first.
        uint64_t orderABCount = std::numeric_limits<uint64_t>::max();
        uint64_t orderBACount = std::numeric_limits<uint64_t>::max();

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

        void computeConnectedComponents();
        vector< vector<vertex_descriptor> > connectedComponents;
    };

    BubbleGraph bubbleGraph;
    void createBubbleGraph();

    // Write the BubbleGraph in graphviz format, coloring
    // the bubbles by discordant ratio and the edges by ambiguity.
    void writeBubbleGraphGraphviz(const string& fileName) const;

    // Write a single component of the BubbleGraph in html/svg format.
    // To compute sfdp layout, only consider edges
    // for which relativePhase() >= minRelativePhase.
    void writeBubbleGraphComponentHtml(
        uint64_t componentId,
        const vector<BubbleGraph::vertex_descriptor>& component) const;

    // ComponentGraph is used by writeBubbleGraphComponentSvg.
    class ComponentGraphVertex {
    public:
        array<double, 2> position;
    };
    using ComponentGraph = boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, ComponentGraphVertex>;


    // Use the BubbleGraph to flag bad bubbles.
    void flagBadBubbles();
    void removeBadBubbles(double discordantRatioThreshold);

    void flagTerminalBubbles();



    // In the PhasingGraph, each vertex represents an OrientedReadId.
    // Two vertices are joined by an undirected edge
    // if the corresponding OrientedReadIds appear
    // together in at least one bubble.
    // There is no global PhasingGraph object.
    // We construct a separate PhasingGraph for each connected
    // component of the BubbleGraph.
    class PhasingGraphVertex {
    public:
        OrientedReadId orientedReadId;
        int64_t phase = 0;
        double eigenvectorComponent = 0.;
        array<double, 2> position;
    };
    class PhasingGraphEdge {
    public:
        uint64_t sameSideCount;
        uint64_t oppositeSideCount;
        double relativePhase() const
        {
            const double ratio = double(sameSideCount) / double(sameSideCount + oppositeSideCount); // In [0,1]
            return 2. * ratio - 1.; // In (-1, 1)
        }
    };
    using PhasingGraphBaseClass =
        boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS,
        PhasingGraphVertex, PhasingGraphEdge>;

    class PhasingGraph: public PhasingGraphBaseClass {
    public:
        std::map<OrientedReadId, vertex_descriptor> vertexMap;
        void createVertices(const vector<OrientedReadId>&);
        void phaseSpectral(bool debug);
        // Write in html/svg format.
        // To compute sfdp layout, only consider edges
        // for which relativePhase() >= minRelativePhase.
        void writeHtml(
            const string& fileName,
            double minRelativePhase) const;
        void writeEdges(
            const string& fileName) const;
    };
    void createPhasingGraph(
        const vector<OrientedReadId>&,
        PhasingGraph&) const;

    // A predicate used to filter PhasingGraph edges for which
    // relativePhase() >= minrelativePhase.
    class PhasingGraphEdgePredicate {
    public:
        PhasingGraphEdgePredicate(
            const PhasingGraph& phasingGraph,
            double minRelativePhase) :
            phasingGraph(&phasingGraph),
            minRelativePhase(minRelativePhase) {}

        const PhasingGraph* phasingGraph;
        double minRelativePhase;

        bool operator() (const PhasingGraph::edge_descriptor e) const
        {
            return (*phasingGraph)[e].relativePhase() >= minRelativePhase;
        }
    };


    // Top level function for phasing.
    void phase(double minRelativePhase);

    // The component and phase of each oriented read.
    // The phase is 0 or 1.
    // Indexed by OrientedRead::getValue().
    vector< pair<uint32_t, uint32_t> > orientedReadsPhase;

    // Given a connected component of the BubbleGraph,
    // find the OrientedReadIds that appear in it.
    // The OrientedReadIds are returned sorted.
    void findComponentOrientedReads(
        const vector<BubbleGraph::vertex_descriptor>&,
        vector<OrientedReadId>&
        ) const;



    // Backward terminal oriented reads are oriented reads that are near the beginning
    // of a phasing component. The first bubble they encounter is a backward terminal
    // bubble. We store the first ordinal where the oriented read encounters the bubble.
    // Similarly, forward terminal oriented reads are oriented reads that are near the end
    // of a phasing component. The last bubble they encounter is a forward terminal
    // bubble. We store the last ordinal where the oriented read encounters the bubble.
    // The vectors below are indexed by OrientedReadId::getValue().
    // The contain std::numeric_limits<uint32_t>::max() for oriented reads
    // that are not backward or forward terminal.
    vector<uint32_t> backwardTerminalOrdinal;
    vector<uint32_t> forwardTerminalOrdinal;
    void findTerminalOrdinals();


    // Functions used to decide if an alignment should be used.
public:
    bool allowAlignment(const AlignmentData&) const;
private:

    const Assembler& assembler;
    bool debug;

};

#endif

