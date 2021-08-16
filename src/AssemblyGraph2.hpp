#ifndef SHASTA_ASSEMBLY_GRAPH2_HPP
#define SHASTA_ASSEMBLY_GRAPH2_HPP

// Assembly graph for assembly mode2.



// Shasta.
#include "Marker.hpp"
#include "MarkerGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "array.hpp"
#include <map>
#include "string.hpp"
#include "vector.hpp"



namespace shasta {
    class AssemblyGraph2;
    class AssemblyGraph2Vertex;
    class AssemblyGraph2Edge;
    class MarkerGraph;

    using AssemblyGraph2BaseClass =
        boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
        AssemblyGraph2Vertex, AssemblyGraph2Edge>;

    using MarkerGraphPath = vector<MarkerGraph::EdgeId>;

    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }
}



class shasta::AssemblyGraph2Vertex {
public:
    MarkerGraph::VertexId markerGraphVertexId;

    AssemblyGraph2Vertex(MarkerGraph::VertexId markerGraphVertexId) :
        markerGraphVertexId(markerGraphVertexId) {}

    // The bubble chains that begin/end at this vertex.
    vector<vector<AssemblyGraph2BaseClass::edge_descriptor> const *> bubbleChainsBeginningHere;
    vector<vector<AssemblyGraph2BaseClass::edge_descriptor> const *> bubbleChainsEndingHere;
};



class shasta::AssemblyGraph2Edge {
public:

    // Id used for gfa output.
    uint64_t id;



    // Each assembly graph edge corresponds to
    // a set of paths in the marker graph.
    // This way it can describe a bubble in the marker graph.

    // Class to describe a single branch.
    class Branch {
    public:
        MarkerGraphPath path;
        bool containsSecondaryEdges;
        Branch(const MarkerGraphPath& path, bool containsSecondaryEdges) :
            path(path),
            containsSecondaryEdges(containsSecondaryEdges) {}

        // Assembled sequence.
        // This excludes the first and last k/2 RLE bases.
        vector<Base> rawSequence;

        // Sequence to be written to gfa.
        vector<Base> gfaSequence;

        // The distinct oriented reads present on edges of this branch.
        // Sorted.
        vector<OrientedReadId> orientedReadIds;

        // Minimum and average coverage on the marker graph graph
        // edges of this branch.
        uint64_t minimumCoverage;
        uint64_t averageCoverage;

        // Fill in orientedReads and average/minimum coverage.
        void storeReadInformation(const MarkerGraph&);
    };
    vector<Branch> branches;

    // The strongest branch.
    uint64_t strongestBranchId;
    void findStrongestBranch();

    // Functions that remove some of the branches.
    void removeAllSecondaryBranches();
    void removeAllBranchesExceptStrongest();

    // Store read information on all branches.
    void storeReadInformation(const MarkerGraph&);

    // This constructor creates an edge without any paths.
    AssemblyGraph2Edge(uint64_t id) : id(id) {}

    // This constructor creates an edge with a single path.
    AssemblyGraph2Edge(
        uint64_t id,
        const MarkerGraphPath& path,
        bool containsSecondaryEdges) :
        id(id), branches(1, Branch(path, containsSecondaryEdges)) {}

    uint64_t ploidy() const {
        return branches.size();
    }

    bool isBubble() const
    {
        return ploidy() > 1;
    }

    bool isDegenerateBubble() const
    {
        return
            (ploidy() == 2)
            and
            (branches[0].rawSequence == branches[1].rawSequence);
    }

    // Construct a string to id each of the markerGraphPaths.
    string pathId(uint64_t branchId) const
    {
        string s = to_string(id);
        if(isBubble()) {
            s.append("." + to_string(branchId));
        }
        return s;
    }

    uint64_t maximumPathLength() const
    {
        uint64_t length = 0;
        for(const Branch& branch: branches) {
            length = max(length, uint64_t(branch.path.size()));
        }
        return length;
    }

    // Return the number of raw bases of sequence identical between
    // all branches at the beginning/end.
    uint64_t countCommonPrefixBases() const;
    uint64_t countCommonSuffixBases() const;

    // The number of raw sequence bases transfered
    // in each direction for gfa output.
    uint64_t backwardTransferCount = 0;
    uint64_t forwardTransferCount = 0;

    // Figure out if this is a bubble is caused by copy number
    // differences in repeats of period up to maxPeriod.
    // If this is the case, stores the shortest period for which this is true.
    // Otherwise, stores 0 as the period.
    void computeCopyNumberDifferencePeriod(uint64_t maxPeriod);
    uint64_t period = 0;
    string color(uint64_t branchId) const;

    // This flag is set for bubbles classified as bad and
    // removed from the bubble graph.
    bool isBad = false;

    // Phasing information - only set for bubbles not marker
    // as bad by the above flag.
    uint64_t componentId = std::numeric_limits<uint64_t>::max();
    uint64_t phase = std::numeric_limits<uint64_t>::max();

    // Bubble chain information.
    // If this edge is part of a bubble chain, this stores
    // pair(bubble chain pointer, position in bubble chain).
    // Otherwise, it stores pair(0, 0).
    pair<const vector<AssemblyGraph2BaseClass::edge_descriptor>*, uint64_t> bubbleChain = {0, 0};
};



class shasta::AssemblyGraph2 : public AssemblyGraph2BaseClass {
public:

    // The constructor creates an edge for each linear path
    // in the marker graph. Therefore, immediately after construction,
    // each edge has a single MarkerGraphPath (no bubbles).
    AssemblyGraph2(
        uint64_t k, // Marker length
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&);

    void writeCsv(const string& baseName) const;
    void writeVerticesCsv(const string& baseName) const;
    void writeEdgesCsv(const string& baseName) const;
    void writeEdgeDetailsCsv(const string& baseName) const;

    // Gfa output includes a csv file with the given base name.
    void writeGfa(
        const string& baseName,
        bool writeSequence = true);
    // Haploid gfa and no csv file.
    void writeHaploidGfa(
        const string& baseName,
        bool writeSequence = true);

    // Hide a AssemblyGraph2BaseClass::Base.
    using Base = shasta::Base;

private:

    // Some shorthands for frequently used types.
    using G = AssemblyGraph2;
    using V = AssemblyGraph2Vertex;
    using E = AssemblyGraph2Edge;

    // Some Assembler data that we need.
    uint64_t k;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MarkerGraph& markerGraph;

    // Map that gives us the vertex descriptor corresponding to
    // each marker graph vertex.
    std::map<MarkerGraph::VertexId, vertex_descriptor> vertexMap;

    // Initial creation of vertices and edges.
    void create();

    // Remove secondary edges making sure to not introduce any dead ends.
    void cleanupSecondaryEdges();

    // Handle superbubbles.
    void handleSuperbubbles(uint64_t edgeLengthThreshold);

    // Remove degenerate branches
    void removeDegenerateBranches();

    // Get the vertex descriptor for the vertex corresponding to
    // a given MarkerGraph::VertexId, creating the vertex if necessary.
    vertex_descriptor getVertex(MarkerGraph::VertexId);

    uint64_t nextEdgeId = 0;

    // Create a new edge corresponding to the given path.
    // Also create the vertices if necessary.
    edge_descriptor addEdge(
        const MarkerGraphPath&,
        bool containsSecondaryEdges);

    // Merge consecutive non-bubbles, when possible.
     void merge(
        bool storeReadInformation,  // If true, store read information for merged edges.
        bool assemble               // If true, assemble merged edges.
        );
    edge_descriptor merge(
        const vector<edge_descriptor>&,
        bool storeReadInformation,  // If true, store read information for merged edges.
        bool assemble               // If true, assemble merged edges.
        );

    // Find linear chains of adjacent non-bubbles.
    // Used by merge.
    void findNonBubbleLinearChains(vector< vector<edge_descriptor> >&) const;



    // Predicate used to select non-bubble edges.
    class IsNonBubbleEdge {
    public:
        IsNonBubbleEdge(const G& g) : g(&g) {}
        IsNonBubbleEdge() : g(0) {}

        const G* g;

        bool operator() (const edge_descriptor e) const
        {
            return not (*g)[e].isBubble();
        }
    };



    // Assemble sequence for every marker graph path of every edge.
    void assemble();

    // Assemble sequence for every marker graph path of a given edge.
    void assemble(edge_descriptor);

    // Store GFA sequence in each edge.
    void storeGfaSequence();

    void hetSnpStatistics(
        uint64_t& transitionCount,
        uint64_t& transversionCount
    ) const;

    // Finds edges that form bubbles, then combine
    // each of them into a single edge with multiple paths.
    void gatherBubbles();
    edge_descriptor createBubble(
        vertex_descriptor v0,
        vertex_descriptor v1,
        const vector<edge_descriptor>&);

    void removeSecondaryBubbles();

    // Find/remove bubbles caused by copy number changes in repeats
    // with period up to maxPeriod.
    void findCopyNumberBubbles(uint64_t maxPeriod);
    void removeCopyNumberBubbles();

    // Remove bubbles marked isBad during phasing.
    // Only keep the strongest branch for each.
    void removeBadBubbles();

    // For each edge, compute the number of raw sequence bases
    // transfered in each direction for gfa output.
    void countTransferredBases();

private:

    // Store read information on all edges.
    void storeReadInformation();


    // Linear chains of bubbles in the AssemblyGraph2.
    vector< vector<edge_descriptor> > bubbleChains;
    void findBubbleChains();

    // Return the gfa id of a bubble chain.
    string bubbleChainId(const vector<edge_descriptor>&) const;

    // Compute the gfa sequence of a bubble chain
    // by concatenating gfa sequence of the strongest branch of
    // each of tis edges.
    void computeBubbleChainGfaSequence(
        const vector<edge_descriptor>&,
        vector<Base>&
        ) const;



    // A superbubble is a subset of the AssemblyGraph2.
    // Each vertex corresponds to a vertex of the AssemblyGraph2.
    // Each edge corresponds to a branch of the AssemblyGraph2.
    class SuperbubbleEdge {
    public:
        AssemblyGraph2::edge_descriptor ae;
        uint64_t branchId;
        SuperbubbleEdge(
            AssemblyGraph2::edge_descriptor ae,
            uint64_t branchId) :
            ae(ae), branchId(branchId) {}
    };
    using SuperbubbleBaseClass =
        boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
        AssemblyGraph2::vertex_descriptor, SuperbubbleEdge>;
    class Superbubble : public SuperbubbleBaseClass  {
    public:
        Superbubble(
            const AssemblyGraph2&,
            const vector<AssemblyGraph2::vertex_descriptor>&
        );
        vector<Superbubble::vertex_descriptor> entrances;
        vector<Superbubble::vertex_descriptor> exits;
        bool isSimpleLinearChain() const;

        // Enumerate paths from the entrance to the exit.
        // There must be exactly one entrance and one exit.
        vector< vector<edge_descriptor> > paths;
        void enumeratePaths();

        // Return the number of distinct AssemblyGraph2 edges
        // that begin/end at a given vertex.
        uint64_t originalInDegree(vertex_descriptor) const;
        uint64_t originalOutDegree(vertex_descriptor) const;

        void write(ostream&, const AssemblyGraph2&) const;
    };



    // Graph used for diploid phasing of the bubbles.
    // Each vertex corresponds to a diploid bubble in the Assembly2Graph.

    class BubbleGraphVertex {
    public:
        AssemblyGraph2::edge_descriptor e;
        uint64_t id;    // The same as the id of the AssemblyGraph2Edge.
        BubbleGraphVertex(
            AssemblyGraph2::edge_descriptor,
            const AssemblyGraph2Edge&);

        // The vertex stores OrientedReadIds that appear on one side but not
        // on the opposite side of the bubble.
        // Stored as pairs (OrientedReadId, side).
        vector< pair<OrientedReadId, uint64_t> > orientedReadIds;

        // The connected component this vertex belongs to.
        uint64_t componentId;

        // The phase assigned to this vertex (bubble)
        // It is only meaningful within each connected component.
        static const uint64_t invalidPhase = std::numeric_limits<uint64_t>::max();
        uint64_t phase = invalidPhase;

        // Layout position for display.
        array<double, 2> position;



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
            for(uint64_t i=0; i<2; i++) {
                for(uint64_t j=0; j<2; j++) {
                    matrix[i][j] = 0;
                }
            }
        }

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
            return diagonalCount() + offDiagonalCount();
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
    };



    using BubbleGraphBaseClass =
        boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS,
        BubbleGraphVertex, BubbleGraphEdge>;

    class BubbleGraph: public BubbleGraphBaseClass {
    public:
        // A table that, for each OrientedReadId, contains a list of
        // pairs(vertex, side) that the OrientedReadId appears on.
        // Indexed by OrientedReadId::getValue().
        vector< vector< pair<BubbleGraph::vertex_descriptor, uint64_t> > > orientedReadsTable;
        void createOrientedReadsTable(uint64_t readCount);
        void createEdges();
        void writeEdgesCsv(const string& fileName) const;
        void removeWeakEdges(uint64_t minReadCount);
        double discordantRatio(vertex_descriptor) const;
        void removeWeakVertices(
            double discordantRatioThreshold,
            vector<AssemblyGraph2::edge_descriptor>& badBubbles);

        // Pre-phasing.
        void writeHtml(const string& fileName);
        void writeGraphviz(const string& fileName) const;

        // Post-phasing.
        void writeHtml(
            const string& fileName,
            const std::set<edge_descriptor>& treeEdges,
            bool onlyShowUnhappyEdges);

        void computeConnectedComponents();
        vector< vector<BubbleGraph::vertex_descriptor> > connectedComponents;

        // Create a new BubbleGraph from a given connected component.
        void extractComponent(uint64_t componentId, BubbleGraph&) const;

        // Return true if the give edge has relative phase consistent
        // with the phases assigned to its two vertices.
        bool edgeIsHappy(edge_descriptor e) const;
    };
    BubbleGraph bubbleGraph;
    void createBubbleGraph(uint64_t readCount);
    void cleanupBubbleGraph(
        uint64_t minReadCount,
        double discordantRatioThreshold,
        double ambiguityThreshold);



    // A predicate used to filter BubbleGraph edges for which
    // relativePhase() >= minRelativePhase.
    // This is only used for html output of the BubbleGraph.
    class BubbleGraphEdgePredicate1 {
    public:
        BubbleGraphEdgePredicate1(
            const BubbleGraph& bubbleGraph,
            double minRelativePhase) :
            bubbleGraph(&bubbleGraph),
            minRelativePhase(minRelativePhase) {}

        const BubbleGraph* bubbleGraph;
        double minRelativePhase;

        bool operator() (const BubbleGraph::edge_descriptor e) const
        {
            return (*bubbleGraph)[e].relativePhase() >= minRelativePhase;
        }
    };



    // A predicate used to filter BubbleGraph spanning tree edges
    // and edges for which relativePhase() >= minRelativePhase.
    // This is only used for graphviz output of each connectec component, post-phasing.
    class BubbleGraphEdgePredicate2 {
    public:
        BubbleGraphEdgePredicate2(
            const BubbleGraph& bubbleGraph,
            double minRelativePhase,
            const std::set<BubbleGraph::edge_descriptor>& treeEdges) :
            bubbleGraph(&bubbleGraph),
            minRelativePhase(minRelativePhase),
            treeEdges(&treeEdges) {}

        const BubbleGraph* bubbleGraph;
        double minRelativePhase;
        const std::set<BubbleGraph::edge_descriptor>* treeEdges;

        bool operator() (const BubbleGraph::edge_descriptor e) const
        {
            return
                (*bubbleGraph)[e].relativePhase() >= minRelativePhase
                or
                treeEdges->find(e) != treeEdges->end();
        }
    };



    // Use each connected component of the bubble graph to phase the bubbles.
    void phase();

};



#endif

