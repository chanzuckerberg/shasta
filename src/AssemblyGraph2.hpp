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
        Branch(const MarkerGraphPath& path) : path(path) {}

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

    // Store read information on all branches.
    void storeReadInformation(const MarkerGraph&);

    // The reverse complement of this edge.
    // It contains the reverse complements of the bubbles of this edge,
    // in the same order.
    AssemblyGraph2BaseClass::edge_descriptor reverseComplement;

    // This constructor creates an edge without any paths.
    AssemblyGraph2Edge(uint64_t id) : id(id) {}

    // This constructor creates an edge with a single path.
    AssemblyGraph2Edge(uint64_t id, const MarkerGraphPath& path) :
        id(id), branches(1, Branch(path)) {}

    uint64_t ploidy() const {
        return branches.size();
    }

    bool isBubble() const
    {
        return ploidy() > 1;
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
    string colorByPeriod(uint64_t branchId) const;

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

    // This writes a gfa and a csv file with the given base name.
    void writeGfa(
        const string& baseName,
        bool writeSequence) const;
    void writeGfaBothStrands(
        const string& baseName,
        bool writeSequence) const;

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

    // Get the vertex descriptor for the vertex corresponding to
    // a given MarkerGraph::VertexId, creating the vertex if necessary.
    vertex_descriptor getVertex(MarkerGraph::VertexId);

    uint64_t nextEdgeId = 0;

    // Create a new edge corresponding to the given path.
    // Also create the vertices if necessary.
    edge_descriptor addEdge(const MarkerGraphPath&);

    // Assemble sequence for every marker graph path of every edge.
    void assemble();

    // Assemble sequence for every marker graph path of a given edge.
    void assemble(edge_descriptor);

    // Store GFA sequence in each edge.
    void storeGfaSequence();

    // Finds edges that form bubbles, then combine
    // each of them into a single edge with multiple paths.
    void gatherBubbles();
    edge_descriptor createBubble(
        vertex_descriptor v0,
        vertex_descriptor v1,
        const vector<edge_descriptor>&);

    // Find bubbles caused by copy number changes in repeats
    // with period up to maxPeriod.
    void findCopyNumberBubbles(uint64_t maxPeriod);

    // For each edge, compute the number of raw sequence bases
    // transfered in each direction for gfa output.
    void countTransferredBases();

    // Return true if an edge has id less than its reverse complement.
    bool idIsLessThanReverseComplement(edge_descriptor) const;

    // Return true if an edge has id greater than its reverse complement.
    bool idIsGreaterThanReverseComplement(edge_descriptor) const;

private:
    void checkReverseComplementEdges() const;

    // Store read information on all edges.
    void storeReadInformation();



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
        void removeWeakVertices(double discordantRatioThreshold);
        void writeHtml(const string& fileName);
        void computeConnectedComponents();
        vector< vector<BubbleGraph::vertex_descriptor> > connectedComponents;
    };
    BubbleGraph bubbleGraph;
    void createBubbleGraph(uint64_t readCount);



    // A predicate used to filter BubbleGraph edges for which
    // relativePhase() >= minRelativePhase.
    // This is only used for html output of the BubbleGraph.
    class BubbleGraphEdgePredicate {
    public:
        BubbleGraphEdgePredicate(
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

};



#endif

