#ifndef SHASTA_ASSEMBLY_GRAPH2_HPP
#define SHASTA_ASSEMBLY_GRAPH2_HPP

// Assembly graph for assembly mode 2.

// Shasta.
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>
#include "string.hpp"
#include "vector.hpp"



namespace shasta {
    class AssemblyGraph2;
    class AssemblyGraph2Vertex;
    class AssemblyGraph2Edge;
    class AssemblyGraph2Statistics;
    class BubbleChain;
    class MarkerGraph;
    class Mode2AssemblyOptions;
    class ReadFlags;

    using AssemblyGraph2BaseClass =
        boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
        AssemblyGraph2Vertex, AssemblyGraph2Edge>;

    using MarkerGraphPath = vector<MarkerGraph::EdgeId>;

    class AssemblerInfo;

    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }
}



// Linear chains of bubbles in the AssemblyGraph2.
class shasta::BubbleChain {
public:
    vector<AssemblyGraph2BaseClass::edge_descriptor> edges;

    // The edges of a BubbleChain are partitioned in PhasingRegions.
    class PhasingRegion {
    public:
        uint64_t firstPosition;
        uint64_t lastPosition;
        bool isPhased;
        uint64_t componentId = std::numeric_limits<uint64_t>::max();   // Only valid if isPhased is true.
    };
    vector<PhasingRegion> phasingRegions;
};



class shasta::AssemblyGraph2Vertex {
public:
    MarkerGraph::VertexId markerGraphVertexId;

    AssemblyGraph2Vertex(MarkerGraph::VertexId markerGraphVertexId) :
        markerGraphVertexId(markerGraphVertexId) {}

    // The bubble chains that begin/end at this vertex.
    vector<BubbleChain const *> bubbleChainsBeginningHere;
    vector<BubbleChain const *> bubbleChainsEndingHere;
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
        Branch() {}

        // Assembled sequence.
        // This excludes the first and last k/2 RLE bases.
        vector<Base> rawSequence;

        // Sequence to be written to gfa.
        vector<Base> gfaSequence;

        // The distinct oriented reads present on edges of this branch.
        // Sorted.
        vector<OrientedReadId> orientedReadIds;

        // Minimum and sum of coverage on the marker graph graph
        // edges of this branch.
        uint64_t minimumCoverage = 0;
        uint64_t coverageSum = 0;
        uint64_t averageCoverage() const
        {
            return uint64_t(std::round(double(coverageSum) / double(path.size())));
        }

        // Fill in orientedReads and average/minimum coverage.
        void storeReadInformation(const MarkerGraph&);
    };
    vector<Branch> branches;

    // The strongest branch.
    uint64_t getStrongestBranchId() const;

    // Functions that remove some of the branches.
    void removeAllSecondaryBranches();
    void removeAllBranchesExceptStrongest();
    void forceMaximumPloidy(uint64_t);

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

    // This constructor creates an edge with two paths.
    AssemblyGraph2Edge(
        uint64_t id,
        const MarkerGraphPath& path0,
        bool containsSecondaryEdges0,
        const MarkerGraphPath& path1,
        bool containsSecondaryEdges1) :
        id(id), branches({Branch(path0, containsSecondaryEdges0), Branch(path1, containsSecondaryEdges1)}) {}

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

    uint64_t minimumPathLength() const
    {
        uint64_t length = std::numeric_limits<uint64_t>::max();
        for(const Branch& branch: branches) {
            length = min(length, uint64_t(branch.path.size()));
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
    static const uint64_t invalidComponentId = std::numeric_limits<uint64_t>::max();
    uint64_t componentId = invalidComponentId;
    static const uint64_t invalidPhase = std::numeric_limits<uint64_t>::max();
    uint64_t phase = invalidPhase;
    bool isPhased() const
    {
        return componentId != std::numeric_limits<uint64_t>::max();
    }

    // Bubble chain information.
    // If this edge is part of a bubble chain, this stores
    // pair(bubble chain pointer, position in bubble chain).
    // Otherwise, it stores pair(0, 0).
    pair<BubbleChain const *, uint64_t> bubbleChain = {0, 0};
};



class shasta::AssemblyGraph2 :
    public AssemblyGraph2BaseClass,
    public MultithreadedObject<AssemblyGraph2> {
public:

    // The constructor creates an edge for each linear path
    // in the marker graph. Therefore, immediately after construction,
    // each edge has a single MarkerGraphPath (no bubbles).
    AssemblyGraph2(
        uint64_t readRepresentation,
        uint64_t k, // Marker length
        const MemoryMapped::Vector<ReadFlags>& readFlags,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&,
        uint64_t pruneLength,
        const Mode2AssemblyOptions&,
        AssemblyGraph2Statistics&,
        size_t threadCount
        );

    void writeCsv(const string& baseName) const;
    void writeVerticesCsv(const string& baseName) const;
    void writeEdgesCsv(const string& baseName) const;
    void writeEdgeDetailsCsv(const string& baseName) const;

    // GFA and FASTA output.
    // These must be called after storeGfaSequence,
    // but writeDetailed can be caller earlier for some combinations of flags
    // (see writeDetailedEarly).
    void writeDetailed(
        const string& baseName,
        bool writeSequence,
        bool writeSequenceLengthInMarkers,
        bool writeCsv,
        bool writeGfa,
        bool writeFasta) const;
    void writeDetailedEarly(const string& baseName);
    void writeHaploid(
        const string& baseName,
        bool writeSequence,
        bool writeCsv,
        bool writeGfa,
        bool writeFasta,
        AssemblyGraph2Statistics* statistics = 0) const;
    void writePhased(
        const string& baseName,
        bool writeSequence,
        bool writeCsv,
        bool writeGfa,
        bool writeFasta,
        AssemblyGraph2Statistics* statistics = 0) const;

    // Hide a AssemblyGraph2BaseClass::Base.
    using Base = shasta::Base;

private:

    // Some shorthands for frequently used types.
    using G = AssemblyGraph2;
    using V = AssemblyGraph2Vertex;
    using E = AssemblyGraph2Edge;

    // Some Assembler data that we need.
    uint64_t readRepresentation;
    uint64_t k;
    const MemoryMapped::Vector<ReadFlags>& readFlags;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
public:
    uint64_t getReadCount() const
    {
        return markers.size() / 2;
    }
private:
    const MarkerGraph& markerGraph;

    // Map that gives us the vertex descriptor corresponding to
    // each marker graph vertex.
    std::map<MarkerGraph::VertexId, vertex_descriptor> vertexMap;

    // Initial creation of vertices and edges.
    void create();

    // Prune leafs in which all branches are shorter than the specified length.
    void prune(uint64_t pruneLength);



    // Superbubble removal.

    // This creates superbubbles using edges shorter than a length threshold (in markers).
    void handleSuperbubbles0(
        uint64_t edgeLengthThreshold,
        uint64_t maxSuperbubbleSize,
        uint64_t maxSuperbubbleChunkSize,
        uint64_t maxSuperbubbleChunkPathCount);

    // This creates superbubbles using all edges not in bubble chains.
    void handleSuperbubbles1(
        uint64_t maxSuperbubbleSize,
        uint64_t maxSuperbubbleChunkSize,
        uint64_t maxSuperbubbleChunkPathCount);

    class Superbubble;

    // This uses a dominator tree to find choking points
    // and partition the superbubble into chunks,
    // then does path enumeration on individual chunks.
    void handleSuperbubble1(
        Superbubble&,
        uint64_t maxSuperbubbleSize,
        uint64_t maxSuperbubbleChunkSize,
        uint64_t maxSuperbubbleChunkPathCount);



    // Remove short loop-back edges.
    void removeShortLoopbackEdges(uint64_t edgeLengthThreshold);

    // Remove degenerate branches
    void removeDegenerateBranches();

    // Get the vertex descriptor for the vertex corresponding to
    // a given MarkerGraph::VertexId, creating the vertex if necessary.
    vertex_descriptor getVertex(MarkerGraph::VertexId);

    // This is used to generate ids for edges (segments) and bubble chains.
    uint64_t nextId = 0;

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

    // Merge an edge with the previous or following edge, if possible.
    edge_descriptor mergeWithPreviousIfPossible(edge_descriptor);
    edge_descriptor mergeWithFollowingIfPossible(edge_descriptor);

    // Find linear chains of adjacent non-bubbles.
    // Used by merge.
    void findNonBubbleLinearChains(vector< vector<edge_descriptor> >&) const;

    // Find and create new bubbles.
    void createNewBubbles();

    void writePloidyHistogram(ostream&) const;



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
    void assembleParallel(uint64_t threadCount);
    void assembleThreadFunction(size_t threadId);
    class AssembleParallelData {
    public:
        vector<edge_descriptor> allEdges;
    };
    AssembleParallelData assembleParallelData;

    // Assemble sequence for every marker graph path of a given edge.
    void assemble(edge_descriptor);

    // Store GFA sequence in each edge.
    void storeGfaSequence();

    void hetSnpStatistics(
        uint64_t& transitionCount,
        uint64_t& transversionCount,
        uint64_t& nonSnpCount
    ) const;

    // Finds edges that form bubbles, then combine
    // each of them into a single edge with multiple paths.
    void gatherBubbles();
    edge_descriptor createBubble(
        vertex_descriptor v0,
        vertex_descriptor v1,
        const vector<edge_descriptor>&);

    void removeSecondaryBubbles(uint64_t secondaryEdgeCleanupThreshold);
    void removeWeakBranches(uint64_t strongBranchThreshold);

    void forceMaximumPloidy(uint64_t maxPloidy);

    // For each edge, compute the number of raw sequence bases
    // transfered in each direction for gfa output.
    void countTransferredBases();

private:

    // Store read information on all edges.
    void storeReadInformation();
    void storeReadInformationParallel(uint64_t threadCount);
    void storeReadInformationThreadFunction(size_t threadId);
    class StoreReadInformationParallelData {
    public:
        vector<edge_descriptor> allEdges;
    };
    StoreReadInformationParallelData storeReadInformationParallelData;


    // Linear chains of bubbles in the AssemblyGraph2.
    vector<BubbleChain> bubbleChains;
    void findBubbleChains();
    void clearBubbleChains();
    void writeBubbleChains();
    void findPhasingRegions();
    void findPhasingRegions(BubbleChain&);
    void writePhasingRegions();

    // Compute the gfa sequence of a bubble chain
    // by concatenating gfa sequence of the strongest branch of
    // each of this edges.
    void computeBubbleChainGfaSequence(
        const BubbleChain&,
        vector<Base>&
        ) const;

    // Compute the gfa sequence of an unphased region
    // by concatenating gfa sequence of the strongest branch of
    // each of this edges.
    void computeUnphasedRegionGfaSequence(
        const BubbleChain&,
        const BubbleChain::PhasingRegion&,
        vector<Base>&
        ) const;

    // Compute the gfa sequence of an haplotype of a phased region.
    void computePhasedRegionGfaSequence(
        const BubbleChain&,
        const BubbleChain::PhasingRegion&,
        uint64_t haplotype,
        vector<Base>&
        ) const;



    // A superbubble is a subset of the AssemblyGraph2.
    // Each vertex corresponds to a vertex of the AssemblyGraph2.
    // Each edge corresponds to a branch of the AssemblyGraph2.
    class SuperbubbleVertex;
    class SuperbubbleEdge;
    using SuperbubbleBaseClass =
        boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
        SuperbubbleVertex, SuperbubbleEdge>;

    class SuperbubbleVertex {
    public:
        // The AssemblyGraph2 vertex corresponding to this Superbubble vertex.
        AssemblyGraph2::vertex_descriptor av;

        SuperbubbleVertex(AssemblyGraph2::vertex_descriptor av = AssemblyGraph2::null_vertex()) : av(av) {}

        // The immediate dominator of this vertex (its parent
        // in the dominator tree).
        // This remains set to null_vertex for vertices unreachable from
        // the entrance.
        // immediateDominator0 is used for the forward dominator tree,
        // from the entrance to the exit.
        // immediateDominator1 is used for the backward dominator tree,
        // from the exit to the entrance.
        SuperbubbleBaseClass::vertex_descriptor immediateDominator0 =
            SuperbubbleBaseClass::null_vertex();
        SuperbubbleBaseClass::vertex_descriptor immediateDominator1 =
            SuperbubbleBaseClass::null_vertex();

        // And index that gives the position of this vertex on the critical path,
        // if this vertex is a on the critical path.
        // See AssemblyGraph2::handleSuperbubble1 for details.
        uint64_t positionInCriticalPath = std::numeric_limits<uint64_t>::max();
    };

    class SuperbubbleEdge {
    public:
        AssemblyGraph2::edge_descriptor ae;
        uint64_t branchId;
        SuperbubbleEdge(
            AssemblyGraph2::edge_descriptor ae,
            uint64_t branchId) :
            ae(ae), branchId(branchId) {}
        uint64_t chunk = std::numeric_limits<uint64_t>::max();
    };

    class Superbubble : public SuperbubbleBaseClass  {
    public:
        Superbubble(
            const AssemblyGraph2&,
            const vector<AssemblyGraph2::vertex_descriptor>&,
            uint64_t edgeLengthThreshold
        );
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
        void enumeratePaths(
            vertex_descriptor entrance,
            vertex_descriptor exit);

        // Return the number of distinct AssemblyGraph2 edges
        // that begin/end at a given vertex.
        uint64_t originalInDegree(vertex_descriptor) const;
        uint64_t originalOutDegree(vertex_descriptor) const;

        void writeGraphviz(ostream&, const AssemblyGraph2&) const;
        void writeGraphviz1(ostream&, const AssemblyGraph2&) const;

        // The critical path.
        // See AssemblyGraph2::handleSuperbubble1 for details.
        vector<vertex_descriptor> criticalPath;
        void computeCriticalPath();

        // Find the chunk that each edge belongs to.
        // This must be called after the dominator trees
        // and the critical path are computed.
        void findChunks();
        void findChunk(edge_descriptor);
        vector< vector<edge_descriptor> > chunkEdges;
    };



    // Iteratively remove bad bubbles using the PhasingGraph.
    void removeBadBubblesIterative(
        uint64_t minConcordantReadCount,
        uint64_t maxDiscordantReadCount,
        double minLogP,
        double epsilon,
        uint64_t superbubbleRemovalEdgeLengthThreshold,
        uint64_t maxSuperbubbleSize,
        uint64_t maxSuperbubbleChunkSize,
        uint64_t maxSuperbubbleChunkPathCount,
        uint64_t pruneLength,
        uint64_t componentSizeThreshold,
        size_t threadCount);

    // Hierarchical phasing using the PhasingGraph.
    void hierarchicalPhase(
        uint64_t minConcordantReadCount,
        uint64_t maxDiscordantReadCount,
        double minLogP,
        double epsilon,
        size_t threadCount);

    // Renumber component to make them contiguous starting at 0.
    void renumberComponents();

};



#endif

