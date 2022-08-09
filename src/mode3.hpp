#ifndef SHASTA_MODE3_HPP
#define SHASTA_MODE3_HPP

/*******************************************************************************

Class mode3::AssemblyGraph is the class used for Mode 3 assembly.
Using GFA terminology, the graph consists of Segments and Links.

A Segment corresponds to a linear sequence of edges, without branches,
in the marker graph.

If an oriented read enters segment 1 immediately after exiting segment 0,
we say that there is a transition 0->1. If there is a sufficient
number of transitions 0->1, we create a link 0->1.

*******************************************************************************/

// Shasta.
#include "hashArray.hpp"
#include "invalid.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "array.hpp"
#include "tuple.hpp"
#include "unordered_map"
#include "vector.hpp"



namespace shasta {
    namespace mode3 {
        class AssemblyGraph;
        class AssemblyGraphJourneyEntry;
        class MarkerGraphJourneyEntry;
        class AssemblyGraphJourneyInterval;
        class AssemblyPath;

    }

    // Some forward declarations of classes in the shasta namespace.
    class CompressedMarker;
    class MarkerGraph;

    extern template class MultithreadedObject<mode3::AssemblyGraph>;
}



// The marker graph journey of an oriented read is the sequence
// of marker graph edges it encounters.
// (An oriented read encounters a marker graph edge
// if the oriented read appears in the marker intervals for the edge).
// The marker graph journey of an oriented read is not necessarily
// a path in the marker graph because the oriented read
// can "skip" marker graph edges due to errors.
// In other places in Shasta, journeys are called "pseudopaths".
// We describe the marker graph journey of each oriented read as a sequence
// of MarkerGraphJourneyEntry objects.
// The MarkerGraphJourneyEntry identifies a marker graph edge
// by the segmentId  in the assembly graph and the position in that segment
// (that is, the first marker graph in the segment is at
// position 0, and so on).
class shasta::mode3::MarkerGraphJourneyEntry {
public:
    uint64_t segmentId;
    uint32_t position;
    array<uint32_t, 2> ordinals;

    bool operator<(const MarkerGraphJourneyEntry& that) const
    {
        return ordinals[0] < that.ordinals[0];
    }
    bool operator==(const MarkerGraphJourneyEntry& that) const
    {
        return
            tie(segmentId, position, ordinals) ==
            tie(that.segmentId, that.position, that.ordinals);
    }
};



// The assembly graph journey of an oriented read is the sequence
// of assembly graph segments (vertices) it encounters.
// The journey on an oriented read is not necessarily
// a path in the assembly graph because the oriented read
// can "skip" segments due to errors.
// We store the assembly graph journey of each oriented read as a sequence
// of AssemblyGraphJourneyEntry objects.
// The AssemblyGraphJourneyEntry stores the segment id and
// the first and last MarkerGraphJourneyEntry objects
// on the segment for the given oriented read.
// Indexed by OrientedReadId::getValue().
// Note a segmentId can appear more than once in the assembly
// graph journey of an oriented read. This can happen
// if the oriented read "goes around" in a tangle caused by repeats.
class shasta::mode3::AssemblyGraphJourneyEntry {
public:
    uint64_t segmentId;

    // The first and last MarkerGraphJourneyEntry that contributed to this
    // AssemblyGraphJourneyEntry.
    array<MarkerGraphJourneyEntry, 2> markerGraphJourneyEntries;
};



// A portion of the assembly graph journey of an oriented read.
class shasta::mode3::AssemblyGraphJourneyInterval {
public:

    OrientedReadId orientedReadId;

    // The first and last position in the assembly graph journey
    // for this oriented read.
    uint64_t first;
    uint64_t last;


    bool operator<(const AssemblyGraphJourneyInterval& that) const
    {
        return tie(orientedReadId, first) < tie(that.orientedReadId, that.first);
    }

};



// An assembly path in the AssemblyGraph
class shasta::mode3::AssemblyPath {
public:

    // The segments on the path.
    // The bool is true for reference segments.
    // The first and last segment are always reference segments.
    // A reference segment is one that is believed to be exclusive
    // to the sequence copy described by this path (that is,
    // it does not appear in other copies or haplotypes).
    vector< pair<uint64_t, bool> > segments;

    void clear()
    {
        segments.clear();
    }
};



// The AssemblyGraph is used to store the Mode 3 assembly graph,
// when it no longer needs to be changed,
// in memory mapped data structures.
class shasta::mode3::AssemblyGraph :
    public MultithreadedObject<AssemblyGraph> {
public:

    // Initial construction.
    AssemblyGraph(
        const string& largeDataFileNamePrefix,
        size_t largeDataPageSize,
        size_t threadCount,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&);

    // Constructor from binary data.
    AssemblyGraph(
        const string& largeDataFileNamePrefix,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&);

    // Data and functions to handle memory mapped data.
    const string& largeDataFileNamePrefix;
    size_t largeDataPageSize;
    string largeDataName(const string&) const;
    template<class T> void createNew(T& t, const string& name)
    {
        t.createNew(largeDataName(name), largeDataPageSize);
    }
    template<class T> void accessExistingReadOnly(T& t, const string& name)
    {
        t.accessExistingReadOnly(largeDataName(name));
    }

    // References to Assembler objects.
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MarkerGraph& markerGraph;

    uint64_t readCount() const
    {
        return markers.size() / 2;
    }

    // Each  linear chain of marker graph edges generates a segment.
    // The marker graph path corresponding to each segment is stored
    // indexed by segment id.
    MemoryMapped::VectorOfVectors<MarkerGraphEdgeId, uint64_t> paths;
    void createSegmentPaths();

    // Average marker graph edge coverage for all segments.
    MemoryMapped::Vector<float> segmentCoverage;
    void computeSegmentCoverage();

    // Keep track of the segment and position each marker graph edge corresponds to.
    // For each marker graph edge, store in the marker graph edge table
    // the corresponding segment id and position in the path, if any.
    // Indexed by the edge id in the marker graph.
    // This is needed when computing assembly graph journeys.
    MemoryMapped::Vector< pair<uint64_t, uint32_t> > markerGraphEdgeTable;
    void computeMarkerGraphEdgeTable(size_t threadCount);
    void computeMarkerGraphEdgeTableThreadFunction(size_t threadId);



    // The marker graph journeys of all oriented reads.
    // Indexed by OrientedReadId::getValue().
    // This is only stored temporarily and used to compute assembly graph journeys.
    MemoryMapped::VectorOfVectors<MarkerGraphJourneyEntry, uint64_t> markerGraphJourneys;
    void computeMarkerGraphJourneys(size_t threadCount);
    void computeMarkerGraphJourneysPass1(size_t threadId);
    void computeMarkerGraphJourneysPass2(size_t threadId);
    void computeMarkerGraphJourneysPass12(uint64_t pass);
    void sortMarkerGraphJourneys(size_t threadId);



    // The assembly graph journeys of all oriented reads.
    // Indexed by OrientedReadId::getValue().
    MemoryMapped::VectorOfVectors<AssemblyGraphJourneyEntry, uint64_t> assemblyGraphJourneys;
    void computeAssemblyGraphJourneys();
    void computeAssemblyGraphJourney(
        const span<MarkerGraphJourneyEntry> markerGraphJourney,
        vector<AssemblyGraphJourneyEntry>& assemblyGraphJourney);

    // Store appearances of segments in assembly graph journeys.
    // For each segment, store pairs (orientedReadId, position in assembly graph journey).
    // Indexed by the segmentId.
    // For each segment, they are sorted.
    MemoryMapped::VectorOfVectors<pair<OrientedReadId, uint64_t>, uint64_t>
        assemblyGraphJourneyInfos;
    void computeAssemblyGraphJourneyInfos();



    // A transition is a sequence of two consecutive positions
    // in the assembly graph journey of an oriented reads.
    // In other words, it describes the transition of an oriented read
    // from a segment to the next segment it encounters.
    // Transitions are used to create edges in the AssemblyGraph (gfa links).
    // Indexed by the linkId. For each link, they are sorted.
    class Transition : public array<MarkerGraphJourneyEntry, 2> {
    public:
        Transition(const array<MarkerGraphJourneyEntry, 2>& x) : array<MarkerGraphJourneyEntry, 2>(x) {}
        Transition() {}
    };
    using SegmentPair = pair<uint64_t, uint64_t>;
    using Transitions = vector< pair<OrientedReadId, Transition> >;
    std::map<SegmentPair, Transitions> transitionMap;
    void findTransitions(std::map<SegmentPair, Transitions>& transitionMap);



    // The links.
    class Link {
    public:
        uint64_t segmentId0;
        uint64_t segmentId1;

        // Flag to indicate whether the two segments are adjacent.
        // This is set if the last marker graph vertex of segmentId0
        // is the same as the first marker graph vertex of segmentId1.
        // In that case the separation will be set to 0.
        // However, the separation is just an estimate, so it
        // could be 0 even when the segments are ot adjacent.
        bool segmentsAreAdjacent;

        // Estimated separation in markers.
        int32_t separation;


        Link(
            uint64_t segmentId0 = 0,
            uint64_t segmentId1 = 0) :
            segmentId0(segmentId0),
            segmentId1(segmentId1) {}
    };
    MemoryMapped::Vector<Link> links;
    void createLinks(
        const std::map<SegmentPair, Transitions>& transitionMap,
        uint64_t minCoverage);

    // The transitions for each link.
    // Indexed by linkId.
    MemoryMapped::VectorOfVectors< pair<OrientedReadId, Transition>, uint64_t> transitions;
    uint64_t linkCoverage(uint64_t linkId) const
    {
        return transitions.size(linkId);
    }

    // The links for each source or target segments.
    // Indexed by segment id.
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> linksBySource;
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> linksByTarget;
    void createConnectivity();


    // Flag back-segments.
    // This does not do a full blown search for locally strongly connected components.
    // A segment is marked as a back-segment if:
    // - It has only a single incoming link.
    // - It has a single outgoing link.
    // - The incoming and outgoing links both connect to/from the same segment.
    void flagBackSegments();
    MemoryMapped::Vector<bool> isBackSegment;



    // Get the children or parents of a given segment.
    // Only use links with at least a specified coverage.
    void getChildren(
        uint64_t segmentId,
        uint64_t minimumLinkCoverage,
        vector<uint64_t>&
        ) const;
    void getParents(
        uint64_t segmentId,
        uint64_t minimumLinkCoverage,
        vector<uint64_t>&
        ) const;
    void getChildrenOrParents(
        uint64_t segmentId,
        uint64_t direction, // 0=forward (children), 1=backward (parents).
        uint64_t minimumLinkCoverage,
        vector<uint64_t>&
        ) const;


    // Find descendants of a given segment, up to a given distance in the graph.
    void findDescendants(
        uint64_t segmentId,
        uint64_t maxDistance,
        vector<uint64_t>& segmentIds
        ) const;

    // BFS with given begin/end.
    // Does a BFS which starts at segmentIdA.
    // and ends when segmentIdB is encountered.
    // The BFS if forward if direction is 0
    // and backward if direction is 1.
    // Computes a vector of all the segments encountered,
    // excluding segmentIdA and segmentIdB,
    // in the order in which they are encountered in the BFS.
    void targetedBfs(
        uint64_t segmentIdA,
        uint64_t segmentIdB,
        uint64_t direction,
        vector<uint64_t>& segments
        ) const;

    void writeGfa(const string& baseName) const;

    // Find the distinct oriented reads that appear on the path
    // of a segment. Also return the average edge coverage for the path.
    double findOrientedReadsOnSegment(
        uint64_t segmentId,
        vector<OrientedReadId>&) const;



    // Get information about the oriented reads that appear on the
    // marker graph path of a segment.
    class SegmentOrientedReadInformation {
    public:

        // The oriented reads on this segment,
        // each storage with an average offset relative to the segment.
        class Info {
        public:
            OrientedReadId orientedReadId;

            // The average offset, in markers, between the
            // beginning of this oriented read and the
            // beginning of the segment.
            int32_t averageOffset;
        };
        vector<Info> infos;
    };
    void getOrientedReadsOnSegment(
        uint64_t segmentId,
        SegmentOrientedReadInformation&) const;

    // Oriented read information for each segment.
    // This is only stored when needed.
    vector<SegmentOrientedReadInformation> segmentOrientedReadInformation;
    void storeSegmentOrientedReadInformation(size_t threadCount);
    void storeSegmentOrientedReadInformationThreadFunction(size_t threadId);



    // Estimate the offset between two segments.
    // Takes as input SegmentOrientedReadInformation objects
    // for the two segments.
    // Common oriented reads between the two segments are used
    // to estimate the average offset, in markers,
    // between the beginning of the segments.
    // The number of common oriented reads
    // is computed and stored in the last argument.
    // If that is zero, the computed offset is not valid.
    void estimateOffset(
        const SegmentOrientedReadInformation& info0,
        const SegmentOrientedReadInformation& info1,
        int64_t& offset,
        uint64_t& commonOrientedReadCount
        ) const;



    // Analyze a pair of segments for common oriented reads,
    // offsets, missing reads, etc.
    class SegmentPairInformation {
    public:

        // The total number of oriented reads present in each segment.
        array<uint64_t, 2> totalCount = {0, 0};

        // The number of oriented reads present in both segments.
        // If this is zero, the rest of the information is not valid.
        uint64_t commonCount = 0;

        // The offset of segment 1 relative to segment 0, in markers.
        int64_t offset = invalid<int64_t>;

        // The number of oriented reads present in each segment
        // but missing from the other segment,
        // and which should have been present based on the above estimated offset.
        array<uint64_t, 2> unexplainedCount = {0, 0};

        // The number of oriented reads that appear in only one
        // of the two segments, but based on the estimated offset
        // are too short to appear in the other segment.
        array<uint64_t, 2> shortCount = {0, 0};

        // Check that the above counts are consistent.
        void check() const
        {
            for(uint64_t i=0; i<2; i++) {
                SHASTA_ASSERT(commonCount + unexplainedCount[i] + shortCount[i] ==
                    totalCount[i]);
            }
        }

        // This computes the fraction of unexplained oriented reads,
        // without counting the short ones.
        double unexplainedFraction(uint64_t i) const
        {
            // return double(unexplainedCount[i]) / double(totalCount[i]);
            return double(unexplainedCount[i]) / double(commonCount + unexplainedCount[i]);
        }
        double maximumUnexplainedFraction() const
        {
            return max(unexplainedFraction(0), unexplainedFraction(1));
        }

        // Jaccard similarity, without counting the short reads.
        double jaccard() const
        {
            return double(commonCount) / double(commonCount + unexplainedCount[0] + unexplainedCount[1]);
        }
    };
    void analyzeSegmentPair(
        uint64_t segmentId0,
        uint64_t segmentId1,
        const SegmentOrientedReadInformation& info0,
        const SegmentOrientedReadInformation& info1,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        SegmentPairInformation&
        ) const;

    // Count the number of common oriented reads between a segment and a link,
    // without counting oriented reads that appear more than once on the
    // segment or on the link.
    void analyzeSegmentLinkPair(
        uint64_t segmentId,
        uint64_t linkId,
        uint64_t& commonOrientedReadCount
    ) const;



    // Find segment pairs a sufficient number of common reads
    // and with low unexplained fraction (in both directions)
    // between segmentId0 and one of its descendants within the specified distance.
    // This requires the vector segmentOrientedReadInformation above to be
    // available.
    void findSegmentPairs(
        uint64_t segmentId0,
        uint64_t maxDistance,
        uint64_t minCommonReadCount,
        double maxUnexplainedFraction,
        vector<uint64_t>& segmentIds1
    ) const;



    // Cluster the segments based on read composition.
    // We find segment pairs a sufficient number of common reads
    // and with low unexplained fraction (in both directions).
    void clusterSegments(size_t threadCount, uint64_t minClusterSize);
    class ClusterSegmentsData {
    public:

        // The segment pairs found by each thread.
        // In each pair, the lower number segment comes first.
        vector< vector< pair<uint64_t, uint64_t> > > threadPairs;
    };
    ClusterSegmentsData clusterSegmentsData;
    void clusterSegmentsThreadFunction1(size_t threadId);
    void addClusterPairs(size_t threadId, uint64_t segmentId0);
    MemoryMapped::Vector<uint64_t> clusterIds;



    // Analyze a subgraph of the assembly graph.

    // Classes used in analyzeSubgraph.
    class AnalyzeSubgraphClasses {
    public:

        // A JourneySnippet describes a sequence of consecutive positions
        // of the assembly graph journey of an oriented read.
        // An OrientedReadId can have than more one JourneySnippet in a given subgraph,
        // but this is not common. It can happen if the assembly graph contains a cycle.
        class JourneySnippet {
        public:

            // The OrientedReadId this refers to.
            OrientedReadId orientedReadId;

            // The sequence of segments encountered.
            vector<uint64_t> segmentIds;

            // The first and last position of this snippet
            // in the assembly graph journey of this OrientedReadId.
            uint64_t firstPosition;
            uint64_t lastPosition() const
            {
                return firstPosition + segmentIds.size() - 1;
            }
        };

        // A Cluster is a set of JourneySnippet's.
        class Cluster {
        public:

            // The snippets in this cluster.
            vector<JourneySnippet> snippets;
            uint64_t coverage() const
            {
                return snippets.size();
            }

            // The segments visited by the snippets of this cluster,
            // each stored with its coverage (number of snippets);
            vector< pair<uint64_t, uint64_t > > segments;
            vector<uint64_t> getSegments() const;

            // Remove segments with coverage less than the specified value.
            void cleanupSegments(uint64_t minClusterCoverage);

            // Construct the segments given the snippets.
            void constructSegments();
        };



        // The SnippetGraph is used by analyzeSubgraph2.
        // A vertex represents a set of snippets and stores
        // the corresponding snippet indexes.
        // An edge x->y is created if there is at least one snippet in y
        // that is an approximate subset of a snippet in x.
        // Strongly connected components are condensed, so after that
        //the graph is guaranteed to have no cycles.
        class SnippetGraphVertex {
        public:
            vector<uint64_t> snippetIndexes;
            uint64_t clusterId = invalid<uint64_t>;
            SnippetGraphVertex() {}
            SnippetGraphVertex(uint64_t snippetIndex) :
                snippetIndexes(1, snippetIndex) {}
        };
        using SnippetGraphBaseClass =
            boost::adjacency_list<boost::setS, boost::listS, boost::bidirectionalS, SnippetGraphVertex>;
        class SnippetGraph : public SnippetGraphBaseClass {
        public:
            uint64_t clusterCount = 0;
            void findDescendants(const vertex_descriptor, vector<vertex_descriptor>&) const;
            void writeGraphviz(const string& fileName) const;
        };
    };



    void analyzeSubgraph(
        const vector<uint64_t>& segmentIds,
        vector<AnalyzeSubgraphClasses::Cluster>&,
        bool debug) const;
    template<uint64_t N> void analyzeSubgraphTemplate(
        const vector<uint64_t>& segmentIds,
        vector<AnalyzeSubgraphClasses::Cluster>&,
        bool debug) const;

    // Given a segment, use a BFS to move in the specified direction until
    // we find a segment with sufficiently high Jaccard similarity
    // and number of common reads.
    // This returns invalid<uint64_t> if no such segment is found
    // within the specified distance.
    uint64_t findSimilarSegmentBfs(
        uint64_t segmentId,
        uint64_t direction, // 0 = forward, 1 = backward
        uint64_t maxDistance,
        uint64_t minCommon,
        double minJaccard) const;

    // Given a segment, move in the specified direction,
    // in order of increasing distance in markers, until
    // we find a segment with sufficiently high Jaccard similarity
    // and number of common reads.
    // This returns invalid<uint64_t> if no such segment is found
    // within the specified distance.
    uint64_t findSimilarSegment(
        uint64_t segmentId,
        uint64_t direction,     // 0 = forward, 1 = backward
        uint64_t maxDistance,   // In markers
        uint64_t minLinkCoverage,
        int32_t minLinkSeparation,
        uint64_t minCommon,
        double maxUnexplainedFraction,
        double minJaccard,
        vector<uint64_t>& segments) const;

    // Create an assembly path starting at a given segment.
    void createAssemblyPath(
        uint64_t segmentId,
        uint64_t direction,    // 0 = forward, 1 = backward
        AssemblyPath&
        ) const;
    void createAssemblyPath1(
        uint64_t segmentId,
        uint64_t direction,    // 0 = forward, 1 = backward
        vector<uint64_t>& path // The segmentId's of the path.
        ) const;
    void createAssemblyPath2(
        uint64_t segmentId,
        uint64_t direction,    // 0 = forward, 1 = backward
        vector<uint64_t>& path // The segmentId's of the path.
        ) const;
    void createAssemblyPath3(
        uint64_t segmentId,
        uint64_t direction,    // 0 = forward, 1 = backward
        AssemblyPath&
        ) const;



    // Compute link separation given a set of Transitions.
    template<class Container> static double linkSeparation(
        const Container& transitions,
        uint64_t pathLength0)
    {
        double averageLinkSeparation = 0.;

        for(const pair<OrientedReadId, Transition>& p: transitions) {
            const Transition& transition = p.second;
            const MarkerGraphJourneyEntry& entry0 = transition[0];
            const MarkerGraphJourneyEntry& entry1 = transition[1];

            SHASTA_ASSERT(entry1.ordinals[0] >= entry0.ordinals[1]);

            const int64_t linkSeparation =
                int64_t(entry1.ordinals[0] - entry0.ordinals[1]) -
                int64_t(pathLength0 - 1 - entry0.position) -
                int64_t(entry1.position);
            averageLinkSeparation += double(linkSeparation);
        }
        averageLinkSeparation /= double(transitions.size());

        return averageLinkSeparation;
    }
};




#endif

