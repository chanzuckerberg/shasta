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
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "array.hpp"
#include "unordered_map"
#include "vector.hpp"



namespace shasta {
    namespace mode3 {
        class AssemblyGraph;
    }

    // Some forward declarations of classes in the shasta namespace.
    class CompressedMarker;
    class MarkerGraph;
}



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
    // This is needed when computing pseudopaths.
    MemoryMapped::Vector< pair<uint64_t, uint32_t> > markerGraphEdgeTable;
    void computeMarkerGraphEdgeTable(size_t threadCount);
    void computeMarkerGraphEdgeTableThreadFunction(size_t threadId);



    // Pseudopaths for all oriented reads.
    // The pseudopath of an oriented read is the sequence of
    // marker graph edges it encounters.
    // For each entry in the pseudopath, the marker graph
    // edge is identified by the segmentId in the AssemblyGraph
    // and the edge position in the segment (e. g. the first marker graph edge
    // in the segment is at position 0).
    // Indexed by OrientedReadId::getValue().
    // It gets removed when no longer needed.
    class PseudoPathEntry {
    public:
        uint64_t segmentId;
        uint32_t position;
        array<uint32_t, 2> ordinals;

        bool operator<(const PseudoPathEntry& that) const
        {
            return ordinals[0] < that.ordinals[0];
        }
    };
    MemoryMapped::VectorOfVectors<PseudoPathEntry, uint64_t> pseudoPaths;
    void computePseudoPaths(size_t threadCount);
    void computePseudoPathsPass1(size_t threadId);
    void computePseudoPathsPass2(size_t threadId);
    void computePseudoPathsPass12(uint64_t pass);
    void sortPseudoPaths(size_t threadId);



    // The compressed pseudopath of an oriented read
    // is the sequence of segmentIds it encounters.
    // Note a segmentId can appear more than once on the compressed
    // pseudopath of an oriented read.
    class CompressedPseudoPathEntry {
    public:
        uint64_t segmentId;

        // The first and last PseudoPathEntry's that contributed to this
        // CompressedPseudoPathEntry.
        array<PseudoPathEntry, 2> pseudoPathEntries;
    };
    // Indexed by OrientedReadId::getValue().
    MemoryMapped::VectorOfVectors<CompressedPseudoPathEntry, uint64_t> compressedPseudoPaths;
    void computeCompressedPseudoPaths();
    void computeCompressedPseudoPath(
        const span<PseudoPathEntry> pseudoPath,
        vector<CompressedPseudoPathEntry>& compressedPseudoPath);



    // A pseudopath transition occurs when the pseudopath of an oriented read
    // moves from a segment to a different segment.
    // Transitions are used to create edges in the AssemblyGraph (gfa links).
    class Transition : public array<PseudoPathEntry, 2> {
    public:
        Transition(const array<PseudoPathEntry, 2>& x) : array<PseudoPathEntry, 2>(x) {}
        Transition() {}
    };

    // Find pseudopath transitions and store them keyed by the pair of segments.
    using SegmentPair = pair<uint64_t, uint64_t>;
    using Transitions = vector< pair<OrientedReadId, Transition> >;
    std::map<SegmentPair, Transitions> transitionMap;
    void findTransitions(std::map<SegmentPair, Transitions>& transitionMap);



    // The links.
    class Link {
    public:
        uint64_t segmentId0;
        uint64_t segmentId1;

        Link(
            uint64_t segmentId0 = 0,
            uint64_t segmentId1 = 0,
            uint64_t coverage = 0) :
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

    // Find descendants of a given segment, up to a given distance in the graph.
    void findDescendants(
        uint64_t segmentId,
        uint64_t maxDistance,
        vector<uint64_t>& segmentIds
        ) const;

    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;

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
        int64_t offset = std::numeric_limits<int64_t>::max();


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
            return double(unexplainedCount[i]) / double(totalCount[i]);
        }
        double maximumUnexplainedFraction() const
        {
            return max(unexplainedFraction(0), unexplainedFraction(1));
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



    // Compute link separation given a set of Transitions.
    template<class Container> static double linkSeparation(
        const Container& transitions,
        uint64_t pathLength0)
    {
        double averageLinkSeparation = 0.;

        for(const pair<OrientedReadId, Transition>& p: transitions) {
            const Transition& transition = p.second;
            const auto& pseudoPathEntry0 = transition[0];
            const auto& pseudoPathEntry1 = transition[1];

            SHASTA_ASSERT(pseudoPathEntry1.ordinals[0] >= pseudoPathEntry0.ordinals[1]);

            const int64_t linkSeparation =
                int64_t(pseudoPathEntry1.ordinals[0] - pseudoPathEntry0.ordinals[1]) -
                int64_t(pathLength0 - 1 - pseudoPathEntry0.position) -
                int64_t(pseudoPathEntry1.position);
            averageLinkSeparation += double(linkSeparation);
        }
        averageLinkSeparation /= double(transitions.size());

        return averageLinkSeparation;
    }
};




#endif

