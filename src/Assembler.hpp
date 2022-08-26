#ifndef SHASTA_ASSEMBLER_HPP
#define SHASTA_ASSEMBLER_HPP

// Shasta.
#include "Alignment.hpp"
#include "AlignmentCandidates.hpp"
#include "AssemblyGraph2Statistics.hpp"
#include "HttpServer.hpp"
#include "Kmer.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "MemoryMappedObject.hpp"
#include "MultithreadedObject.hpp"
#include "ReadGraph.hpp"
#include "ReadId.hpp"
#include "shastaTypes.hpp"

// Standard library.
#include "memory.hpp"
#include "string.hpp"
#include "tuple.hpp"
#include "utility.hpp"

namespace shasta {

    class Assembler;
    class AssemblerInfo;
    class AssemblyGraph;
    class Alignment;
    class AlignmentData;
    class AlignmentGraph;
    class AlignmentInfo;
    class AlignOptions;
    class AssemblerOptions;
    class AssembledSegment;
    class AssemblyGraph2;
    class CompressedAssemblyGraph;
    class ConsensusCaller;
    class Histogram2;
    class InducedAlignment;
    class LocalAssemblyGraph;
    class LocalAlignmentCandidateGraph;
    class LocalAlignmentGraph;
    class LocalMarkerGraph;
    class LocalReadGraph;
    class LocalReadGraphTriangles;
    class LocalMarkerGraphRequestParameters;
    class LongBaseSequences;
    class MarkerConnectivityGraph;
    class MarkerConnectivityGraphVertexMap;
    class Mode2AssemblyOptions;
    class OrientedReadPair;
    class Reads;
    class ReferenceOverlapMap;


    namespace MemoryMapped {
        class ByteAllocator;
    }

    // Write an html form to select strand.
    void writeStrandSelection(
        ostream&,               // The html stream to write the form to.
        const string& name,     // The selection name.
        bool select0,           // Whether strand 0 is selected.
        bool select1);          // Whether strand 1 is selected.

    namespace Align4 {
        class MatrixEntry;
        class Options;
    }

    namespace mode3 {
        class AssemblyGraph;
    }

    extern template class MultithreadedObject<Assembler>;
}

namespace spoa {
    class AlignmentEngine;
    class Graph;
}

class DisjointSets;


// Sanity check that we are compiling on x86_64 or aarch64
#if !__x86_64__ && !__aarch64__
#error "Shasta can only be built on an x86_64 machine (64-bit Intel/AMD) or an ARM64 machine. "
#endif



// Class used to store various pieces of assembler information in shared memory.
class shasta::AssemblerInfo {
public:

    // The read representation used: 0 = raw sequence, 1 = RLE sequence
    uint64_t readRepresentation;

    // The length of k-mers used to define markers.
    size_t k;

    // The page size in use for this run.
    size_t largeDataPageSize;

    // Assembly mode (0=haploid, 2=phased).
    uint64_t assemblyMode;

    // Read graph creation method.
    uint64_t readGraphCreationMethod;


    // Statistics on the number of reads discarded on input.
    // These are incremented during each call to addReadsFromFasta.

    // The number of reads and raw bases discarded because the read
    // contained invalid bases.
    uint64_t discardedInvalidBaseReadCount = 0;
    uint64_t discardedInvalidBaseBaseCount = 0; // Only counts the valid bases in those reads.

    // The number of reads and raw bases discarded because the read length
    // was less than minReadLength.
    uint64_t discardedShortReadReadCount = 0;
    uint64_t discardedShortReadBaseCount = 0;

    // The number of reads and raw bases discarded because the read
    // contained repeat counts greater than 255.
    uint64_t discardedBadRepeatCountReadCount = 0;
    uint64_t discardedBadRepeatCountBaseCount = 0;



    // Statistics for the reads kept in the assembly
    // and not discarded on input.
    // These are computed and stored by histogramReadLength.
    size_t readCount = 0;
    size_t baseCount = 0;
    size_t readN50 = 0;
    uint64_t minReadLength = 0;

    // Other read statistics.
    size_t palindromicReadCount = 0;
    size_t chimericReadCount = 0;
    uint64_t isolatedReadCount = 0;
    uint64_t isolatedReadBaseCount = 0;

    // Alignment criteria actually used.
    // For readGraph creation method 0, they are the values specified
    // by the command line options and/or configuration.
    // For readGraph creation method 2, they are dynamically selected
    // based on alignments statistics.
    double actualMinAlignedFraction = 0;
    uint64_t actualMinAlignedMarkerCount = 0;
    uint64_t actualMaxDrift = 0;
    uint64_t actualMaxSkip = 0;
    uint64_t actualMaxTrim = 0;

    // Marker graph statistics.
    size_t markerGraphVerticesNotIsolatedCount = 0;
    size_t markerGraphEdgesNotRemovedCount = 0;
    uint64_t markerGraphMinCoverageUsed = 0;

    // Assembly graph statistics.
    size_t assemblyGraphAssembledEdgeCount = 0;
    size_t totalAssembledSegmentLength = 0;
    size_t longestAssembledSegmentLength = 0;
    size_t assembledSegmentN50 = 0;

    // Mode 2 assembly statistics.
    AssemblyGraph2Statistics assemblyGraph2Statistics;

    // Performance information.
    double assemblyElapsedTimeSeconds = 0.;
    double averageCpuUtilization;
    uint64_t peakMemoryUsage = 0ULL;
    uint64_t threadCount = 0;
    uint64_t virtualCpuCount = 0;
    uint64_t totalAvailableMemory = 0;

    inline string peakMemoryUsageForSummaryStats() {
        return peakMemoryUsage > 0 ? to_string(peakMemoryUsage) : "Not determined.";
    }
};



class shasta::Assembler :
    public MultithreadedObject<Assembler>
    , public HttpServer {
public:


    /***************************************************************************

    The constructors specify the file name prefix for binary data files.
    If this is a directory name, it must include the final "/".

    The constructor also specifies the page size for binary data files.
    Typically, for a large run binary data files will reside in a huge page
    file system backed by 2MB pages.
    1GB huge pages are also supported.
    The page sizes specified here must be equal to, or be an exact multiple of,
    the actual size of the pages backing the data.

    ***************************************************************************/

    // Constructor to be called one to create a new run.
    Assembler(
        const string& largeDataFileNamePrefix,
        bool createNew,
        uint64_t readRepresentation, // 0 = raw sequence, 1 = RLE sequence. Only used if createNew.
        size_t largeDataPageSize);

    // Add reads.
    // The reads in the specified file are added to those already previously present.
    void addReads(
        const string& fileName,
        uint64_t minReadLength,
        bool noCache,
        size_t threadCount);

    // Create a histogram of read lengths.
    void histogramReadLength(const string& fileName);


    // Functions related to markers.
    // See the beginning of Marker.hpp for more information.
    void findMarkers(size_t threadCount);
    void accessMarkers();
    void writeMarkers(ReadId, Strand, const string& fileName);
    vector<KmerId> getMarkers(ReadId, Strand);
    void writeMarkerFrequency();

    // Write the reads that overlap a given read.
    void writeOverlappingReads(ReadId, Strand, const string& fileName);

    // Compute a marker alignment of two oriented reads.
    void alignOrientedReads(
        ReadId, Strand,
        ReadId, Strand,
        size_t maxSkip,  // Maximum ordinal skip allowed.
        size_t maxDrift, // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency
    );

    // Compute marker alignments of an oriented read with all reads
    // for which we have an Overlap.
    void alignOverlappingOrientedReads(
        ReadId, Strand,
        size_t maxSkip,   // Maximum ordinal skip allowed.
        size_t maxDrift,  // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency,
        size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
        size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
    );



    // Compute an alignment for each alignment candidate.
    // Store summary information for the ones that are good enough,
    // without storing details of the alignment.
    void computeAlignments(
        const AlignOptions&,

        // Number of threads. If zero, a number of threads equal to
        // the number of virtual processors is used.
        size_t threadCount
    );
    void accessAlignmentData();
    void accessAlignmentDataReadWrite();
    void writeAlignmentDetails() const;


    // Loop over all alignments in the read graph
    // to create vertices of the global marker graph.
    // Throw away vertices with coverage (number of markers)
    // less than minCoverage or more than maxCoverage.
    // Also throw away "bad" vertices - that is, vertices
    // with more than one marker on the same oriented read.
    void createMarkerGraphVertices(

        // Minimum coverage (number of markers) for a vertex
        // of the marker graph to be kept.
        size_t minCoverage,

        // Maximum coverage (number of markers) for a vertex
        // of the marker graph to be kept.
        size_t maxCoverage,

        // Minimum coverage per strand (number of markers required
        // on each strand) for a vertex of the marker graph to be kept.
        uint64_t minCoveragePerStrand,

        // Flag that specifies whether to allow more than one marker on the
        // same oriented read id on a single marker graph vertex.
        bool allowDuplicateMarkers,

        // These two are used by PeakFinder in the automatic selection
        // of minCoverage when minCoverage is set to 0.
        double peakFinderMinAreaFraction,
        uint64_t peakFinderAreaStartIndex,

        // Number of threads. If zero, a number of threads equal to
        // the number of virtual processors is used.
        size_t threadCount
    );



    // Find the vertex of the global marker graph that contains a given marker.
    // The marker is specified by the ReadId and Strand of the oriented read
    // it belongs to, plus the ordinal of the marker in the oriented read.
    MarkerGraphVertexId getGlobalMarkerGraphVertex(
        ReadId,
        Strand,
        uint32_t ordinal) const;

    // Find the markers contained in a given vertex of the global marker graph.
    // Returns the markers as tuples(read id, strand, ordinal).
    vector< tuple<ReadId, Strand, uint32_t> >
        getGlobalMarkerGraphVertexMarkers(MarkerGraph::VertexId) const;

    // Compute marker graph vertex coverage statistics by KmerId.
    void vertexCoverageStatisticsByKmerId() const;



    // Approximate transitive reduction of the marker graph.
    // This does the following, in this order:
    // - All edges with coverage less than or equal to lowCoverageThreshold
    //   are marked wasRemovedByTransitiveReduction.
    // - All edges with coverage 1 and a marker skip
    //   greater than edgeMarkerSkipThreshold
    //   are marked wasRemovedByTransitiveReduction.
    // - Edges with coverage greater than lowCoverageThreshold
    //   and less then highCoverageThreshold are processed in
    //   ordered of increasing coverage:
    //   * For each such edge A->B, we look for a path of length
    //     at most maxDistance starting at A and ending at B  that does not use
    //     edge A->B and also does not use any
    //     edges already marked wasRemovedByTransitiveReduction.
    //   * If such a path is found, the edge is marked
    //     wasRemovedByTransitiveReduction.
    // - Edges with coverage highCoverageThreshold or greater
    //   are left untouched.
    // The marker graph is guaranteed to be strand symmetric
    // when this begins, and we have to guarantee that it remains
    // strand symmetric when this ends.
    // To achieve this, we always process the two edges
    // in a reverse complemented pair together.
    void transitiveReduction(
        size_t lowCoverageThreshold,
        size_t highCoverageThreshold,
        size_t maxDistance,
        size_t edgeMarkerSkipThreshold);



    // Approximate reverse transitive reduction of the marker graph.
    // The goal is to remove local back-edges.
    // This works similarly to transitive reduction,
    // but in the opposite direction.
    // This does the following:
    // - Edges with coverage greater than lowCoverageThreshold
    //   and less then highCoverageThreshold are processed in
    //   ordered of increasing coverage:
    //   * For each such edge A->B, we look for a path of length
    //     at most maxDistance starting at B and ending at A
    //     that does not use edge A->B and also does not use any
    //     edges already marked wasRemovedByTransitiveReduction.
    //   * If such a path is found, the edge is marked
    //     wasRemovedByTransitiveReduction.
    void reverseTransitiveReduction(
        size_t lowCoverageThreshold,
        size_t highCoverageThreshold,
        size_t maxDistance);



private:

    // Data filled in by the constructor.
    string largeDataFileNamePrefix;
    size_t largeDataPageSize;

    // Function to construct names for binary objects.
    string largeDataName(const string& name) const
    {
        if(largeDataFileNamePrefix.empty()) {
            return "";  // Anonymous;
        } else {
            return largeDataFileNamePrefix + name;
        }
    }



    // Various pieces of assembler information stored in shared memory.
    // See class AssemblerInfo for more information.
public:
    MemoryMapped::Object<AssemblerInfo> assemblerInfo;
    uint64_t getMarkerGraphMinCoverageUsed() const
    {
        return assemblerInfo->markerGraphMinCoverageUsed;
    }
private:



    // Reads in RLE representation.
    shared_ptr<Reads> reads;
public:
    const Reads& getReads() const {
        SHASTA_ASSERT(reads);
        return *reads;
    }

    uint64_t adjustCoverageAndGetNewMinReadLength(uint64_t desiredCoverage);

    // Write a csv file with summary information for each read.
public:
    void writeReadsSummary();

    void computeReadIdsSortedByName();


private:



    // Table of all k-mers of length k.
    // Among all 4^k k-mers of length k, we choose a subset
    // that we call "markers".
    // The value of k used is stored in assemblerInfo.
    // The k-mer table is a vector of 4^k pairs,
    // indexed by k-mer id as computed using Kmer::id(k).
    // The markers are selected at the beginning of an assembly
    // and never changed, and selected in such a way that,
    // if (and only if) a k-mer is a marker, its reverse complement
    // is also a marker. That is, for all permitted values of i, 0 <= i < 4^k:
    // kmerTable[i].isMarker == kmerTable[kmerTable[i].reverseComplementKmerId].isMarker
    MemoryMapped::Vector<KmerInfo> kmerTable;
    void checkKmersAreOpen() const;

public:
    void accessKmers();
    void writeKmers(const string& fileName) const;

    // Select marker k-mers randomly.
    void randomlySelectKmers(
        size_t k,           // k-mer length.
        double probability, // The probability that a k-mer is selected as a marker.
        int seed            // For random number generator.
    );



    // Select marker k-mers randomly, but excluding
    // the ones that have high frequency in the reads.
    void selectKmersBasedOnFrequency(

        // k-mer length.
        size_t k,

        // The desired marker density
        double markerDensity,

        // Seed for random number generator.
        int seed,

        // Exclude k-mers enriched by more than this amount.
        // Enrichment is the ratio of k-mer frequency in reads
        // over what a random distribution would give.
        double enrichmentThreshold,

        size_t threadCount
    );



    // In this version, marker k-mers are selected randomly, but excluding
    // any k-mer that is over-enriched even in a single oriented read.
    void selectKmers2(

        // k-mer length.
        size_t k,

        // The desired marker density
        double markerDensity,

        // Seed for random number generator.
        int seed,

        // Exclude k-mers enriched by more than this amount,
        // even in a single oriented read.
        // Enrichment is the ratio of k-mer frequency in reads
        // over what a random distribution would give.
        double enrichmentThreshold,

        size_t threadCount
    );
private:

    class SelectKmers2Data {
    public:

        double enrichmentThreshold;

        // The number of times each k-mer appears in an oriented read.
        // Indexed by KmerId.
        MemoryMapped::Vector<uint64_t> globalFrequency;

        // The number of oriented reads that each k-mer is
        // over-enriched in by more than a factor enrichmentThreshold.
        // Indexed by KmerId.
        MemoryMapped::Vector<ReadId> overenrichedReadCount;

    };
    SelectKmers2Data selectKmers2Data;
    void selectKmers2ThreadFunction(size_t threadId);



    // In this version, marker k-mers are selected randomly, but excluding
    // k-mers that appear repeated at short distances in any oriented read.
    // More precisely, for each k-mer we compute the minimum distance
    // (in RLE bases) at which any two copies of that k-mer appear in any oriented read.
    // K-mers for which this minimum distance is less than distanceThreshold
    // are not used as markers. Marker k-mers are selected randomly among the
    // remaining k-mers, until the desired marker density is achieved.
public:
    void selectKmers4(

        // k-mer length.
        uint64_t k,

        // The desired marker density
        double markerDensity,

        // Seed for random number generator.
        uint64_t seed,

        // Exclude k-mers that appear in any read in two copies,
        // with the two copies closer than this distance (in RLE bases).
        uint64_t distanceThreshold,

        size_t threadCount
    );
private:
    void selectKmers4ThreadFunction(size_t threadId);
    class SelectKmers4Data {
    public:

        // The number of times each k-mer appears in an oriented read.
        // Indexed by KmerId.
        MemoryMapped::Vector<uint64_t> globalFrequency;

        // The minimum distance at which two copies of each k-mer
        // appear in any oriented read.
        // Indexed by KmerId.
        MemoryMapped::Vector< pair<std::mutex, uint32_t> > minimumDistance;

    };
    SelectKmers4Data selectKmers4Data;



    // Read the k-mers from file.
public:
    void readKmersFromFile(uint64_t k, const string& fileName);

private:
    void computeKmerFrequency(size_t threadId);
    void initializeKmerTable();



    // The markers on all oriented reads. Indexed by OrientedReadId::getValue().
    MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t> markers;
    void checkMarkersAreOpen() const;

    // Get markers sorted by KmerId for a given OrientedReadId.
    void getMarkersSortedByKmerId(
        OrientedReadId,
        vector<MarkerWithOrdinal>&) const;

    // Given a marker by its OrientedReadId and ordinal,
    // return the corresponding global marker id.
    MarkerId getMarkerId(OrientedReadId, uint32_t ordinal) const;
    MarkerId getReverseComplementMarkerId(OrientedReadId, uint32_t ordinal) const;
    MarkerId getMarkerId(const MarkerDescriptor& m) const
    {
        return getMarkerId(m.first, m.second);
    }
    MarkerId getReverseComplementMarkerId(const MarkerDescriptor& m) const
    {
        return getReverseComplementMarkerId(m.first, m.second);
    }

    // Inverse of the above: given a global marker id,
    // return its OrientedReadId and ordinal.
    // This requires a binary search in the markers toc.
    // This could be avoided, at the cost of storing
    // an additional 4 bytes per marker.
public:
    pair<OrientedReadId, uint32_t> findMarkerId(MarkerId) const;
private:



    // Pairs (KmerId, ordinal), sorted by KmerId, for each oriented read.
    // Indexed by orientedReadId.getValue().
    // Used by alignment method 4.
    MemoryMapped::VectorOfVectors< pair<KmerId, uint32_t>, uint64_t> sortedMarkers;
public:
    void computeSortedMarkers(uint64_t threadCount);
    bool accessSortedMarkers();
private:
    void computeSortedMarkersThreadFunction1(size_t threadId);
    void computeSortedMarkersThreadFunction2(size_t threadId);



    // Given a MarkerId, compute the MarkerId of the
    // reverse complemented marker.
    MarkerId findReverseComplement(MarkerId) const;


    // Flag palindromic reads.
public:
    void flagPalindromicReads(
        uint32_t maxSkip,
        uint32_t maxDrift,
        uint32_t maxMarkerFrequency,
        double alignedFractionThreshold,
        double nearDiagonalFractionThreshold,
        uint32_t deltaThreshold,
        size_t threadCount);
private:
    void flagPalindromicReadsThreadFunction(size_t threadId);
    class FlagPalindromicReadsData {
    public:
        uint32_t maxSkip;
        uint32_t maxDrift;
        uint32_t maxMarkerFrequency;
        double alignedFractionThreshold;
        double nearDiagonalFractionThreshold;
        uint32_t deltaThreshold;
    };
    FlagPalindromicReadsData flagPalindromicReadsData;


    // Check if an alignment between two reads should be suppressed,
    // bases on the setting of command line option
    // --Align.sameChannelReadAlignment.suppressDeltaThreshold.
    bool suppressAlignment(ReadId, ReadId, uint64_t delta);

    // Remove all alignment candidates for which suppressAlignment
    // returns false.
public:
    void suppressAlignmentCandidates(uint64_t delta, size_t threadCount);
private:
    class SuppressAlignmentCandidatesData {
    public:
        uint64_t delta;
        MemoryMapped::Vector<bool> suppress; // For each alignment candidate.
    };
    SuppressAlignmentCandidatesData suppressAlignmentCandidatesData;
    void suppressAlignmentCandidatesThreadFunction(size_t threadId);



    // Alignment candidates found by the LowHash algorithm.
    // They all have readId0<readId1.
    // They are interpreted with readId0 on strand 0.
    AlignmentCandidates alignmentCandidates;

public:
    void writeAlignmentCandidates(bool useReadName=false, bool verbose=false) const;
private:


    // Use the LowHash (modified MinHash) algorithm to find candidate alignments.
    // Use as features sequences of m consecutive special k-mers.
public:
    void findAlignmentCandidatesLowHash0(
        size_t m,                       // Number of consecutive k-mers that define a feature.
        double hashFraction,            // Low hash threshold.
        // Iteration control. See MinHashOptions for details.
        size_t minHashIterationCount,
        double alignmentCandidatesPerRead,

        size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for lowHash.
        size_t minBucketSize,           // The minimum size for a bucket to be used.
        size_t maxBucketSize,           // The maximum size for a bucket to be used.
        size_t minFrequency,            // Minimum number of lowHash hits for a pair to become a candidate.
        size_t threadCount
    );
    void findAlignmentCandidatesLowHash1(
        size_t m,                       // Number of consecutive k-mers that define a feature.
        double hashFraction,            // Low hash threshold.
        size_t minHashIterationCount,
        size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for lowHash.
        size_t minBucketSize,           // The minimum size for a bucket to be used.
        size_t maxBucketSize,           // The maximum size for a bucket to be used.
        size_t minFrequency,            // Minimum number of lowHash hits for a pair to become a candidate.
        size_t threadCount
    );
    void markAlignmentCandidatesAllPairs();
    void accessAlignmentCandidates();
    void accessAlignmentCandidateTable();
    vector<OrientedReadPair> getAlignmentCandidates() const;
    void computeCandidateTable();

private:
    void checkAlignmentCandidatesAreOpen() const;


    // LowHash statistics for read.
    // For each read we count the number of times a low hash
    // puts the read in:
    // - A sparse bucket (bucketSize < minFrequency), index 0.
    // - A good bucket (minFrequency <= bucketSize <= maxFrequency), index 1.
    // - A crowded bucket (bucketSize > maxFrequency), index 2.
    MemoryMapped::Vector< array<uint64_t, 3> > readLowHashStatistics;
    void accessReadLowHashStatistics();

    bool createLocalAlignmentCandidateGraph(
        vector<OrientedReadId>& starts,
        uint32_t maxDistance,               // How far to go from starting oriented read.
        bool allowChimericReads,
        double timeout,                     // Or 0 for no timeout.
        bool inGoodAlignmentsRequired,      // Only add an edge to the local graph if it's in the "good" alignments
        bool inReadgraphRequired,           // Only add an edge to the local graph if it's in the ReadGraph
        LocalAlignmentCandidateGraph& graph
    );

    // This method is used as an alternative to createLocalAlignmentCandidateGraph, in the case that the user
    // wants to see only the edges that are inferred from the PAF, and none others. Coloring/labelling w.r.t.
    // the different subgroups still applies (candidate, good alignment, read graph)
    bool createLocalReferenceGraph(
        vector<OrientedReadId>& starts,
        uint32_t maxDistance,           // How far to go from starting oriented read.
        bool allowChimericReads,
        double timeout,                 // Or 0 for no timeout.
        LocalAlignmentCandidateGraph& graph
    );
public:
    // Construct a subgraph around a read in the candidate graph
    // and write all subgraph reads to a Fasta file.
    void writeLocalAlignmentCandidateReads(
            ReadId readId,
            Strand strand,
            uint32_t maxDistance,
            bool allowChimericReads,
            bool allowCrossStrandEdges,
            bool allowInconsistentAlignmentEdges
    );
private:
    // Compute a marker alignment of two oriented reads.
    void alignOrientedReads(
        OrientedReadId,
        OrientedReadId,
        size_t maxSkip,     // Maximum ordinal skip allowed.
        size_t maxDrift,    // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency
    );
    // This lower level version takes as input vectors of
    // markers already sorted by kmerId.
    void alignOrientedReads(
        const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
        size_t maxSkip,  // Maximum ordinal skip allowed.
        size_t maxDrift, // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency
    );
    // This version allows reusing the AlignmentGraph and Alignment
    void alignOrientedReads(
        const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
        size_t maxSkip,             // Maximum ordinal skip allowed.
        size_t maxDrift,            // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency,
        bool debug,
        AlignmentGraph&,
        Alignment&,
        AlignmentInfo&
    );
public:
    void analyzeAlignmentMatrix(ReadId, Strand, ReadId, Strand);
private:


    // Alternative alignment functions with 1 suffix (SeqAn).
    void alignOrientedReads1(
        OrientedReadId,
        OrientedReadId,
        int matchScore,
        int mismatchScore,
        int gapScore);
    void alignOrientedReads1(
        OrientedReadId,
        OrientedReadId,
        int matchScore,
        int mismatchScore,
        int gapScore,
        Alignment&,
        AlignmentInfo&);
public:
    void alignOrientedReads1(
        ReadId, Strand,
        ReadId, Strand,
        int matchScore,
        int mismatchScore,
        int gapScore);
private:


    // Alternative alignment function with 3 suffix (SeqAn, banded).
    void alignOrientedReads3(
        OrientedReadId,
        OrientedReadId,
        int matchScore,
        int mismatchScore,
        int gapScore,
        double downsamplingFactor,
        int bandExtend,
        int maxBand,
        Alignment&,
        AlignmentInfo&);



    // Member functions that use alignment algorithm 4.
public:

    // Python-callable version.
    void alignOrientedReads4(
        ReadId, Strand,
        ReadId, Strand,
        uint64_t deltaX,
        uint64_t deltaY,
        uint64_t minEntryCountPerCell,
        uint64_t maxDistanceFromBoundary,
        uint64_t minAlignedMarkerCount,
        double minAlignedFraction,
        uint64_t maxSkip,
        uint64_t maxDrift,
        uint64_t maxTrim,
        uint64_t maxBand,
        int64_t matchScore,
        int64_t mismatchScore,
        int64_t gapScore) const;

    // Align two reads using alignment method 4.
    // If debug is true, detailed output to html is produced.
    // Otherwise, html is not used.
    void alignOrientedReads4(
        OrientedReadId,
        OrientedReadId,
        const Align4::Options&,
        MemoryMapped::ByteAllocator&,
        Alignment&,
        AlignmentInfo&,
        bool debug) const;

    // Intermediate level version used by the http server.
    void alignOrientedReads4(
        OrientedReadId,
        OrientedReadId,
        uint64_t deltaX,
        uint64_t deltaY,
        uint64_t minEntryCountPerCell,
        uint64_t maxDistanceFromBoundary,
        uint64_t minAlignedMarkerCount,
        double minAlignedFraction,
        uint64_t maxSkip,
        uint64_t maxDrift,
        uint64_t maxTrim,
        uint64_t maxBand,
        int64_t matchScore,
        int64_t mismatchScore,
        int64_t gapScore,
        Alignment&,
        AlignmentInfo&
        ) const;

private:


    // Create a local alignment graph starting from a given oriented read
    // and walking out a given distance on the global alignment graph.
    // An alignment graph is an undirected graph in which each vertex
    // represents an oriented read. Two vertices are joined by an
    // undirected edge if we have a found a good alignment,
    // and stored it, between the corresponding oriented reads.
    bool createLocalAlignmentGraph(
        OrientedReadId,
        size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
        size_t maxTrim,                 // Maximum left/right trim (expressed in bases).
        uint32_t distance,              // How far to go from starting oriented read.
        double timeout,                 // Or 0 for no timeout.
        LocalAlignmentGraph&
    );

    // Compute marker alignments of an oriented read with all reads
    // for which we have an Overlap.
    void alignOverlappingOrientedReads(
        OrientedReadId,
        size_t maxSkip,                 // Maximum ordinal skip allowed.
        size_t maxDrift,                // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency,
        size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
        size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
    );



    // Count the common marker near a given ordinal offset for
    // two oriented reads. This can be used to check
    // whether an alignmnent near the specified ordinal offset exists.
public:
    uint32_t countCommonMarkersNearOffset(
        OrientedReadId,
        OrientedReadId,
        int32_t offset,
        int32_t offsetTolerance
    );
    uint32_t countCommonMarkersWithOffsetIn(
        OrientedReadId,
        OrientedReadId,
        int32_t minOffset,
        int32_t maxOffset
    );
private:



    // The good alignments we found.
    // They are stored with readId0<readId1 and with strand0==0.
    // The order in compressedAlignments matches that in alignmentData.
    MemoryMapped::Vector<AlignmentData> alignmentData;
    MemoryMapped::VectorOfVectors<char, uint64_t> compressedAlignments;

    void checkAlignmentDataAreOpen() const;
public:
    void accessCompressedAlignments();
private:



    // Get the stored alignments involving a given oriented read.
    // This performs swaps and reverse complementing as necessary,
    // to return alignments in which the first oriented read is
    // the one specified as the argument.
    class StoredAlignmentInformation {
    public:
        uint64_t alignmentId;
        OrientedReadId orientedReadId;
        Alignment alignment;
    };
    void getStoredAlignments(
        OrientedReadId,
        vector<StoredAlignmentInformation>&) const;
    void getStoredAlignments(
        OrientedReadId,
        const vector<OrientedReadId>&,
        vector<StoredAlignmentInformation>&) const;



    // The alignment table stores the AlignmentData that each oriented read is involved in.
    // Stores, for each OrientedReadId, a vector of indexes into the alignmentData vector.
    // Indexed by OrientedReadId::getValue(),
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> alignmentTable;
    void computeAlignmentTable();



    // Private functions and data used by computeAlignments.
    void computeAlignmentsThreadFunction(size_t threadId);
    class ComputeAlignmentsData {
    public:

        // Not owned.
        const AlignOptions* alignOptions = 0;

        // The AlignmentInfo found by each thread.
        vector< vector<AlignmentData> > threadAlignmentData;

        // Compressed alignments corresponding to the AlignmentInfo found by each thread.
        vector< shared_ptr< MemoryMapped::VectorOfVectors<char, uint64_t> > > threadCompressedAlignments;
    };
    ComputeAlignmentsData computeAlignmentsData;



    // Find in the alignment table the alignments involving
    // a given oriented read, and return them with the correct
    // orientation (this may involve a swap and/or reverse complement
    // of the AlignmentInfo stored in the alignmentTable).
    vector< pair<OrientedReadId, AlignmentInfo> >
        findOrientedAlignments(OrientedReadId, bool inReadGraphOnly) const;



    // Analyze the stored alignments involving a given oriented read.
private:
    void analyzeAlignments1(ReadId, Strand) const;


    // Read graph and related functions and data.
    // For more information, see comments in ReadGraph.hpp.
    ReadGraph readGraph;
public:
    void createReadGraph(
        uint32_t maxAlignmentCount,
        uint32_t maxTrim);

    void createReadGraph2(
        uint32_t maxAlignmentCount,
        double markerCountPercentile,
        double alignedFractionPercentile,
        double maxSkipPercentile,
        double maxDriftPercentile,
        double maxTrimPercentile);

    void setReadGraph2Criteria(
            double markerCountPercentile,
            double alignedFractionPercentile,
            double maxSkipPercentile,
            double maxDriftPercentile,
            double maxTrimPercentile);

    bool passesReadGraph2Criteria(const AlignmentInfo& info) const;
    void accessReadGraph();
    void accessReadGraphReadWrite();
    void checkReadGraphIsOpen() const;
    void removeReadGraphBridges(uint64_t maxDistance);
    void analyzeReadGraph();
    void readGraphClustering();
    void writeReadGraphEdges(bool useReadName=false) const;



    // Triangle and least square analysis of the read graph
    // to flag inconsistent alignments.
    void flagInconsistentAlignments(
        uint64_t triangleErrorThreshold,
        uint64_t leastSquareErrorThreshold,
        uint64_t leastSquareMaxDistance,
        size_t threadCount);
private:
    void flagInconsistentAlignmentsThreadFunction1(size_t threadId);
    void flagInconsistentAlignmentsThreadFunction2(size_t threadId);
    class FlagInconsistentAlignmentsData {
    public:

        // Arguments of flagInconsistentAlignments, stored here
        // to make them visible to the threads.
        uint64_t triangleErrorThreshold;
        uint64_t leastSquareErrorThreshold;
        uint64_t leastSquareMaxDistance;

        // The alignment offset for each edge of the read graph,
        // oriented with the lowest OrientedReadId first.
        MemoryMapped::Vector<int32_t> edgeOffset;

        // The inconsistent read graph edge ids found by each thread.
        vector< vector<uint64_t> > threadEdgeIds;
    };
    FlagInconsistentAlignmentsData flagInconsistentAlignmentsData;
public:



    // Functions and data used with read creation for iterative assembly.
    void createReadGraphUsingPseudoPaths(
        int64_t matchScore,
        int64_t mismatchScore,
        int64_t gapScore,
        double mismatchSquareFactor,
        double minScore,
        uint64_t maxAlignmentCount,
        size_t threadCount);
    class CreateReadGraphsingPseudoPathsAlignmentData {
    public:
        uint64_t alignedMarkerCount = 0;

        // Pseudo-path alignment information.
        uint64_t weakMatchCount = 0;
        uint64_t strongMatchCount = 0;
        uint64_t mismatchCount = 0;
    };
    class CreateReadGraphUsingPseudoPathsData {
    public:
        int64_t matchScore;
        int64_t mismatchScore;
        int64_t gapScore;

        // The pseudopaths of all oriented reads.
        // Indexed by OrientedReadId::getValue().
        vector< vector<AssemblyGraphEdgeId> > pseudoPaths;

        // Vector to store information about each alignment.
        vector<CreateReadGraphsingPseudoPathsAlignmentData> alignmentInfos;
    };
    CreateReadGraphUsingPseudoPathsData createReadGraphUsingPseudoPathsData;

    // Thread function used to compute pseudoPaths.
    void createReadGraphUsingPseudoPathsThreadFunction1(size_t threadId);
    // Thread functions used to align pseudopaths.
    void createReadGraphUsingPseudoPathsThreadFunction2(size_t threadId);



#if 0
    // Functions and data for the version that uses mini-assemblies.
private:
    void createReadGraph2ThreadFunction(size_t threadId);
    void createReadGraph2LowLevel(ReadId);
    class CreateReadGraph2Data {
    public:
        vector<bool> keepAlignment;
    };
    CreateReadGraph2Data createReadGraph2Data;
public:
#endif



    // Approximate strand separation in the read graph.
    void flagCrossStrandReadGraphEdges1(int maxDistance, size_t threadCount);
private:
    void flagCrossStrandReadGraphEdges1ThreadFunction(size_t threadId);
    class FlagCrossStrandReadGraphEdges1Data {
    public:
        size_t maxDistance;
        vector<bool> isNearStrandJump;
    };
    FlagCrossStrandReadGraphEdges1Data flagCrossStrandReadGraphEdges1Data;
public:


    // Strict strand separation in the read graph.
    void flagCrossStrandReadGraphEdges2();



    // Create the ReadGraph given a bool vector that specifies which
    // alignments should be used in the read graph.
    void createReadGraphUsingSelectedAlignments(vector<bool>& keepAlignment);

    // Add alignments to avoid coverage holes.
    void fixCoverageHoles(vector<bool>& keepAlignment) const;


public:
    // Use the read graph to flag chimeric reads.
    void flagChimericReads(size_t maxDistance, size_t threadCount);
private:
    class FlagChimericReadsData {
    public:
        size_t maxDistance;
    };
    FlagChimericReadsData flagChimericReadsData;
    void flagChimericReadsThreadFunction(size_t threadId);


    // Create a local subgraph of the global read graph,
    // starting at a given vertex and extending out to a specified
    // distance (number of edges).
    bool createLocalReadGraph(
            OrientedReadId start,
            uint32_t maxDistance,   // How far to go from starting oriented read.
            bool allowChimericReads,
            bool allowCrossStrandEdges,
            bool allowInconsistentAlignmentEdges,
            double timeout,         // Or 0 for no timeout.
            LocalReadGraph&);

    // Create a local subgraph of the global read graph,
    // starting at any number of  given vertexes and extending out to a specified
    // distance (number of edges).
    bool createLocalReadGraph(
        const vector<OrientedReadId>& starts,
        uint32_t maxDistance,   // How far to go from starting oriented read.
        bool allowChimericReads,
        bool allowCrossStrandEdges,
        bool allowInconsistentAlignmentEdges,
        double timeout,         // Or 0 for no timeout.
        LocalReadGraph&);

    // Triangle analysis of the local read graph.
    // Returns a vector of triangles and their alignment residuals,
    // sorted by decreasing residual.
    void triangleAnalysis(
        LocalReadGraph&,
        LocalReadGraphTriangles&) const;

    // Singular value decomposition analysis of the local read graph.
    void leastSquareAnalysis(
        LocalReadGraph&,
        vector<double>& singularValues) const;

public:



    // Write a FASTA file containing all reads that appear in
    // the local read graph.
    void writeLocalReadGraphReads(
        ReadId,
        Strand,
        uint32_t maxDistance,
        bool allowChimericReads,
        bool allowCrossStrandEdges,
        bool allowInconsistentAlignmentEdges);


    // Compute connected components of the read graph.
    // This just writes a csv file and has no other side effects
    // (nothing is stored).
    void computeReadGraphConnectedComponents() const;



    // Private functions and data used by createMarkerGraphVertices.
private:
    void createMarkerGraphVerticesThreadFunction1(size_t threadId);
    void createMarkerGraphVerticesThreadFunction2(size_t threadId);
    void createMarkerGraphVerticesThreadFunction21(size_t threadId);
    void createMarkerGraphVerticesThreadFunction3(size_t threadId);
    void createMarkerGraphVerticesThreadFunction4(size_t threadId);
    void createMarkerGraphVerticesThreadFunction5(size_t threadId);
    void createMarkerGraphVerticesThreadFunction45(int);
    void createMarkerGraphVerticesThreadFunction6(size_t threadId);
    void createMarkerGraphVerticesThreadFunction7(size_t threadId);
    void createMarkerGraphVerticesDebug1(uint64_t stage);
    class CreateMarkerGraphVerticesData {
    public:

        // Parameters.
        uint64_t minCoveragePerStrand;
        bool allowDuplicateMarkers;

        // The total number of oriented markers.
        uint64_t orientedMarkerCount;

        // Disjoint sets data structures.
        shared_ptr<DisjointSets> disjointSetsPointer;

        // The disjoint set that each oriented marker was assigned to.
        // See createMarkerGraphVertices for details.
        MemoryMapped::Vector<MarkerGraph::VertexId> disjointSetTable;

        // Work area used for multiple purposes.
        // See createMarkerGraphVertices for details.
        MemoryMapped::Vector<MarkerGraph::VertexId> workArea;

        // The markers in each disjoint set with coverage in the requested range.
        MemoryMapped::VectorOfVectors<MarkerId, MarkerGraph::VertexId> disjointSetMarkers;

        // Flag disjoint sets that contain more than one marker on the same oriented read.
        MemoryMapped::Vector<bool> isBadDisjointSet;

    };
    CreateMarkerGraphVerticesData createMarkerGraphVerticesData;



    void checkMarkerGraphVerticesAreAvailable() const;

    // Check for consistency of globalMarkerGraphVertex and globalMarkerGraphVertices.
    void checkMarkerGraphVertices(
        size_t minCoverage,
        size_t maxCoverage);



    // Marker graph.
public:
    MarkerGraph markerGraph;
    void removeMarkerGraph()
    {
        markerGraph.remove();
    }

    // Find the reverse complement of each marker graph vertex.
    void findMarkerGraphReverseComplementVertices(size_t threadCount);
    void accessMarkerGraphVertices(bool readWriteAccess = false);
    void accessMarkerGraphReverseComplementVertex(bool readWriteAccess = false);
    void removeMarkerGraphVertices();
private:
    void findMarkerGraphReverseComplementVerticesThreadFunction1(size_t threadId);
    void findMarkerGraphReverseComplementVerticesThreadFunction2(size_t threadId);



    // Given a marker graph vertex, follow all of the contributing oriented
    // reads to their next vertex, but without moving forward more than
    // maxSkip markers.
    // In the returned vector, each entry correspond to a marker in the given vertex
    // (in the same order) and gives the next VertexId for that oriented read.
    // The next VertexId can be invalidVertexId if the oriented read has no vertices
    // past the starting VertexId.
    void findNextMarkerGraphVertices(
        MarkerGraphVertexId,
        uint32_t maxSkip,
        vector<MarkerGraphVertexId>&) const;



    // Clean up marker graph vertices that have duplicate markers
    // (more than one marker on the same oriented reads).
    // Such vertices are only generated when using --MarkerGraph.allowDuplicateMarkers.
public:
    void cleanupDuplicateMarkers(
        uint64_t threadCount,
        uint64_t minCoverage,
        uint64_t minCoveragePerStrand,
        double pattern1Threshold,
        bool pattern1CreateNewVertices,
        bool pattern2CreateNewVertices);
private:
    void cleanupDuplicateMarkersThreadFunction(size_t threadId);
    void cleanupDuplicateMarkersPattern1(
        MarkerGraph::VertexId,
        uint64_t minCoverage,
        uint64_t minCoveragePerStrand,
        bool createNewVertices,
        vector<MarkerDescriptor>&,
        vector<bool>& isDuplicateOrientedReadId,
        bool debug,
        ostream& out);
    void cleanupDuplicateMarkersPattern2(
        MarkerGraph::VertexId,
        uint64_t minCoverage,
        uint64_t minCoveragePerStrand,
        bool createNewVertices,
        vector<MarkerDescriptor>&,
        vector<bool>& isDuplicateOrientedReadId,
        bool debug,
        ostream& out);
    class CleanupDuplicateMarkersData {
    public:
        uint64_t minCoverage;
        uint64_t minCoveragePerStrand;
        double pattern1Threshold;
        bool pattern1CreateNewVertices;
        bool pattern2CreateNewVertices;
        uint64_t badVertexCount;    // Total number of vertices with duplicate markers.
        uint64_t pattern1Count;
        uint64_t pattern2Count;

        MarkerGraph::VertexId nextVertexId;

        // Get the next vertex id, then increment it in thread safe way.
        MarkerGraph::VertexId getAndIncrementNextVertexId()
        {
            return __sync_fetch_and_add(&nextVertexId, 1);
        }

    };
    CleanupDuplicateMarkersData cleanupDuplicateMarkersData;



    // Create marker graph edges.
public:
    void createMarkerGraphEdges(size_t threadCount);
    void accessMarkerGraphEdges(bool accessEdgesReadWrite, bool accessConnectivityReadWrite = false);
    void checkMarkerGraphEdgesIsOpen() const;
    void accessMarkerGraphConsensus();
private:
    void createMarkerGraphEdgesThreadFunction0(size_t threadId);
    void createMarkerGraphEdgesThreadFunction1(size_t threadId);
    void createMarkerGraphEdgesThreadFunction2(size_t threadId);
    void createMarkerGraphEdgesThreadFunction12(size_t threadId, size_t pass);
    void createMarkerGraphEdgesBySourceAndTarget(size_t threadCount);
    class CreateMarkerGraphEdgesData {
    public:
        vector< shared_ptr< MemoryMapped::Vector<MarkerGraph::Edge> > > threadEdges;
        vector< shared_ptr< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> > > threadEdgeMarkerIntervals;
    };
    CreateMarkerGraphEdgesData createMarkerGraphEdgesData;



    // "Strict" version of createMarkerGraphEdges.
    // Differences from createMarkerGraphEdges:
    // - Will only create edges in which all contributing oriented reads have
    //   exactly the same RLE sequence. If more than one distinct RLE sequence
    //   is present, the edge is split into two parallel edges.
    // - Enforces minEdgeCoverage and minEdgeCoveragePerStrand.
    //   An edge is not generated if the total number of oriented
    //   reads on the edge is less than minEdgeCoverage,
    //   of it the number of oriented reads on each strand is less
    //   than minEdgeCoveragePerStrand.
    // - The main loop is written differently - it loops over reads
    //   rather than marker graph vertices.
    // Because of these strict criteria, this version generates frequent breaks
    // in contiguity that must later be fixed by other means.
public:
    void createMarkerGraphEdgesStrict(
        uint64_t minEdgeCoverage,
        uint64_t minEdgeCoveragePerStrand,
        size_t threadCount);
private:
    void createMarkerGraphEdgesStrictPass1(size_t threadId);
    void createMarkerGraphEdgesStrictPass2(size_t threadId);
    void createMarkerGraphEdgesStrictPass12(size_t threadId, uint64_t pass);
    void createMarkerGraphEdgesStrictPass3(size_t threadId);
    class CreateMarkerGraphEdgesStrictData {
    public:
        uint64_t minEdgeCoverage;
        uint64_t minEdgeCoveragePerStrand;

        // Marker intervals, stored in a VectorOfVectors
        // indexed by the source vertex id, vertexId0.
        // Each interval stores the target vertex, vertexId1.
        class MarkerIntervalInfo {
        public:
            MarkerGraph::VertexId vertexId1;
            OrientedReadId orientedReadId;
            uint32_t ordinal0;
            uint32_t ordinal1;
            bool operator<(const MarkerIntervalInfo& that) const
            {
                return tie(vertexId1, orientedReadId, ordinal0, ordinal1) <
                    tie(that.vertexId1, that.orientedReadId, that.ordinal0, that.ordinal1);
            }
        };
        MemoryMapped::VectorOfVectors<MarkerIntervalInfo, uint64_t> markerIntervalInfos;

        // The edges and corresponding marker intervals found by each thread.
        vector< shared_ptr< MemoryMapped::Vector<MarkerGraph::Edge> > > threadEdges;
        vector< shared_ptr< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> > > threadEdgeMarkerIntervals;

        // Class used in createMarkerGraphEdgesStrictPass3.
        class MarkerIntervalInfo3 {
        public:
            MarkerInterval markerInterval;
            uint32_t overlap;       // The number of overlapping bases between the markers.
            vector<Base> sequence;  // RLE
            bool operator<(const MarkerIntervalInfo3& that) const
            {
                return tie(overlap, sequence, markerInterval) <
                    tie(that.overlap, that.sequence, that.markerInterval);
            }
        };
    };
    CreateMarkerGraphEdgesStrictData createMarkerGraphEdgesStrictData;



    // Write out the sets of parallel marker graph edges.
    // Only createMarkerGraphedgesStrict can create parallel edges.
public:
    void writeParallelMarkerGraphEdges() const;



    // Function createMarkerGraphSecondaryEdges can be called after createMarkerGraphEdgesStrict
    // to create a minimal amount of additional non-strict edges (secondary edges)
    // sufficient to restore contiguity.
    void createMarkerGraphSecondaryEdges(
        uint32_t secondaryEdgeMaxSkip,
        size_t threadCount);
private:
    void createMarkerGraphSecondaryEdges(
        uint32_t secondaryEdgeMaxSkip,
        bool aggressive,
        size_t threadCount);
public:



    // Cluster the oriented reads on a marker graph edge based on their sequence.
    // This returns a vector of connected components.
    // Each connected component is an index into the marker intervals for the edge.
    vector< vector<uint64_t> > clusterMarkerGraphEdgeOrientedReads(
        MarkerGraphEdgeId,
        double errorRateThreshold,
        bool debug) const;



    // Use clusterMarkerGraphEdgeOrientedReads to split secondary marker graph edges
    // where necessary.
    void splitMarkerGraphSecondaryEdges(
        double errorRateThreshold,
        uint64_t minCoverage,
        size_t threadCount);
    void splitMarkerGraphSecondaryEdgesThreadFunction(size_t threadId);
    class SplitMarkerGraphSecondaryEdgesData {
    public:
        double errorRateThreshold;
        uint64_t minCoverage;
        uint64_t initialSecondaryCount;
        uint64_t splitCount;
        uint64_t createdCount;

        // The new edges that were created by each thread.
        class Edge {
        public:
            MarkerGraphVertexId source;
            MarkerGraphVertexId target;
            vector<MarkerInterval> markerIntervals;
        };
        vector< vector<Edge> > threadEdges;
    };
    SplitMarkerGraphSecondaryEdgesData splitMarkerGraphSecondaryEdgesData;



    // Set marker graph edge flags to specified values for all marker graph edges.
    // Specify any value other than 0 or 1 leaves that flag unchanged.
    // Only useful for debugging.
    void setMarkerGraphEdgeFlags(
        uint8_t wasRemovedByTransitiveReduction,
        uint8_t wasPruned,
        uint8_t isSuperBubbleEdge,
        uint8_t isLowCoverageCrossEdge,
        uint8_t wasAssembled);



    // Find the reverse complement of each marker graph edge.
public:
    void findMarkerGraphReverseComplementEdges(size_t threadCount);
    void accessMarkerGraphReverseComplementEdge();
private:
    void findMarkerGraphReverseComplementEdgesThreadFunction1(size_t threadId);
    void findMarkerGraphReverseComplementEdgesThreadFunction2(size_t threadId);


    // Check that the marker graph is strand symmetric.
    // This can only be called after both findMarkerGraphReverseComplementVertices
    // and findMarkerGraphReverseComplementEdges have been called.
public:
    void checkMarkerGraphIsStrandSymmetric(size_t threadCount = 0);
private:
    void checkMarkerGraphIsStrandSymmetricThreadFunction1(size_t threadId);
    void checkMarkerGraphIsStrandSymmetricThreadFunction2(size_t threadId);



public:

    // Prune leaves from the strong subgraph of the global marker graph.
    void pruneMarkerGraphStrongSubgraph(size_t iterationCount);

private:



    // Private access functions for the global marker graph.
    // See the public section for some more that are callable from Python.

    // Find the vertex of the global marker graph that contains a given marker.
    // The marker is specified by the ReadId and Strand of the oriented read
    // it belongs to, plus the ordinal of the marker in the oriented read.
    // If the marker is not contained in any vertex, return
    // MarkerGraph::invalidVertexId.
    MarkerGraph::VertexId getGlobalMarkerGraphVertex(
        OrientedReadId,
        uint32_t ordinal) const;

    // Get pairs (ordinal, marker graph vertex id) for all markers of an oriented read.
    // The pairs are returned sorted by ordinal.
    void getMarkerGraphVertices(
        OrientedReadId,
        vector< pair<uint32_t, MarkerGraph::VertexId> >&);

    // Find the markers contained in a given vertex of the global marker graph.
    // The markers are stored as pairs(oriented read id, ordinal).
    void getGlobalMarkerGraphVertexMarkers(
        MarkerGraph::VertexId,
        vector< pair<OrientedReadId, uint32_t> >&) const;

    void getGlobalMarkerGraphVertexChildren(
        MarkerGraphVertexId,
        vector< pair<MarkerGraphVertexId, vector<MarkerInterval> > >&,
        vector< pair<MarkerGraphVertexId, MarkerInterval> >& workArea
        ) const;

    // Given two marker graph vertices, get the marker intervals
    // that a possible edge between the two vertices would have.
    void getMarkerIntervals(
        MarkerGraphVertexId,
        MarkerGraphVertexId,
        vector<MarkerInterval>&
        ) const;

    // Return true if a vertex of the global marker graph has more than
    // one marker for at least one oriented read id.
    bool isBadMarkerGraphVertex(MarkerGraph::VertexId) const;

    // Write csv files with detailed marker graph information.
    void debugWriteMarkerGraph(const string& fileNamePrefix = "") const;

    // Write a csv file with information on all marker graph vertices for which
    // isBadMarkerGraphVertex returns true.
public:
    void writeBadMarkerGraphVertices() const;
private:

    // Find out if a vertex is a forward or backward leaf of the pruned
    // strong subgraph of the marker graph.
    // A forward leaf is a vertex with out-degree 0.
    // A backward leaf is a vertex with in-degree 0.
    bool isForwardLeafOfMarkerGraphPrunedStrongSubgraph(MarkerGraph::VertexId) const;
    bool isBackwardLeafOfMarkerGraphPrunedStrongSubgraph(MarkerGraph::VertexId) const;

    // Given an edge of the pruned strong subgraph of the marker graph,
    // return the next/previous edge in the linear chain the edge belongs to.
    // If the edge is the last/first edge in its linear chain, return MarkerGraph::invalidEdgeId.
    MarkerGraphEdgeId nextEdgeInMarkerGraphPrunedStrongSubgraphChain(MarkerGraphEdgeId) const;
    MarkerGraphEdgeId previousEdgeInMarkerGraphPrunedStrongSubgraphChain(MarkerGraphEdgeId) const;

    // Return the out-degree or in-degree (number of outgoing/incoming edges)
    // of a vertex of the pruned strong subgraph of the marker graph.
    size_t markerGraphPrunedStrongSubgraphOutDegree(MarkerGraph::VertexId) const;
    size_t markerGraphPrunedStrongSubgraphInDegree (MarkerGraph::VertexId) const;

    // Return true if an edge disconnects the local subgraph.
    bool markerGraphEdgeDisconnectsLocalStrongSubgraph(
        MarkerGraphEdgeId edgeId,
        size_t maxDistance,

        // Work areas, to reduce memory allocation activity.

        // Each of these two must be sized maxDistance+1.
        array<vector< vector<MarkerGraphEdgeId> >, 2>& verticesByDistance,

        // Each of these two must be sized globalMarkerGraphVertices.size()
        // and set to all false on entry.
        // It is left set to all false on exit, so it can be reused.
        array<vector<bool>, 2>& vertexFlags
        ) const;



    // Each oriented read corresponds to a path in the marker graph.
    // This function computes a subset of that path
    // covering the specified range of marker ordinals for the given
    // oriented read.
    void computeOrientedReadMarkerGraphPath(
        OrientedReadId,
        uint32_t firstOrdinal,
        uint32_t lastOrdinal,
        vector<MarkerGraphEdgeId>& path,
        vector< pair<uint32_t, uint32_t> >& pathOrdinals
        ) const;

    // Create the marker connectivity graph starting with a given marker.
    void createMarkerConnectivityGraph(
        OrientedReadId,
        uint32_t ordinal,
        bool useReadGraphAlignmentsOnly,
        MarkerConnectivityGraph&) const;
    void createMarkerConnectivityGraph(
        OrientedReadId,
        uint32_t ordinal,
        bool useReadGraphAlignmentsOnly,
        MarkerConnectivityGraph&,
        MarkerConnectivityGraphVertexMap&) const;

    // Compute an alignment between two oriented reads
    // induced by the marker graph. See InducedAlignment.hpp for more
    // information.
    void computeInducedAlignment(
        OrientedReadId,
        OrientedReadId,
        InducedAlignment&
    );

    // Compute induced alignments between an oriented read orientedReadId0
    // and the oriented reads stored sorted in orientedReadIds1.
    void computeInducedAlignments(
        OrientedReadId orientedReadId0,
        const vector<OrientedReadId>& orientedReadIds1,
        vector<InducedAlignment>& inducedAlignments);

    // Fill in compressed ordinals of an InducedAlignment.
    void fillCompressedOrdinals(
        OrientedReadId,
        OrientedReadId,
        InducedAlignment&);

    // Find the markers aligned to a given marker.
    // This is slow and cannot be used during assembly.
    void findAlignedMarkers(
        OrientedReadId, uint32_t ordinal,
        bool useReadGraphAlignmentsOnly,
        vector< pair<OrientedReadId, uint32_t> >&) const;



    // Extract a local subgraph of the global marker graph.
    bool extractLocalMarkerGraph(
        OrientedReadId,
        uint32_t ordinal,
        uint64_t distance,
        int timeout,                 // Or 0 for no timeout.
        uint64_t minVertexCoverage,
        uint64_t minEdgeCoverage,
        bool useWeakEdges,
        bool usePrunedEdges,
        bool useSuperBubbleEdges,
        bool useLowCoverageCrossEdges,
        bool useRemovedSecondaryEdges,
        LocalMarkerGraph&
        );
    bool extractLocalMarkerGraph(
        MarkerGraph::VertexId,
        uint64_t distance,
        int timeout,                 // Or 0 for no timeout.
        uint64_t minVertexCoverage,
        uint64_t minEdgeCoverage,
        bool useWeakEdges,
        bool usePrunedEdges,
        bool useSuperBubbleEdges,
        bool useLowCoverageCrossEdges,
        bool useRemovedSecondaryEdges,
        LocalMarkerGraph&
        );

    // Compute consensus sequence for a vertex of the marker graph.
    void computeMarkerGraphVertexConsensusSequence(
        MarkerGraph::VertexId,
        vector<Base>& sequence,
        vector<uint32_t>& repeatCounts
        );



    // Class used to store spoa details needed by the http server.
    // See computeMarkerGraphEdgeConsensusSequenceUsingSpoa for details.
    class ComputeMarkerGraphEdgeConsensusSequenceUsingSpoaDetail {
    public:

        // If there is a very long marker interval,
        // the shortest sequence is used as consensus.
        // In that case, this flag is set and nothing else is stored.
        bool hasLongMarkerInterval;

        // Assembly mode: 1=overlapping bases, 2=intervening bases.
        int assemblyMode;   // 1 or 2.

        // Data stored when hasLongMarkerInterval is set.
        size_t iShortest;

        // Data stored for assembly mode 1.



        // Data stored for assembly mode 2.

        // The alignment for each distinct sequence.
        // Indexed by distinct sequence index.
        vector<string> msa;

        // The consensus, including gap bases.
        vector<AlignedBase> alignedConsensus;
        vector<uint8_t> alignedRepeatCounts;

        // The indexes of oriented reads that have each of the distinct sequences.
        // Indexed by distinct sequence index (same as the index used
        // for the msa vector above).
        vector< vector<size_t> > distinctSequenceOccurrences;

        // The alignment row corresponding to each oriented read.
        vector<int> alignmentRow;
    };



    // Get the RLE sequence implied by a MarkerInterval.
    // If the markers overlap, returns the number of
    // overlapping RLE bases in overlappingRleBaseCount
    // and empty rleSequence.
    // Otherwise, returns zero overlappingRleBaseCount
    // and the intervening sequence in rleSequence
    // (which can be empty if the two markers are exactly adjacent).
    void getMarkerIntervalRleSequence(
        const MarkerInterval&,
        uint64_t& overlappingRleBaseCount,
        vector<Base>& rleSequence) const;



    // Use spoa to compute consensus sequence for an edge of the marker graph.
    // This does not include the bases corresponding to the flanking markers.
    void computeMarkerGraphEdgeConsensusSequenceUsingSpoa(
        MarkerGraphEdgeId,
        uint32_t markerGraphEdgeLengthThresholdForConsensus,
        const std::unique_ptr<spoa::AlignmentEngine>& spoaAlignmentEngine,
        spoa::Graph& spoaAlignmentGraph,
        vector<Base>& sequence,
        vector<uint32_t>& repeatCounts,
        uint8_t& overlappingBaseCount,
        ComputeMarkerGraphEdgeConsensusSequenceUsingSpoaDetail&,
        vector< pair<uint32_t, CompressedCoverageData> >* coverageData // Optional
        );



    // Simplify the marker graph.
    // The first argument is a number of marker graph edges.
    // See the code for detail on its meaning and how it is used.
public:
    void simplifyMarkerGraph(
        const vector<size_t>& maxLength, // One value for each iteration.
        bool debug);
private:
    void simplifyMarkerGraphIterationPart1(
        size_t iteration,
        size_t maxLength,
        bool debug);
    void simplifyMarkerGraphIterationPart2(
        size_t iteration,
        size_t maxLength,
        bool debug);



    // Create a coverage histogram for vertices and edges of the
    // marker graph. This counts all vertices that are not isolated
    // (are connected to no edges that are not marked removed)
    // and all edges that are not marked as removed.
    // Output is to csv files.
public:
    void computeMarkerGraphCoverageHistogram();


    // In the assembly graph, each vertex corresponds to a linear chain
    // of edges in the pruned strong subgraph of the marker graph.
    // A directed vertex A->B is created if the last marker graph vertex
    // of the edge chain corresponding to A coincides with the
    // first marker graph vertex of the edge chain corresponding to B.
    shared_ptr<AssemblyGraph> assemblyGraphPointer;
    void removeAssemblyGraph()
    {
        assemblyGraphPointer.reset();
    }
    void createAssemblyGraphVertices();
    void accessAssemblyGraphVertices();
    void createAssemblyGraphEdges();
    void accessAssemblyGraphEdgeLists();
    void accessAssemblyGraphEdges();
    void accessAssemblyGraphOrientedReadsByEdge();
    void writeAssemblyGraph(const string& fileName) const;
    void pruneAssemblyGraph(uint64_t pruneLength);

    // Gather and write out all reads that contributed to
    // each assembly graph edge.
    void gatherOrientedReadsByAssemblyGraphEdge(size_t threadCount);
    void writeOrientedReadsByAssemblyGraphEdge();
private:
    void gatherOrientedReadsByAssemblyGraphEdgePass1(size_t threadId);
    void gatherOrientedReadsByAssemblyGraphEdgePass2(size_t threadId);
    void gatherOrientedReadsByAssemblyGraphEdgePass(int);

    // Extract a local assembly graph from the global assembly graph.
    // This returns false if the timeout was exceeded.
    bool extractLocalAssemblyGraph(
        AssemblyGraphEdgeId,
        int distance,
        double timeout,
        LocalAssemblyGraph&) const;
public:
    void colorGfaBySimilarityToSegment(
        AssemblyGraphEdgeId,
        uint64_t minVertexCount,
        uint64_t minEdgeCount);


    // Compute consensus repeat counts for each vertex of the marker graph.
    void assembleMarkerGraphVertices(size_t threadCount);
    void accessMarkerGraphVertexRepeatCounts();
private:
    void assembleMarkerGraphVerticesThreadFunction(size_t threadId);
public:



    // Optional computation of coverage data for marker graph vertices.
    // This is only called if Assembly.storeCoverageData in shasta.conf is True.
    void computeMarkerGraphVerticesCoverageData(size_t threadCount);
private:
    void computeMarkerGraphVerticesCoverageDataThreadFunction(size_t threadId);
    class ComputeMarkerGraphVerticesCoverageDataData {
    public:

        // The results computed by each thread.
        // For each threadId:
        // - threadVertexIds[threadId] contains the vertex ids processed by each thread.
        // - threadVertexCoverageData[threadId] contains the coverage data for those vertices.
        vector< shared_ptr<
            MemoryMapped::Vector<MarkerGraph::VertexId> > > threadVertexIds;
        vector< shared_ptr<
            MemoryMapped::VectorOfVectors<pair<uint32_t, CompressedCoverageData>, uint64_t> > >
            threadVertexCoverageData;
    };
    ComputeMarkerGraphVerticesCoverageDataData computeMarkerGraphVerticesCoverageDataData;



    // Find the set of assembly graph edges encountered on a set
    // of edges in the marker graph. The given marker graph edges
    // could form a path, but don't have to.
    void findAssemblyGraphEdges(
        const vector<MarkerGraphEdgeId>& markerGraphEdges,
        vector<AssemblyGraphEdgeId>& assemblyGraphEdges
        ) const;



    // Pseudo-paths.
    // An oriented read corresponds to a path (sequence of adjacent edges)
    // in the marker graph, which
    // can be computed via computeOrientedReadMarkerGraphPath.
    // That path encounters a sequence of assembly graph edges,AssemblyGraph::EdgeId
    // which is not necessarily a path in the assembly graph
    // because not all marker graph edges belong to an assembly graph edge.
    // We call this sequence the pseudo-path of an oriented read in the assembly graph.
public:
    class PseudoPathEntry {
    public:
        AssemblyGraphEdgeId segmentId;

        // The first and last ordinal on the oriented read
        // where this assembly graph edge (segment) is encountered.
        uint32_t firstOrdinal;
        uint32_t lastOrdinal;

        // The first and last position in the assembly graph edge
        // (segment) where this oriented read is encountered.
        uint32_t firstPosition;
        uint32_t lastPosition;

        // The number of marker graph edges on the oriented read
        // where this assembly graph edge (segment) is encountered.
        uint32_t markerGraphEdgeCount;
    };
    using PseudoPath = vector<PseudoPathEntry>;
    void computePseudoPath(
        OrientedReadId,

        // The marker graph path computed using computeOrientedReadMarkerGraphPath.
        // This is computed by this function - it does not neet to be filled in
        // in advance.
        vector<MarkerGraphEdgeId>& path,
        vector< pair<uint32_t, uint32_t> >& pathOrdinals,

        // The pseudo-path computed by this function.
        PseudoPath&) const;
    void writePseudoPath(ReadId, Strand) const;
    static void getPseudoPathSegments(const PseudoPath&, vector<AssemblyGraphEdgeId>&);



    // Detangle the AssemblyGraph.
    void detangle();    // detangleMethod 1
    void detangle2(     // detangleMethod 2
        uint64_t diagonalReadCountMin,
        uint64_t offDiagonalReadCountMax,
        double detangleOffDiagonalRatio
         );



    // CompressedAssemblyGraph.
    // Note that we have no persistent version of this.
    // It must be created from scratch each time.
    void createCompressedAssemblyGraph();
    shared_ptr<CompressedAssemblyGraph> compressedAssemblyGraph;
    void colorCompressedAssemblyGraph(const string&);


public:
    // Mark as isLowCoverageCrossEdge all low coverage cross edges
    // of the assembly graph and the corresponding marker graph edges.
    // These edges are then considered removed.
    // An edge v0->v1 of the assembly graph is a cross edge if:
    // - in-degree(v0)=1, out-degree(v0)>1
    // - in-degree(v1)>1, out-degree(v1)=1
    // A cross edge is marked as isCrossEdge if its average edge coverage
    // is <= crossEdgeCoverageThreshold.
    void removeLowCoverageCrossEdges(uint32_t crossEdgeCoverageThreshold);



    // Assemble consensus sequence and repeat counts for each marker graph edge.
    void assembleMarkerGraphEdges(
        size_t threadCount,

        // This controls when we give up trying to compute consensus for long edges.
        uint32_t markerGraphEdgeLengthThresholdForConsensus,

        // Request storing detailed coverage information in binary format.
        bool storeCoverageData,

        // Request assembling all edges (used by Mode 2 assembly)
        bool assembleAllEdges
        );
private:
    void assembleMarkerGraphEdgesThreadFunction(size_t threadId);
    class AssembleMarkerGraphEdgesData {
    public:

        // The arguments to assembleMarkerGraphEdges, stored here so
        // they are accessible to the threads.
        uint32_t markerGraphEdgeLengthThresholdForConsensus;
        bool storeCoverageData;
        bool assembleAllEdges;

        // The results computed by each thread.
        // For each threadId:
        // threadEdgeIds[threadId] contains the edge ids processed by each thread.
        // threadEdgeConsensusSequence[threadId]  and
        // threadEdgeConsensusOverlappingBaseCount[threadId] contains the corresponding
        // consensus sequence and repeat counts.
        // These are temporary data which are eventually gathered into
        // MarkerGraph::edgeConsensus and MarkerGraph::edgeConsensusOverlappingBaseCount
        // before assembleMarkerGraphEdges completes.
        // See their definition for more details about their meaning.
        vector< shared_ptr< MemoryMapped::Vector<MarkerGraphEdgeId> > > threadEdgeIds;
        vector< shared_ptr< MemoryMapped::VectorOfVectors<pair<Base, uint8_t>, uint64_t> > > threadEdgeConsensus;
        vector< shared_ptr< MemoryMapped::Vector<uint8_t> > > threadEdgeConsensusOverlappingBaseCount;

        vector< shared_ptr<
            MemoryMapped::VectorOfVectors<pair<uint32_t, CompressedCoverageData>, uint64_t> > >
            threadEdgeCoverageData;
    };
    AssembleMarkerGraphEdgesData assembleMarkerGraphEdgesData;

    // Access coverage data for vertices and edges of the marker graph.
    // This is only available if the run had Assembly.storeCoverageData set to True
    // in shasta.conf.
public:
    void accessMarkerGraphCoverageData();
private:


    // Assemble sequence for an edge of the assembly graph.
    void assembleAssemblyGraphEdge(
        AssemblyGraphEdgeId,
        bool storeCoverageData,
        AssembledSegment&);
    // Lower level version that works on a generic marker graph path.
    void assembleAssemblyGraphEdge(
        const span<const MarkerGraphEdgeId>&,
        bool storeCoverageData,
        AssembledSegment&);
public:
    AssembledSegment assembleAssemblyGraphEdge(
        AssemblyGraphEdgeId,
        bool storeCoverageData);


    // Assemble sequence for all edges of the assembly graph.
    void assemble(
        size_t threadCount,
        uint32_t storeCoverageDataCsvLengthThreshold);
    void accessAssemblyGraphSequences();
    void computeAssemblyStatistics();
private:
    class AssembleData {
    public:
        uint32_t storeCoverageDataCsvLengthThreshold;

        // The results created by each thread.
        // All indexed by threadId.
        vector< vector<AssemblyGraphEdgeId> > edges;
        vector< shared_ptr<LongBaseSequences> > sequences;
        vector< shared_ptr<MemoryMapped::VectorOfVectors<uint8_t, uint64_t> > > repeatCounts;
        void allocate(size_t threadCount);
        void free();
    };
    AssembleData assembleData;
    void assembleThreadFunction(size_t threadId);



    // Assemble the RLE sequence of a path of the marker graph, under the assumption
    // that, for each edge, all oriented reads have exactly the same sequence.
    // This will be the case if edges were created by Assembler::createMarkerGraphEdgesStrict.
public:
    void assembleMarkerGraphPathRleStrict(
        span<const MarkerGraphEdgeId> path,
        vector<Base>& rleSequence
    ) const;
    // Same, but for an assembly graph edge.
    void assembleAssemblyGraphEdgeRleStrict(
        AssemblyGraphEdgeId,
        vector<Base>& rleSequence
    ) const;



    // Write the assembly graph in GFA 1.0 format defined here:
    // https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
    void writeGfa1(const string& fileName);
    void writeGfa1BothStrands(const string& fileName);
    void writeGfa1BothStrandsNoSequence(const string& fileName);
private:
    // Construct the CIGAR string given two vectors of repeat counts.
    // Used by writeGfa1.
    static void constructCigarString(
        const span<uint8_t>& repeatCounts0,
        const span<uint8_t>& repeatCounts1,
        string&
        );

public:

    // Write assembled sequences in FASTA format.
    void writeFasta(const string& fileName);



    // Write a csv file that can be used to color the double-stranded GFA
    // in Bandage based on the presence of two oriented reads
    // on each assembly graph edge.
    // Red    =  only oriented read id 0 is present
    // Blue   =  only oriented read id 1 is present
    // Purple =  both oriented read id 0 and oriented read id 1 are present
    // Grey   =  neither oriented read id 0 nor oriented read id 1 are present
    void colorGfaWithTwoReads(
        ReadId readId0, Strand strand0,
        ReadId readId1, Strand strand1,
        const string& fileName
        ) const;



    // Color key segments in the gfa file.
    // A segment (assembly graph edge) v0->v1 is a key segment if in-degree(v0)<2 and
    // out_degree(v)<2, that is, there is no uncertainty on what preceeds
    // and follows the segment.
    void colorGfaKeySegments(const string& fileName) const;



    // Write a csv file describing the marker graph path corresponding to an
    // oriented read and the corresponding pseudo-path on the assembly graph.
    void writeOrientedReadPath(ReadId, Strand, const string& fileName) const;



    // Analyze pseudo-paths of oriented reads.
    void alignPseudoPaths(ReadId, Strand, ReadId, Strand);



    // Data and functions used for the http server.
    // This function puts the server into an endless loop
    // of processing requests.
    void writeHtmlBegin(ostream&) const;
    void writeHtmlEnd(ostream&) const;
    void writeAssemblySummary(ostream&);
    void writeAssemblySummaryBody(ostream&);
    void writeAssemblySummaryJson(ostream&);
    void writeAssemblyIndex(ostream&) const;
    static void writeStyle(ostream& html);


    void writeNavigation(ostream&) const;
    void writeNavigation(
        ostream& html,
        const string& title,
        const vector<pair <string, string> >&) const;

    static void writePngToHtml(
        ostream& html,
        const string& pngFileName,
        const string useMap = ""
        );
    static void writeGnuPlotPngToHtml(
        ostream& html,
        int width,
        int height,
        const string& gnuplotCommands);

    void fillServerFunctionTable();
    void processRequest(
        const vector<string>& request,
        ostream&,
        const BrowserInformation&) override;
    void exploreSummary(const vector<string>&, ostream&);
    void exploreRead(const vector<string>&, ostream&);
    void exploreReadRaw(const vector<string>&, ostream&);
    void exploreReadRle(const vector<string>&, ostream&);
    void blastRead(const vector<string>&, ostream&);
    void exploreAlignmentCandidateGraph(const vector<string>& request, ostream& html);
    void exploreAlignments(const vector<string>&, ostream&);
    void exploreAlignmentCoverage(const vector<string>&, ostream&);
    void exploreAlignment(const vector<string>&, ostream&);
    void alignSequencesInBaseRepresentation(const vector<string>&, ostream&);
    void alignSequencesInMarkerRepresentation(const vector<string>&, ostream&);
    void exploreAlignmentGraph(const vector<string>&, ostream&);
    void exploreReadGraph(const vector<string>&, ostream&);
    void exploreUndirectedReadGraph(const vector<string>&, ostream&);
    void exploreDirectedReadGraph(const vector<string>&, ostream&);
    void exploreCompressedAssemblyGraph(const vector<string>&, ostream&);
    static bool parseCommaSeparatedReadIDs(string& commaSeparatedReadIds, vector<OrientedReadId>& readIds, ostream& html);
    static void addScaleSvgButtons(ostream&, uint64_t sizePixels);
    class HttpServerData {
    public:
        shared_ptr<LocalAlignmentCandidateGraph> referenceOverlapGraph;

        using ServerFunction = void (Assembler::*) (
            const vector<string>& request,
            ostream&);
        std::map<string, ServerFunction> functionTable;
        string docsDirectory;
        string referenceFastaFileName = "reference.fa";

        const AssemblerOptions* assemblerOptions = 0;

        void createGraphEdgesFromOverlapMap(const ReferenceOverlapMap& overlapMap);

    };
    HttpServerData httpServerData;

    // For the display of the alignment candidate graph, we can optionally
    // specify a PAF file containing alignments of reads to the reference.
    // Persistent data structures from loading the PAF are stored as
    // members of HttpServerData
    void loadAlignmentsPafFile(const string& alignmentsPafFileAbsolutePath);

    // Display alignments in an html table.
    void displayAlignments(
        OrientedReadId,
        const vector< pair<OrientedReadId, AlignmentInfo> >&,
        bool showIsInReadGraphFlag,
        ostream&) const;
    void displayAlignment(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        const AlignmentInfo& alignment,
        ostream&) const;


    // Functions and data used by the http server
    // for display of the local marker graph.
    void exploreMarkerGraph(const vector<string>&, ostream&);
    void getLocalMarkerGraphRequestParameters(
        const vector<string>&,
        LocalMarkerGraphRequestParameters&) const;
    void exploreMarkerGraphVertex(const vector<string>&, ostream&);
    void exploreMarkerGraphEdge(const vector<string>&, ostream&);
    void exploreMarkerCoverage(const vector<string>&, ostream&);
    void exploreMarkerGraphInducedAlignment(const vector<string>&, ostream&);
    void followReadInMarkerGraph(const vector<string>&, ostream&);
    void exploreMarkerConnectivity(const vector<string>&, ostream&);
    void renderEditableAlignmentConfig(
        const int method,
        const uint64_t maxSkip,
        const uint64_t maxDrift,
        const uint32_t maxMarkerFrequency,
        const uint64_t minAlignedMarkerCount,
        const double minAlignedFraction,
        const uint64_t maxTrim,
        const int matchScore,
        const int mismatchScore,
        const int gapScore,
        const double downsamplingFactor,
        int bandExtend,
        int maxBand,
        uint64_t align4DeltaX,
        uint64_t align4DeltaY,
        uint64_t align4MinEntryCountPerCell,
        uint64_t align4MaxDistanceFromBoundary,
        ostream& html
    );
    void writeColorPicker(ostream& html, string svgId);

    void writeMakeAllTablesCopyable(ostream&) const;

    // Do bulk sampling of reads and accumulate stats about their alignments
    void assessAlignments(const vector<string>& request, ostream& html);
    void sampleReads(vector<OrientedReadId>& sample, uint64_t n);
    void sampleReads(vector<OrientedReadId>& sample, uint64_t n, uint64_t minLength, uint64_t maxLength);
    void sampleReadsFromDeadEnds(
            vector<OrientedReadId>& sample,
            vector<bool>& isLeftEnd,
            uint64_t n);

    void sampleReadsFromDeadEnds(
            vector<OrientedReadId>& sample,
            vector<bool>& isLeftEnd,
            uint64_t n,
            uint64_t minLength,
            uint64_t maxLength);

    void countDeadEndOverhangs(
            const vector<pair<OrientedReadId, AlignmentInfo> >& allAlignmentInfo,
            const vector<bool>& isLeftEnd,
            Histogram2& overhangLengths,
            uint32_t minOverhang);

    // Compute all alignments for a given read.
    // This can be slow for large assemblies,
    // and therefore the computation in multithreaded.
    void computeAllAlignments(const vector<string>&, ostream&);
    void computeAllAlignmentsThreadFunction(size_t threadId);
    class ComputeAllAlignmentsData {
    public:
        OrientedReadId orientedReadId0;
        size_t minMarkerCount;
        size_t maxSkip;
        size_t maxDrift;
        uint32_t maxMarkerFrequency;
        size_t minAlignedMarkerCount;
        double minAlignedFraction;
        size_t maxTrim;
        int method;
        int matchScore;
        int mismatchScore;
        int gapScore;
        double downsamplingFactor;
        int bandExtend;
        int maxBand;
        uint64_t align4DeltaX;
        uint64_t align4DeltaY;
        uint64_t align4MinEntryCountPerCell;
        uint64_t align4MaxDistanceFromBoundary;
        // The alignments found by each thread.
        vector< vector< pair<OrientedReadId, AlignmentInfo> > > threadAlignments;
    };
    ComputeAllAlignmentsData computeAllAlignmentsData;


    // Access all available assembly data, without thorwing an exception
    // on failures.
public:
    void accessAllSoft();

    // Store assembly time.
    void storeAssemblyTime(
        double elapsedTimeSeconds,
        double averageCpuUtilization);

    void storePeakMemoryUsage(uint64_t peakMemoryUsage);

    // Functions and data used by the http server
    // for display of the local assembly graph.
private:
    void exploreAssemblyGraph(const vector<string>&, ostream&);
    class LocalAssemblyGraphRequestParameters {
    public:
        AssemblyGraphEdgeId edgeId;
        bool edgeIdIsPresent;
        uint32_t maxDistance;
        bool maxDistanceIsPresent;
        bool useDotLayout;
        bool showVertexLabels;
        bool showEdgeLabels;
        uint32_t sizePixels;
        bool sizePixelsIsPresent;
        double timeout;
        bool timeoutIsPresent;
        void writeForm(ostream&, AssemblyGraphEdgeId edgeCount) const;
        bool hasMissingRequiredParameters() const;
    };
    void getLocalAssemblyGraphRequestParameters(
        const vector<string>&,
        LocalAssemblyGraphRequestParameters&) const;
    void exploreAssemblyGraphEdge(const vector<string>&, ostream&);
    void exploreAssemblyGraphEdgesSupport(const vector<string>&, ostream&);


    // Http server functionality specific to mode 3 assembly.
    void exploreMode3AssemblyGraph(const vector<string>&, ostream&);
    void exploreMode3AssemblyGraphSegment(const vector<string>&, ostream&);
    void exploreMode3AssemblyGraphSegmentPair(const vector<string>&, ostream&);
    void exploreMode3AssemblyGraphLink(const vector<string>&, ostream&);
    void exploreMode3MetaAlignment(const vector<string>&, ostream&);
    void exploreMode3AssemblyPath(const vector<string>&, ostream&);
    void exploreMode3LinkAssembly(const vector<string>&, ostream&);



    // Set up the ConsensusCaller used to compute the "best"
    // base and repeat count at each assembly position.
    // The argument to setupConsensusCaller specifies
    // the consensus caller to be used.
    // It can be one of the following:
    // - Modal
    //   Selects the SimpleConsensusCaller.
    // - Median
    //   Selects the MedianConsensusCaller.
    // - Bayesian:fileName
    //   Selects the SimpleBayesianConsensusCaller,
    //   using fileName as the configuration file.
    //   Filename must be an absolute path (it must begin with "/").
public:
    void setupConsensusCaller(const string&);
private:
    shared_ptr<ConsensusCaller> consensusCaller;
public:



    // Assembly graph for mode 2 assembly.
    shared_ptr<AssemblyGraph2> assemblyGraph2Pointer;
    void createAssemblyGraph2(
        uint64_t pruneLength,
        const Mode2AssemblyOptions&,
        size_t threadCount,
        bool debug);


    // Mode 3 assembly.
    void mode3Assembly(
        size_t threadCount);
    shared_ptr<mode3::AssemblyGraph> assemblyGraph3Pointer;
    void accessMode3AssemblyGraph();
    void analyzeMode3Subgraph(const vector<uint64_t>& segmentIds);
    void createMode3PathGraph();



public:
    void test();
};

#endif
