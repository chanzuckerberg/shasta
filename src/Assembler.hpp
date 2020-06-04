#ifndef SHASTA_ASSEMBLER_HPP
#define SHASTA_ASSEMBLER_HPP

// Shasta.
#include "Alignment.hpp"
#include "AlignmentCandidates.hpp"
#include "AssembledSegment.hpp"
#include "AssemblyGraph.hpp"
#include "ConflictReadGraph.hpp"
#include "Coverage.hpp"
#include "DirectedReadGraph.hpp"
#include "dset64-gccAtomic.hpp"
#include "HttpServer.hpp"
#include "InducedAlignment.hpp"
#include "Kmer.hpp"
#include "LongBaseSequence.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "MemoryMappedObject.hpp"
#include "MultithreadedObject.hpp"
#include "OrientedReadPair.hpp"
#include "PhasingData.hpp"
#include "ReadGraph.hpp"
#include "ReadFlags.hpp"
#include "ReadId.hpp"

// Standard library.
#include "memory.hpp"
#include "string.hpp"
#include "tuple.hpp"


namespace shasta {

    class Assembler;
    class AssemblerInfo;
    class Alignment;
    class AlignmentGraph;
    class AlignmentInfo;
    class AssemblerOptions;
    class AssembledSegment;
    class CompressedAssemblyGraph;
    class ConsensusCaller;
    class LocalAssemblyGraph;
    class LocalAlignmentGraph;
    class LocalReadGraph;
    class SegmentGraph;

#ifdef SHASTA_HTTP_SERVER
    class LocalMarkerGraph;
#endif

    // Write an html form to select strand.
    void writeStrandSelection(
        ostream&,               // The html stream to write the form to.
        const string& name,     // The selection name.
        bool select0,           // Whether strand 0 is selected.
        bool select1);          // Whether strand 1 is selected.
}


// Sanity check that we are compiling on x86_64.
#if !__x86_64__
#error "Shasta can only be built on an x86_64 machine (64-bit Intel/AMD)"
#endif



// Class used to store various pieces of assembler information in shared memory.
class shasta::AssemblerInfo {
public:

    // The length of k-mers used to define markers.
    size_t k;

    // The page size in use for this run.
    size_t largeDataPageSize;



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

    // Other read statistics.
    size_t palindromicReadCount = 0;
    size_t chimericReadCount = 0;
    uint64_t isolatedReadCount = 0;
    uint64_t isolatedReadBaseCount = 0;

    // Marker graph statistics.
    size_t markerGraphVerticesNotIsolatedCount = 0;
    size_t markerGraphEdgesNotRemovedCount = 0;

    // Assembly graph statistics.
    size_t assemblyGraphAssembledEdgeCount = 0;
    size_t totalAssembledSegmentLength = 0;
    size_t longestAssembledSegmentLength = 0;
    size_t assembledSegmentN50 = 0;

    // Performance statistics.
    double assemblyElapsedTimeSeconds = 0.;
    double averageCpuUtilization;

};



class shasta::Assembler :
    public MultithreadedObject<Assembler>
#ifdef SHASTA_HTTP_SERVER
    , public HttpServer
#endif
    {
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
        size_t largeDataPageSize);

    // Add reads.
    // The reads in the specified file are added to those already previously present.
    void addReads(
        const string& fileName,
        size_t minReadLength,
        bool noCache,
        size_t threadCount);

    // Create a histogram of read lengths.
    void histogramReadLength(const string& fileName);

    // Function to write one or all reads in Fasta format.
    void writeReads(const string& fileName);
    void writeRead(ReadId, const string& fileName);
    void writeOrientedRead(ReadId, Strand, const string& fileName);



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

        // Alignment method.
        int alignmentMethod,

        // Marker frequency threshold.
        // When computing an alignment between two oriented reads,
        // marker kmers that appear more than this number of times
        // in either of the two oriented reads are discarded
        // (in both oriented reads).
        // Change to size_t when conversion completed.
        uint32_t maxMarkerFrequency,

        // The maximum ordinal skip to be tolerated between successive markers
        // in the alignment.
        size_t maxSkip,

        // The maximum relative ordinal drift to be tolerated between successive markers
        // in the alignment.
        size_t maxDrift,

        // Minimum number of alignment markers for an alignment to be used.
        size_t minAlignedMarkerCount,

        // The minimum fraction of aligned markers for an alignment to be used.
        double minAlignedFraction,

        // Maximum left/right trim (in bases) for an alignment to be used.
        size_t maxTrim,

        // Scores used to compute method 1 alignments.
        int matchScore,
        int mismatchScore,
        int gapScore,

        // Parameters for alignment method 3.
        double downsamplingFactor,
        int bandExtend,

        // If true, discard containment alignments.
        bool suppressContainments,

        // If true, store good alignments in a compressed format.
        bool storeAlignments,

        // Number of threads. If zero, a number of threads equal to
        // the number of virtual processors is used.
        size_t threadCount
    );
    void accessAlignmentData();


    // Loop over all alignments in the read graph
    // to create vertices of the global marker graph.
    // Throw away vertices with coverage (number of markers)
    // less than minCoverage or more than maxCoverage.
    // Also throw away "bad" vertices - that is, vertices
    // with more than one marker on the same oriented read.
    void createMarkerGraphVertices(

        // The method to be used to compute alignments.
        int alignMethod,

        // The maximum frequency of marker k-mers to be used in
        // computing alignments.
        uint32_t maxMarkerFrequency,

        // The maximum ordinal skip to be tolerated between successive markers
        // in the alignment.
        size_t maxSkip,

        // The maximum ordinal drift to be tolerated between successive markers
        // in the alignment.
        size_t maxDrift,

        // Scores for method 1 alignments.
        int matchScore,
        int mismatchScore,
        int gapScore,

        // Parameters for method 3 alignments.
        double downsamplingFactor,
        int bandExtend,

        // The method used to create the read graph.
        // This affects which alignments are used to create the marker graph.
        int readGraphCreationMethod,

        // Minimum coverage (number of markers) for a vertex
        // of the marker graph to be kept.
        size_t minCoverage,

        // Maximum coverage (number of markers) for a vertex
        // of the marker graph to be kept.
        size_t maxCoverage,

        // Number of threads. If zero, a number of threads equal to
        // the number of virtual processors is used.
        size_t threadCount
    );



    // Python-callable access functions for the global marker graph.
    // See the private section for some more not callable from Python.
    void accessMarkerGraphVertices(bool readWriteAccess = false);

    // Find the vertex of the global marker graph that contains a given marker.
    // The marker is specified by the ReadId and Strand of the oriented read
    // it belongs to, plus the ordinal of the marker in the oriented read.
    MarkerGraph::VertexId getGlobalMarkerGraphVertex(
        ReadId,
        Strand,
        uint32_t ordinal) const;

    // Find the markers contained in a given vertex of the global marker graph.
    // Returns the markers as tuples(read id, strand, ordinal).
    vector< tuple<ReadId, Strand, uint32_t> >
        getGlobalMarkerGraphVertexMarkers(MarkerGraph::VertexId) const;

    // Find the children or parents of a vertex of the global marker graph.
    vector<MarkerGraph::VertexId>
        getGlobalMarkerGraphVertexChildren(
        MarkerGraph::VertexId) const;
    vector<MarkerGraph::VertexId>
        getGlobalMarkerGraphVertexParents(
        MarkerGraph::VertexId) const;



    // Python-callable function to get information about an edge of the
    // global marker graph. Returns an empty vector if the specified
    // edge does not exist.
    class GlobalMarkerGraphEdgeInformation {
    public:
        ReadId readId;
        Strand strand;
        uint32_t ordinal0;
        uint32_t ordinal1;
        uint32_t position0;
        uint32_t position1;
        uint32_t overlappingBaseCount;
        string sequence;
    };
    vector<GlobalMarkerGraphEdgeInformation> getGlobalMarkerGraphEdgeInformation(
        MarkerGraph::VertexId,
        MarkerGraph::VertexId
        );

private:
    // Lower-level, more efficient version of the above
    // (but it returns less information).
    void getGlobalMarkerGraphEdgeInfo(
        MarkerGraph::VertexId,
        MarkerGraph::VertexId,
        vector<MarkerInterval>&
        );
public:



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
    MemoryMapped::Object<AssemblerInfo> assemblerInfo;



    /***************************************************************************

    The reads used for this assembly.
    Indexed by ReadId.

    We use a run-length representation
    (https://en.wikipedia.org/wiki/Run-length_encoding)
    for reads: all repeated bases are removed, and
    for each base we store a repeat base count that says how many
    times that base was repeated in the original read.
    Many assembly phases use only the run-length representation
    (without using the base repeat count).
    This includes the generation of markers, the computation of
    alignments, and the creation of the marker graph.

    For example, suppose we have the following read sequence:

    TAATCATTTTGATGTAAGTCTAAAAATTTCACCTTAATACTTATTTTTCC

    The read is stored like this:

    TATCATGATGTAGTCTATCACTATACTATC
    121114111112111153112221112152

    The first line is stored in reads[readId] and the second line
    is stored in readRepeatCount[readId]. Note that base caller errors
    in the number of times a base is repeated (a g. AAAA versus AAAAA)
    leave the first line unchanged.

    In this representation the sequence never has any repeated bases.
    However, repeats with period 2 or longer are possible, for example TATA
    above.

    In the run-length representation, the read sequence, which has all
    repeated bases removed, is insensitive to base caller errors
    in the number of repetitions of individual bases, which is the
    most common type of error in nanopore sequencing.
    Therefore, it is hoped that the run-length representation results
    in better resilience to base caller errors.

    The base repeat count is stored in one byte per base and so
    it can store base repeat lengths up to 255.
    If a read contains bases repeated 256 or more times,
    it cannot be stored and is discarded on input.
    The stored base repeat count is never zero, and is one
    for most bases.

    The run-length representation is typically around 25% shorter than
    the raw representation, due to the removal of repeated bases.
    This gives some performance benefits for assembly phases that
    don't use the base repeat counts. However, the need to store
    repeat base counts (8 bits per base) increases the memory
    requirement from 2 to 10 bits per base. So overall the
    run-length representation requires more memory for the reads
    than the raw representation.

    Run-length representations that are more economic in memory are possible,
    at the price of additional code complexity and performance cost
    in assembly phases that use the base repeat counts.

    ***************************************************************************/

    LongBaseSequences reads;
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t> readRepeatCounts;
public:
    ReadId readCount() const
    {
        return ReadId(reads.size());
    }
private:
    void checkReadsAreOpen() const;
    void checkReadNamesAreOpen() const;
    void checkReadMetaDataAreOpen() const;
    void checkReadId(ReadId) const;



    // Return a base of an oriented read.
    Base getOrientedReadBase(
        OrientedReadId orientedReadId,
        uint32_t position)
    {
        const auto& read = reads[orientedReadId.getReadId()];
        if(orientedReadId.getStrand() == 0) {
            return read[position];
        } else {
            return read[read.baseCount-1-position].complement();
        }
    }



    // Same as above, but also returns the repeat count.
    pair<Base, uint8_t> getOrientedReadBaseAndRepeatCount(
        OrientedReadId orientedReadId,
        uint32_t position)
    {

        // Extract the read id and strand.
        const ReadId readId = orientedReadId.getReadId();
        const Strand strand = orientedReadId.getStrand();

        // Access the bases and repeat counts for this read.
        const auto& read = reads[readId];
        const auto& counts = readRepeatCounts[readId];

        // Compute the position as stored, depending on strand.
        uint32_t orientedPosition = position;
        if(strand == 1) {
            orientedPosition = uint32_t(read.baseCount) - 1 - orientedPosition;
        }

        // Extract the base and repeat count at this position.
        pair<Base, uint8_t> p = make_pair(read[orientedPosition], counts[orientedPosition]);

        // Complement the base, if necessary.
        if(strand == 1) {
            p.first = p.first.complement();
        }

        return p;
    }



    // Return a vector containing the raw sequence of an oriented read.
    vector<Base> getOrientedReadRawSequence(OrientedReadId);

    // Return the length of the raw sequence of a read.
    // If using the run-length representation of reads, this counts each
    // base a number of times equal to its repeat count.
    size_t getReadRawSequenceLength(ReadId);

    // Get a vector of the raw read positions
    // corresponding to each position in the run-length
    // representation of an oriented read.
    vector<uint32_t> getRawPositions(OrientedReadId) const;

    // The names of the reads from the input fasta or fastq files.
    // Indexed by ReadId.
    // Note that we don't enforce uniqueness of read names.
    // We don't use read names to identify reads.
    // These names are only used as an aid in tracing each read
    // back to its origin.
    MemoryMapped::VectorOfVectors<char, uint64_t> readNames;

    // Read meta data. This is the information following the read name
    // in the header line for fasta and fastq files.
    // Indexed by ReadId.
    MemoryMapped::VectorOfVectors<char, uint64_t> readMetaData;

    // Return a meta data field for a read, or an empty string
    // if that field is missing. This treats the meta data
    // as a space separated sequence of Key=Value,
    // without embedded spaces in each Key=Value pair.
    span<char> getMetaData(ReadId, const string& key);

    // Function to write a read in Fasta format.
    void writeRead(ReadId, ostream&);
    void writeOrientedRead(OrientedReadId, ostream&);
    void writeOrientedRead(OrientedReadId, const string& fileName);

    // Write a csv file with summary information for each read.
public:
    void writeReadsSummary();



    // Read flags.
private:
    MemoryMapped::Vector<ReadFlags> readFlags;
public:
    void initializeReadFlags();
    void accessReadFlags(bool readWriteAccess);
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

    // Inverse of the above: given a global marker id,
    // return its OrientedReadId and ordinal.
    // This requires a binary search in the markers toc.
    // This could be avoided, at the cost of storing
    // an additional 4 bytes per marker.
public:
    pair<OrientedReadId, uint32_t> findMarkerId(MarkerId) const;
private:

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
    void writeAlignmentCandidates() const;
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
    vector<OrientedReadPair> getAlignmentCandidates() const;
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
        Alignment&,
        AlignmentInfo&);


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
    
    void checkAlignmentDataAreOpen();
public:
    void accessCompressedAlignments();
private:

    // The alignment table stores the AlignmentData that each oriented read is involved in.
    // Stores, for each OrientedReadId, a vector of indexes into the alignmentData vector.
    // Indexed by OrientedReadId::getValue(),
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> alignmentTable;
    void computeAlignmentTable();



    // Private functions and data used by computeAlignments.
    void computeAlignmentsThreadFunction(size_t threadId);
    class ComputeAlignmentsData {
    public:

        // Parameters.
        size_t alignmentMethod;
        uint32_t maxMarkerFrequency;
        size_t maxSkip;
        size_t minAlignedMarkerCount;
        double minAlignedFraction;
        size_t maxTrim;
        size_t maxDrift;
        int matchScore;
        int mismatchScore;
        int gapScore;
        double downsamplingFactor;
        int bandExtend;
        bool suppressContainments;
        bool storeAlignments;

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
        findOrientedAlignments(OrientedReadId) const;



    // Read graph and related functions and data.
    // For more information, see comments in ReadGraph.hpp.
    ReadGraph readGraph;
public:
    void createReadGraph(
        uint32_t maxAlignmentCount,
        uint32_t maxTrim);
    void accessReadGraph();
    void accessReadGraphReadWrite();
    void checkReadGraphIsOpen();


    void flagCrossStrandReadGraphEdges(int maxDistance, size_t threadCount);
private:
    void flagCrossStrandReadGraphEdgesThreadFunction(size_t threadId);
    class FlagCrossStrandReadGraphEdgesData {
    public:
        size_t maxDistance;
        vector<bool> isNearStrandJump;
    };
    FlagCrossStrandReadGraphEdgesData flagCrossStrandReadGraphEdgesData;



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
        size_t maxTrim,         // Used to define containment.
        double timeout,         // Or 0 for no timeout.
        LocalReadGraph&);

public:



    // Directed version of the read graph.
    // This is only used if--ReadGraph.creationMethod is 1
    // If --ReadGraph.creationMethod is 0, this is not used
    // and instead we use the undirected read graph defined above.
    DirectedReadGraph directedReadGraph;
    void createDirectedReadGraph(
        uint64_t maxTrim,
        uint64_t containedNeighborCount,
        uint64_t uncontainedNeighborCountPerDirection);
    void accessDirectedReadGraphReadOnly();
    void accessDirectedReadGraphReadWrite();
    void markDirectedReadGraphConflictEdges1();
    void markDirectedReadGraphConflictEdges2(int radius);
    void markDirectedReadGraphConflictEdges3(int radius);



    // Write a FASTA file containing all reads that appear in
    // the local read graph.
    void writeLocalReadGraphReads(
        ReadId,
        Strand,
        uint32_t maxDistance,
        bool allowChimericReads,
        bool allowCrossStrandEdges);


    // Compute connected components of the read graph.
    // This treats chimeric reads as isolated.
    // Components with fewer than minComponentSize are considered
    // small and excluded from assembly by setting the
    // isInSmallComponent for all the reads they contain.
    void computeReadGraphConnectedComponents(size_t minComponentSize);



    // Private functions and data used by createMarkerGraphVertices.
private:
    void createMarkerGraphVerticesThreadFunction1(size_t threadId);
    void createMarkerGraphVerticesThreadFunction2(size_t threadId);
    void createMarkerGraphVerticesThreadFunction3(size_t threadId);
    void createMarkerGraphVerticesThreadFunction4(size_t threadId);
    void createMarkerGraphVerticesThreadFunction5(size_t threadId);
    void createMarkerGraphVerticesThreadFunction45(int);
    void createMarkerGraphVerticesThreadFunction6(size_t threadId);
    void createMarkerGraphVerticesThreadFunction7(size_t threadId);
    class CreateMarkerGraphVerticesData {
    public:

        // Parameters.
        int alignMethod;
        size_t maxSkip;
        size_t maxDrift;
        int matchScore;
        int mismatchScore;
        int gapScore;
        double downsamplingFactor;
        int bandExtend;
        int readGraphCreationMethod;
        uint32_t maxMarkerFrequency;


        // The total number of oriented markers.
        uint64_t orientedMarkerCount;

        // Disjoint sets data structures.
        // MemoryMapped::Vector< std::atomic<DisjointSets::Aint> > disjointSetsData;
        MemoryMapped::Vector<DisjointSets::Aint> disjointSetsData;
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



    void checkMarkerGraphVerticesAreAvailable();

    // Check for consistency of globalMarkerGraphVertex and globalMarkerGraphVertices.
    void checkMarkerGraphVertices(
        size_t minCoverage,
        size_t maxCoverage);



    // Marker graph.
public:
    MarkerGraph markerGraph;

    // Find the reverse complement of each marker graph vertex.
    void findMarkerGraphReverseComplementVertices(size_t threadCount);
    void accessMarkerGraphReverseComplementVertex();
    void removeMarkerGraphVertices();
private:
    void findMarkerGraphReverseComplementVerticesThreadFunction1(size_t threadId);
    void findMarkerGraphReverseComplementVerticesThreadFunction2(size_t threadId);



    // Create marker graph edges.
public:
    void createMarkerGraphEdges(size_t threadCount);
    void accessMarkerGraphEdges(bool accessEdgesReadWrite);
    void checkMarkerGraphEdgesIsOpen();
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

    // Analyze a vertex of the Marker graph.
    void analyzeMarkerGraphVertex(MarkerGraph::VertexId) const;

    // Refine the marker graph by removing vertices in tangle regions,
    // then recreating edges. This must be called after
    // transitive reduction. After this is called, the only
    // two MarkerGraph field filled in are vertices and vertexTable.
    // Everything else has to be recreated.
    void refineMarkerGraph(
        uint64_t refineThreshold,
        size_t threadCount);

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

    // Find the children or parents of a vertex of the global marker graph.
    void getGlobalMarkerGraphVertexChildren(
        MarkerGraph::VertexId,
        vector<MarkerGraph::VertexId>&,
        bool append = false
        ) const;
    void getGlobalMarkerGraphVertexParents(
        MarkerGraph::VertexId,
        vector<MarkerGraph::VertexId>&,
        bool append = false
        ) const;

    // This version also returns the oriented read ids and ordinals
    // that caused a child to be marked as such.
    void getGlobalMarkerGraphVertexChildren(
        MarkerGraph::VertexId,
        vector< pair<MarkerGraph::VertexId, vector<MarkerInterval> > >&,
        vector< pair<MarkerGraph::VertexId, MarkerInterval> >& workArea
        ) const;
    void getGlobalMarkerGraphVertexParents(
        MarkerGraph::VertexId,
        vector< pair<MarkerGraph::VertexId, vector<MarkerInterval> > >&,
        vector< pair<MarkerGraph::VertexId, MarkerInterval> >& workArea
        ) const;

    // Return true if a vertex of the global marker graph has more than
    // one marker for at least one oriented read id.
    bool isBadMarkerGraphVertex(MarkerGraph::VertexId) const;

    // Find out if a vertex is a forward or backward leaf of the pruned
    // strong subgraph of the marker graph.
    // A forward leaf is a vertex with out-degree 0.
    // A backward leaf is a vertex with in-degree 0.
    bool isForwardLeafOfMarkerGraphPrunedStrongSubgraph(MarkerGraph::VertexId) const;
    bool isBackwardLeafOfMarkerGraphPrunedStrongSubgraph(MarkerGraph::VertexId) const;

    // Given an edge of the pruned strong subgraph of the marker graph,
    // return the next/previous edge in the linear chain the edge belongs to.
    // If the edge is the last/first edge in its linear chain, return MarkerGraph::invalidEdgeId.
    MarkerGraph::EdgeId nextEdgeInMarkerGraphPrunedStrongSubgraphChain(MarkerGraph::EdgeId) const;
    MarkerGraph::EdgeId previousEdgeInMarkerGraphPrunedStrongSubgraphChain(MarkerGraph::EdgeId) const;

    // Return the out-degree or in-degree (number of outgoing/incoming edges)
    // of a vertex of the pruned strong subgraph of the marker graph.
    size_t markerGraphPrunedStrongSubgraphOutDegree(MarkerGraph::VertexId) const;
    size_t markerGraphPrunedStrongSubgraphInDegree (MarkerGraph::VertexId) const;

    // Return true if an edge disconnects the local subgraph.
    bool markerGraphEdgeDisconnectsLocalStrongSubgraph(
        MarkerGraph::EdgeId edgeId,
        size_t maxDistance,

        // Work areas, to reduce memory allocation activity.

        // Each of these two must be sized maxDistance+1.
        array<vector< vector<MarkerGraph::EdgeId> >, 2>& verticesByDistance,

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
        vector<MarkerGraph::EdgeId>& path,
        vector< pair<uint32_t, uint32_t> >& pathOrdinals
        ) const;



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

    // Evaluate an induced alignment.
    // Contrary to InducedAlignment::evaluate, this takes into account
    // markers that don't correspond to a marker graph vertex.
    bool evaluateInducedAlignment(
        OrientedReadId orientedReadId0,
        OrientedReadId orientgedReadId1,
        const InducedAlignment&,
        const InducedAlignmentCriteria&,
        vector<uint64_t>& work);



    // Conflict read graph.
    // See ConflictReadGraph.hpp for more information.
public:
    void createConflictReadGraph(
        uint64_t threadCount,
        uint32_t maxOffsetSigma,
        uint32_t maxTrim,
        uint32_t maxSkip,
        uint32_t minAlignedMarkerCount);
    void accessConflictReadGraph();
    // void colorConflictReadGraph();
    void cleanupConflictReadGraph();
private:
    // void createConflictReadGraphThreadFunction1(size_t threadId);
    void createConflictReadGraphThreadFunction2(size_t threadId);
    void addConflictGraphEdges(
        ReadId,
        const InducedAlignmentCriteria&,
        // Work areas.
        vector<OrientedReadId>&,
        vector<OrientedReadId>&,
        vector<InducedAlignment>&,
        vector<uint64_t>&
        );
    class CreateConflictReadGraphData {
    public:
        InducedAlignmentCriteria inducedAlignmentCriteria;
    };
    CreateConflictReadGraphData createConflictReadGraphData;
    ConflictReadGraph conflictReadGraph;



    // Class used by colorConflictReadGraph.
    // It creates a priority_queue containing this type.
    class ColorConflictReadGraphData {
    public:
        DirectedReadGraph::VertexId vertexId;
        uint64_t conflictReadGraphDegree;
        uint64_t directedReadGraphKeptDegree;

        ColorConflictReadGraphData(
            DirectedReadGraph::VertexId,
            const DirectedReadGraph&,
            const ConflictReadGraph&);

        // Used for sorting, so values with low conflict degree come first.
        bool operator <(const ColorConflictReadGraphData& that) const
        {
            return
                conflictReadGraphDegree < that.conflictReadGraphDegree
                or
                (
                    conflictReadGraphDegree == that.conflictReadGraphDegree
                    and
                    directedReadGraphKeptDegree > that.directedReadGraphKeptDegree
                );
        }

        // Used for sorting, so values with low conflict degree
        // are at the top of the queue.
        bool operator>(const ColorConflictReadGraphData& that) const
        {
            return
                conflictReadGraphDegree > that.conflictReadGraphDegree
                or
                (
                    conflictReadGraphDegree == that.conflictReadGraphDegree
                    and
                    directedReadGraphKeptDegree < that.directedReadGraphKeptDegree
                );
        }
    };



#ifdef SHASTA_HTTP_SERVER
    // Extract a local subgraph of the global marker graph.
    bool extractLocalMarkerGraphUsingStoredConnectivity(
        OrientedReadId,
        uint32_t ordinal,
        int distance,
        int timeout,                 // Or 0 for no timeout.
        bool useWeakEdges,
        bool usePrunedEdges,
        bool useSuperBubbleEdges,
        LocalMarkerGraph&
        );
    bool extractLocalMarkerGraphUsingStoredConnectivity(
        MarkerGraph::VertexId,
        int distance,
        int timeout,                 // Or 0 for no timeout.
        bool useWeakEdges,
        bool usePrunedEdges,
        bool useSuperBubbleEdges,
        LocalMarkerGraph&
        );
#endif

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



    // Use spoa to compute consensus sequence for an edge of the marker graph.
    // This does not include the bases corresponding to the flanking markers.
    void computeMarkerGraphEdgeConsensusSequenceUsingSpoa(
        MarkerGraph::EdgeId,
        uint32_t markerGraphEdgeLengthThresholdForConsensus,
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
    void createAssemblyGraphVertices();
    void accessAssemblyGraphVertices();
    void createAssemblyGraphEdges();
    void accessAssemblyGraphEdgeLists();
    void accessAssemblyGraphEdges();
    void accessAssemblyGraphOrientedReadsByEdge();
    void writeAssemblyGraph(const string& fileName) const;
    void findAssemblyGraphBubbles();

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
        AssemblyGraph::EdgeId,
        int distance,
        double timeout,
        LocalAssemblyGraph&) const;
public:
    void colorGfaBySimilarityToSegment(
        AssemblyGraph::EdgeId,
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
        const vector<MarkerGraph::EdgeId>& markerGraphEdges,
        vector<AssemblyGraph::EdgeId>& assemblyGraphEdges
        ) const;



    // Pseudo-paths.
    // An oriented read corresponds to a path (sequence of adjacent edges)
    // in the marker graph, which
    // can be computed via computeOrientedReadMarkerGraphPath.
    // That path encounters a sequence of assembly graph edges,
    // which is not necessarily a path in the assembly graph
    // because not all marker graph edges belong to an assembly graph edge.
    // We call this sequence the pseudo-path of an oriented read in the assembly graph.
    class PseudoPathEntry {
    public:
        AssemblyGraph::EdgeId segmentId;

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
        vector<MarkerGraph::EdgeId>& path,
        vector< pair<uint32_t, uint32_t> >& pathOrdinals,

        // The pseudo-path computed by this function.
        PseudoPath&) const;
public:
    void writePseudoPath(ReadId, Strand) const;



    // Detangle the AssemblyGraph.
    void detangle();



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
        bool storeCoverageData
        );
private:
    void assembleMarkerGraphEdgesThreadFunction(size_t threadId);
    class AssembleMarkerGraphEdgesData {
    public:

        // The arguments to assembleMarkerGraphEdges, stored here so
        // they are accessible to the threads.
        uint32_t markerGraphEdgeLengthThresholdForConsensus;
        bool storeCoverageData;

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
        vector< shared_ptr< MemoryMapped::Vector<MarkerGraph::EdgeId> > > threadEdgeIds;
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
    // Optionally outputs detailed assembly information
    // in html (skipped if the html pointer is 0).
    void assembleAssemblyGraphEdge(
        AssemblyGraph::EdgeId,
        bool storeCoverageData,
        AssembledSegment&);
public:
    AssembledSegment assembleAssemblyGraphEdge(
        AssemblyGraph::EdgeId,
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
        vector< vector<AssemblyGraph::EdgeId> > edges;
        vector< shared_ptr<LongBaseSequences> > sequences;
        vector< shared_ptr<MemoryMapped::VectorOfVectors<uint8_t, uint64_t> > > repeatCounts;
        void allocate(size_t threadCount);
        void free();
    };
    AssembleData assembleData;
    void assembleThreadFunction(size_t threadId);



    // Write the assembly graph in GFA 1.0 format defined here:
    // https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
public:
    void writeGfa1(const string& fileName);
    void writeGfa1BothStrands(const string& fileName);
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



    // The segment graph is used to detangle and determine reachability
    // in the assembly graph.
    void createSegmentGraph();
    void colorGfaBySegmentGraphChain(uint64_t chainId) const;
    shared_ptr<SegmentGraph> segmentGraphPointer;



    // Analyze oriented read paths in the marker graph and in the assembly graph.
    void analyzeOrientedReadPaths(int readGraphCreationMethod) const;

    // Analyze paths of oriented reads that go through a given assembly graph edge (segment).
    void analyzeOrientedReadPathsThroughSegment(
        AssemblyGraph::EdgeId segmentId) const;




    // Data and functions used for the http server.
    // This function puts the server into an endless loop
    // of processing requests.
    void writeHtmlBegin(ostream&, bool navigation=true) const;
    void writeHtmlEnd(ostream&) const;
    void writeAssemblySummary(ostream&);
    void writeAssemblySummaryBody(ostream&);
    void writeAssemblySummaryJson(ostream&);
    void writeAssemblyIndex(ostream&) const;
    void writeNavigation(ostream&) const;
    void writeNavigation(
        ostream& html,
        const string& title,
        const vector<pair <string, string> >&) const;
    static void writeStyle(ostream& html);
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

#ifdef SHASTA_HTTP_SERVER

    void fillServerFunctionTable();
    void processRequest(
        const vector<string>& request,
        ostream&,
        const BrowserInformation&) override;
    void writeMakeAllTablesSelectable(ostream&) const;
    void exploreSummary(const vector<string>&, ostream&);
    void exploreRead(const vector<string>&, ostream&);
    void blastRead(const vector<string>&, ostream&);
    void exploreAlignments(const vector<string>&, ostream&);
    void exploreAlignment(const vector<string>&, ostream&);
    void displayAlignmentMatrix(const vector<string>&, ostream&);
    void exploreAlignmentGraph(const vector<string>&, ostream&);
    void exploreReadGraph(const vector<string>&, ostream&);
    void exploreUndirectedReadGraph(const vector<string>&, ostream&);
    void exploreDirectedReadGraph(const vector<string>&, ostream&);
    void exploreCompressedAssemblyGraph(const vector<string>&, ostream&);
    class HttpServerData {
    public:

        using ServerFunction = void (Assembler::*) (
            const vector<string>& request,
            ostream&);
        std::map<string, ServerFunction> functionTable;
        string docsDirectory;
        string referenceFastaFileName = "reference.fa";

        const AssemblerOptions* assemblerOptions = 0;

    };
    HttpServerData httpServerData;

    // Display alignments in an html table.
    void displayAlignments(
        OrientedReadId,
        const vector< pair<OrientedReadId, AlignmentInfo> >&,
        ostream&) const;
    void displayAlignment(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        const AlignmentInfo& alignment,
        ostream&) const;

    // Functions and data used by the http server
    // for display of the local marker graph.
    void exploreMarkerGraph(const vector<string>&, ostream&);
    class LocalMarkerGraphRequestParameters {
    public:
        MarkerGraph::VertexId vertexId;
        bool vertexIdIsPresent;
        uint32_t maxDistance;
        bool maxDistanceIsPresent;
        bool addLabels;
        bool useDotLayout;  // If true, use dot. If false, use sfdp.
        bool useWeakEdges;
        bool usePrunedEdges;
        bool useSuperBubbleEdges;
        uint32_t sizePixels;
        bool sizePixelsIsPresent;
        double vertexScalingFactor;
        bool vertexScalingFactorIsPresent;
        string vertexScalingFactorString() const;
        double arrowScalingFactor;
        bool arrowScalingFactorIsPresent;
        string arrowScalingFactorString() const;
        int timeout;
        bool timeoutIsPresent;
        void writeForm(ostream&, MarkerGraph::VertexId vertexCount) const;
        bool hasMissingRequiredParameters() const;
    };
    void getLocalMarkerGraphRequestParameters(
        const vector<string>&,
        LocalMarkerGraphRequestParameters&) const;
    void exploreMarkerGraphVertex(const vector<string>&, ostream&);
    void exploreMarkerGraphEdge(const vector<string>&, ostream&);
    void exploreMarkerCoverage(const vector<string>&, ostream&);
    void exploreMarkerGraphInducedAlignment(const vector<string>&, ostream&);
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
        const uint32_t bandExtend,
        ostream& html
    );
#endif



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
        uint32_t bandExtend;
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



    // Functions and data used by the http server
    // for display of the local assembly graph.
private:
    void exploreAssemblyGraph(const vector<string>&, ostream&);
    class LocalAssemblyGraphRequestParameters {
    public:
        AssemblyGraph::EdgeId edgeId;
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
        bool highlightedReadIdIsPresent;
        ReadId highlightedReadId;
        Strand highlightedStrand;
        void writeForm(ostream&, AssemblyGraph::EdgeId edgeCount, bool allowHighlighting) const;
        bool hasMissingRequiredParameters() const;
    };
    void getLocalAssemblyGraphRequestParameters(
        const vector<string>&,
        LocalAssemblyGraphRequestParameters&) const;
    void exploreAssemblyGraphEdge(const vector<string>&, ostream&);
    void exploreAssemblyGraphEdgesSupport(const vector<string>&, ostream&);



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



    // Functions and data structures used for phasing.
public:
    void createPhasingData(
        size_t threadCount,
        double phasingSimilarityThreshold,
        int maxNeighborCount);
    void accessPhasingData();
    // double computePhasingSimilarity(OrientedReadId, OrientedReadId);
    // double computePhasingSimilarity(ReadId, Strand, ReadId, Strand);
    double computePhasingSimilarity(AssemblyGraph::EdgeId, AssemblyGraph::EdgeId);
    uint64_t countCommonInternalOrientedReads(AssemblyGraph::EdgeId, AssemblyGraph::EdgeId);
private:
    PhasingData phasingData;
    void phasingGatherOrientedReads(size_t threadCount);
    void phasingGatherOrientedReadsPass1(size_t threadId);
    void phasingGatherOrientedReadsPass2(size_t threadId);
    void phasingGatherOrientedReadsPass(int pass);
    void phasingGatherAssemblyGraphEdges(size_t threadCount);
    void phasingGatherAssemblyGraphEdgesPass1(size_t threadId);
    void phasingGatherAssemblyGraphEdgesPass2(size_t threadId);
    void phasingGatherAssemblyGraphEdgesPass(int pass);
    void phasingSortAssemblyGraphEdges(size_t threadCount);
    void phasingSortAssemblyGraphEdgesThreadFunction(size_t threadId);
    void phasingWriteBipartiteGraph();
    void phasingFindSimilarForks();

public:
    void test();
};

#endif
