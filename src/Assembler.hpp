#ifndef CZI_SHASTA_ASSEMBLER_HPP
#define CZI_SHASTA_ASSEMBLER_HPP

// Shasta.
#include "Alignment.hpp"
#include "AssembledSegment.hpp"
#include "AssemblyGraph.hpp"
#include "Coverage.hpp"
#include "dset64.hpp"
#include "HttpServer.hpp"
#include "Kmer.hpp"
#include "LongBaseSequence.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "MemoryMappedObject.hpp"
#include "MultitreadedObject.hpp"
#include "OrientedReadPair.hpp"
#include "ReadGraph.hpp"
#include "ReadFlags.hpp"
#include "ReadId.hpp"

#ifndef SHASTA_STATIC_EXECUTABLE
// MarginPhase.
#include "marginPhase/callConsensus.h"
#endif

// Standard library.
#include "memory.hpp"
#include "string.hpp"
#include "tuple.hpp"

namespace ChanZuckerberg {
    namespace shasta {

        // Forward declarations of classes defined here.
        class Assembler;
        class AssemblerInfo;

        // Forward declarations of classes defined elsewhere.
        class Alignment;
        class AlignmentGraph;
        class AlignmentInfo;
        class AssembledSegment;
        class ConsensusCaller;
        class LocalAlignmentGraph;
        class LocalAssemblyGraph;
#ifndef SHASTA_STATIC_EXECUTABLE
        class LocalMarkerGraph;
#endif
        class LocalReadGraph;
        class MarkerInterval;
        namespace MemoryMapped {
            template<class Int, class T> class VectorOfVectors;
        }

        // Write an html form to select strand.
        void writeStrandSelection(
            ostream&,               // The html stream to write the form to.
            const string& name,     // The selection name.
            bool select0,           // Whether strand 0 is selected.
            bool select1);          // Whether strand 1 is selected.
    }
}



// Class used to store various pieces of assembler information in shared memory.
class ChanZuckerberg::shasta::AssemblerInfo {
public:

    // The length of k-mers used to define markers.
    size_t k;

    // The page size in use for this run.
    size_t largeDataPageSize;



    // Statistics on the number of reads discarded on input.
    // These are incremented during each call to addReadsFromFasta.

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



class ChanZuckerberg::shasta::Assembler :
    public MultithreadedObject<Assembler>
#ifndef SHASTA_STATIC_EXECUTABLE
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

    // Destructor.
    ~Assembler();

    // Add reads from a fasta file.
    // The reads are added to those already previously present.
    void addReadsFromFasta(
        const string& fileName,
        size_t minReadLength,
        size_t blockSize,
        size_t threadCountForReading,
        size_t threadCountForProcessing);

    // Create a histogram of read lengths.
    void histogramReadLength(const string& fileName);

    // Function to write one or all reads in Fasta format.
    void writeReads(const string& fileName);
    void writeRead(ReadId, const string& fileName);
    void writeOrientedRead(ReadId, Strand, const string& fileName);

    // Functions related to the k-mer table.
    void accessKmers();
    void writeKmers(const string& fileName) const;
    void randomlySelectKmers(
        size_t k,           // k-mer length.
        double probability, // The probability that a k-mer is selected as a marker.
        int seed            // For random number generator.
    );

    // Functions related to markers.
    // See the beginning of Marker.hpp for more information.
    void findMarkers(size_t threadCount);
    void accessMarkers();
    void writeMarkers(ReadId, Strand, const string& fileName);

    // Use the minHash algorithm to find candidate alignments.
    // Use as features sequences of m consecutive special k-mers.
    void findAlignmentCandidatesMinHash(
        size_t m,                       // Number of consecutive k-mers that define a feature.
        size_t minHashIterationCount,   // Number of minHash iterations.
        size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
        size_t maxBucketSize,           // The maximum size for a bucket to be used.
        size_t minFrequency,            // Minimum number of minHash hits for a pair to become a candidate.
        size_t threadCount
    );
    void findAlignmentCandidatesLowHash(
        size_t m,                       // Number of consecutive k-mers that define a feature.
        double hashFraction,            // Low hash threshold.
        size_t minHashIterationCount,   // Number of lowHash iterations.
        size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for lowHash.
        size_t maxBucketSize,           // The maximum size for a bucket to be used.
        size_t minFrequency,            // Minimum number of lowHash hits for a pair to become a candidate.
        size_t threadCount
    );
    void accessAlignmentCandidates();

    // Write the reads that overlap a given read.
    void writeOverlappingReads(ReadId, Strand, const string& fileName);

    // Compute a marker alignment of two oriented reads.
    void alignOrientedReads(
        ReadId, Strand,
        ReadId, Strand,
        size_t maxSkip, // Maximum ordinal skip allowed.
        uint32_t maxMarkerFrequency
    );

    // Compute marker alignments of an oriented read with all reads
    // for which we have an Overlap.
    void alignOverlappingOrientedReads(
        ReadId, Strand,
        size_t maxSkip,                 // Maximum ordinal skip allowed.
        uint32_t maxMarkerFrequency,
        size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
        size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
    );



#ifndef SHASTA_STATIC_EXECUTABLE
    // Extract a local marker graph from the global marker graph.
    void extractLocalMarkerGraph(

        // The ReadId, Strand, and ordinal that identify the
        // marker corresponding to the start vertex
        // for the local marker graph to be created.
        ReadId,
        Strand,
        uint32_t ordinal,

        // Maximum distance from the start vertex (number of edges in the global marker graph).
        int distance,

        // Minimum coverage for a strong vertex or edge (affects coloring).
        size_t minCoverage
        );
#endif



    // Compute an alignment for each alignment candidate.
    // Store summary information for the ones that are good enough,
    // without storing details of the alignment.
    void computeAlignments(

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

        // Minimum number of alignment markers for an alignment to be used.
        size_t minAlignedMarkerCount,

        // Maximum left/right trim (in bases) for an alignment to be used.
        size_t maxTrim,

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

        // The maximum frequency of marker k-mers to be used in
        // computing alignments.
        uint32_t maxMarkerFrequency,

        // The maximum ordinal skip to be tolerated between successive markers
        // in the alignment.
        size_t maxSkip,

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
    void accessMarkerGraphVertices();

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


#ifndef SHASTA_STATIC_EXECUTABLE
    // Create a local marker graph and return its local assembly path.
    // The local marker graph is specified by its start vertex
    // and maximum distance (number of edges) form the start vertex.
    vector<MarkerGraph::VertexId> getLocalAssemblyPath(
        MarkerGraph::VertexId,
        int maxDistance
        );
#endif

    // Find weak edges in the marker graph.
    void flagMarkerGraphWeakEdges(
        size_t lowCoverageThreshold,
        size_t highCoverageThreshold,
        size_t maxDistance,
        size_t edgeMarkerSkipThreshold);



    // Call this before explore to make the documentation available.
    void setDocsDirectory(const string&);

    // Call this before explore to specify the name of the fasta
    // file containing the reference to be used with Blast commands.
    void setReferenceFastaFileName(const string&);



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

    // Function to write a read in Fasta format.
    void writeRead(ReadId, ostream&);
    void writeOrientedRead(OrientedReadId, ostream&);
    void writeOrientedRead(OrientedReadId, const string& fileName);



    // Read flags.
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
    pair<OrientedReadId, uint32_t> findMarkerId(MarkerId) const;

    // Given a MarkerId, compute the MarkerId of the
    // reverse complemented marker.
    MarkerId findReverseComplement(MarkerId) const;


    // Flag palindromic reads.
public:
    void flagPalindromicReads(
        uint32_t maxSkip,
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
        uint32_t maxMarkerFrequency;
        double alignedFractionThreshold;
        double nearDiagonalFractionThreshold;
        uint32_t deltaThreshold;
    };
    FlagPalindromicReadsData flagPalindromicReadsData;


    // Alignment candidate found by the MinHash algorithm.
    // They all have readId0<readId1.
    MemoryMapped::Vector<OrientedReadPair> alignmentCandidates;
    void checkAlignmentCandidatesAreOpen() const;



    // Compute a marker alignment of two oriented reads.
    void alignOrientedReads(
        OrientedReadId,
        OrientedReadId,
        size_t maxSkip, // Maximum ordinal skip allowed.
        uint32_t maxMarkerFrequency
    );
    // This lower level version takes as input vectors of
    // markers already sorted by kmerId.
    void alignOrientedReads(
        const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
        size_t maxSkip,  // Maximum ordinal skip allowed.
        uint32_t maxMarkerFrequency
    );
    // This version allows reusing the AlignmentGraph and Alignment
    void alignOrientedReads(
        const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
        size_t maxSkip,             // Maximum ordinal skip allowed.
        uint32_t maxMarkerFrequency,
        bool debug,
        AlignmentGraph&,
        Alignment&,
        AlignmentInfo&
    );

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
        uint32_t maxMarkerFrequency,
        size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
        size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
    );



    // The good alignments we found.
    // They are stored with readId0<readId1 and with strand0==0.
    MemoryMapped::Vector<AlignmentData> alignmentData;
    void checkAlignmentDataAreOpen();

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
        uint32_t maxMarkerFrequency;
        size_t maxSkip;
        size_t minAlignedMarkerCount;
        size_t maxTrim;

        // The AlignmentInfo found by each thread.
        vector< vector<AlignmentData> > threadAlignmentData;
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
#if 1
    void createReadGraphNew(
        uint32_t maxAlignmentCount,
        uint32_t maxTrim);
#endif
    void accessReadGraph();
    void accessReadGraphReadWrite();
    void checkReadGraphIsOpen();


    void flagCrossStrandReadGraphEdges();
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
        size_t maxSkip;
        uint32_t maxMarkerFrequency;

        // The total number of oriented markers.
        uint64_t orientedMarkerCount;

        // Disjoint sets data structures.
        MemoryMapped::Vector< std::atomic<DisjointSets::Aint> > disjointSetsData;
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
    MarkerGraph markerGraph;

    // Find the reverse complement of each marker graph vertex.
public:
    void findMarkerGraphReverseComplementVertices(size_t threadCount);
    void accessMarkerGraphReverseComplementVertex();
private:
    void findMarkerGraphReverseComplementVerticesThreadFunction1(size_t threadId);
    void findMarkerGraphReverseComplementVerticesThreadFunction2(size_t threadId);



    // Create marker graph edges.
public:
    void createMarkerGraphEdges(size_t threadCount);
    void accessMarkerGraphEdges(bool accessEdgesReadWrite);
    void checkMarkerGraphEdgesIsOpen();        // The edges and their MarkerIntervals found by each thread.
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


#ifndef SHASTA_STATIC_EXECUTABLE
    // Extract a local marker graph from the global marker graph.
    void extractLocalMarkerGraph(

        // The OrientedReadId and ordinal that identify the
        // marker corresponding to the start vertex
        // for the local marker graph to be created.
        OrientedReadId,
        uint32_t ordinal,

        // Maximum distance from the start vertex (number of edges).
        int distance,

        // Minimum coverage for a strong vertex.
        size_t minCoverage,

        // Minimum consensus for a strong edge.
        size_t minConsensus
        );
    bool extractLocalMarkerGraph(
        OrientedReadId,
        uint32_t ordinal,
        int distance,
        double timeout,                 // Or 0 for no timeout.
        LocalMarkerGraph&
        );
    bool extractLocalMarkerGraph(
        MarkerGraph::VertexId,
        int distance,
        double timeout,                 // Or 0 for no timeout.
        LocalMarkerGraph&
        );

    // Versions of the above that use stored connectivity of the
    // global marker graph instead of creating it on the fly.
    bool extractLocalMarkerGraphUsingStoredConnectivity(
        OrientedReadId,
        uint32_t ordinal,
        int distance,
        double timeout,                 // Or 0 for no timeout.
        bool useWeakEdges,
        bool usePrunedEdges,
        bool useShortCycleEdges,
        bool useBubbleEdges,
        bool useBubbleReplacementEdges,
        bool useSuperBubbleEdges,
        bool useSuperBubbleReplacementEdges,
        LocalMarkerGraph&
        );
    bool extractLocalMarkerGraphUsingStoredConnectivity(
        MarkerGraph::VertexId,
        int distance,
        double timeout,                 // Or 0 for no timeout.
        bool useWeakEdges,
        bool usePrunedEdges,
        bool useShortCycleEdges,
        bool useBubbleEdges,
        bool useBubbleReplacementEdges,
        bool useSuperBubbleEdges,
        bool useSuperBubbleReplacementEdges,
        LocalMarkerGraph&
        );
#endif

    // Compute consensus sequence for a vertex of the marker graph.
    void computeMarkerGraphVertexConsensusSequence(
        MarkerGraph::VertexId,
        vector<Base>& sequence,
        vector<uint32_t>& repeatCounts
        );

    // Compute consensus sequence for an edge of the marker graph.
    // This does not include the bases corresponding to the flanking markers.
#ifndef SHASTA_STATIC_EXECUTABLE
    void computeMarkerGraphEdgeConsensusSequenceUsingSeqan(
        MarkerGraph::EdgeId,
        vector<Base>& sequence,
        vector<uint32_t>& repeatCounts,
        uint8_t& overlappingBaseCount
        );
#endif
    void computeMarkerGraphEdgeConsensusSequenceUsingSpoa(
        MarkerGraph::EdgeId,
        uint32_t markerGraphEdgeLengthThresholdForConsensus,
        vector<Base>& sequence,
        vector<uint32_t>& repeatCounts,
        uint8_t& overlappingBaseCount,
        vector< pair<uint32_t, CompressedCoverageData> >* coverageData // Optional
        );
#ifndef SHASTA_STATIC_EXECUTABLE
    void computeMarkerGraphEdgeConsensusSequenceUsingMarginPhase(
        MarkerGraph::EdgeId,
        vector<Base>& sequence,
        vector<uint32_t>& repeatCounts,
        uint8_t& overlappingBaseCount
        );
#endif



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



    // In the assembly graph, each vertex corresponds to a linear chain
    // of edges in the pruned strong subgraph of the marker graph.
    // A directed vertex A->B is created if the last marker graph vertex
    // of the edge chain corresponding to A coincides with the
    // first marker graph vertex of the edge chain corresponding to B.
    AssemblyGraph assemblyGraph;
public:
    void createAssemblyGraphVertices();
    void accessAssemblyGraphVertices();
    void createAssemblyGraphEdges();
    void accessAssemblyGraphEdgeLists();
    void accessAssemblyGraphEdges();
    void writeAssemblyGraph(const string& fileName) const;
private:

    // Extract a local assembly graph from the global assembly graph.
    // This returns false if the timeout was exceeded.
    bool extractLocalAssemblyGraph(
        AssemblyGraph::EdgeId,
        int distance,
        double timeout,
        LocalAssemblyGraph&) const;
public:


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



public:
    // Assemble consensus sequence and repeat counts for each marker graph edge.
    void assembleMarkerGraphEdges(
        size_t threadCount,

        // This controls when we give up trying to compute consensus for long edges.
        uint32_t markerGraphEdgeLengthThresholdForConsensus,

        // Parameter to control whether we use spoa or marginPhase
        // to compute consensus sequence.
        bool useMarginPhase,

        // Request storing detailed coverage information.
        bool storeCoverageData
        );
    void accessMarkerGraphEdgeConsensus();
private:
    void assembleMarkerGraphEdgesThreadFunction(size_t threadId);
    class AssembleMarkerGraphEdgesData {
    public:

        // The arguments to assembleMarkerGraphEdges, stored here so
        // they are accessible to the threads.
        uint32_t markerGraphEdgeLengthThresholdForConsensus;
        bool useMarginPhase;
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
    void assemble(size_t threadCount);
    void accessAssemblyGraphSequences();
    void computeAssemblyStatistics();
private:
    class AssembleData {
    public:

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
private:
    // Construct the CIGAR string given two vectors of repeat counts.
    // Used by writeGfa1.
    static void constructCigarString(
        const MemoryAsContainer<uint8_t>& repeatCounts0,
        const MemoryAsContainer<uint8_t>& repeatCounts1,
        string&
        );

public:

    // Write assembled sequences in FASTA format.
    void writeFasta(const string& fileName);



    // Analyze assembled sequence by mapping it to
    // a given reference sequence used as ground truth.
    // The first one is python-callable.
    // The edgeId is the same as in FASTA or GFA output
    // Begin and end specify the portion of the assembled
    // segment to be analyzed and they are expressed in
    // raw coordinates (not run-length).
    // The reference sequence is given in its raw
    // representation (not run-length).
    // Assembled sequenceis mapped to reference sequence
    // using SeqAn's implementation of the Needleman-Wunsch
    // algorithm, which is very expensive in time and memory
    // when the sequences are long.
    void analyzeAssembledSequence(
        AssemblyGraph::EdgeId edgeId,   // Same as in fasta, gfa output.
        uint32_t begin,
        uint32_t end,
        const string& referenceSequence
    );
    void analyzeAssembledSequence(
        AssemblyGraph::EdgeId edgeId,
        uint32_t begin,
        uint32_t end,
        const vector<Base>& referenceSequence
    );



    // Data and functions used for the http server.
    // This function puts the server into an endless loop
    // of processing requests.
    void writeHtmlBegin(ostream&, bool navigation=true) const;
    void writeHtmlEnd(ostream&) const;
    void writeAssemblySummary(ostream&);
    void writeAssemblySummaryBody(ostream&);
    void writeNavigation(ostream&) const;
    void writeNavigation(
        ostream& html,
        const string& title,
        const vector<pair <string, string> >&) const;
    static void writeStyle(ostream& html);

#ifndef SHASTA_STATIC_EXECUTABLE

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
    class HttpServerData {
    public:

        using ServerFunction = void (Assembler::*) (
            const vector<string>& request,
            ostream&);
        std::map<string, ServerFunction> functionTable;
        string docsDirectory;
        string referenceFastaFileName = "reference.fa";

    };
    HttpServerData httpServerData;

    // Display alignments in an html table.
    void displayAlignments(
        OrientedReadId,
        const vector< pair<OrientedReadId, AlignmentInfo> >&,
        ostream&);


    // Functions and data used by the http server
    // for display of the local marker graph.
    void exploreMarkerGraph(const vector<string>&, ostream&);
    class LocalMarkerGraphRequestParameters {
    public:
        MarkerGraph::VertexId vertexId;
        bool vertexIdIsPresent;
        uint32_t maxDistance;
        bool maxDistanceIsPresent;
        bool detailed;
        bool useStoredConnectivity;
        bool useWeakEdges;
        bool usePrunedEdges;
        bool useShortCycleEdges;
        bool useBubbleEdges;
        bool useBubbleReplacementEdges;
        bool useSuperBubbleEdges;
        bool useSuperBubbleReplacementEdges;
        bool showVertexId;
        bool showOptimalSpanningTree;
        bool showAssembledSequence;
        uint32_t minCoverage;
        bool minCoverageIsPresent;
        uint32_t sizePixels;
        bool sizePixelsIsPresent;
        double timeout;
        bool timeoutIsPresent;
        string portionToDisplay;
        void writeForm(ostream&, MarkerGraph::VertexId vertexCount) const;
        bool hasMissingRequiredParameters() const;
    };
    void getLocalMarkerGraphRequestParameters(
        const vector<string>&,
        LocalMarkerGraphRequestParameters&) const;
    void showLocalMarkerGraphAlignments(
        ostream& html,
        const LocalMarkerGraph&,
        const LocalMarkerGraphRequestParameters&
        );
    void exploreMarkerGraphVertex(const vector<string>&, ostream&);
    void exploreMarkerGraphEdge(const vector<string>&, ostream&);
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
        uint32_t maxMarkerFrequency;
        size_t minAlignedMarkerCount;
        size_t maxTrim;
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
        void writeForm(ostream&, AssemblyGraph::EdgeId edgeCount) const;
        bool hasMissingRequiredParameters() const;
    };
    void getLocalAssemblyGraphRequestParameters(
        const vector<string>&,
        LocalAssemblyGraphRequestParameters&) const;
    void exploreAssemblyGraphEdge(const vector<string>&, ostream&);



    // The ConsensusCaller used to compute the "best"
    // base and repeat count at each assembly position.
    // The argument to setupConsensusCaller specifies
    // the consensus caller to be used.
    // Supported values are:
    // - SimpleConsensusCaller
    // - SimpleBayesianConsensusCaller
public:
    void setupConsensusCaller(const string&);
private:
    shared_ptr<ConsensusCaller> consensusCaller;



    // Parameters for marginPhase.
    // Read from file MarginPhase.json in the run directory.
public:
    void setupMarginPhase();
private:
    void checkMarginPhaseWasSetup();
#ifndef SHASTA_STATIC_EXECUTABLE
    PolishParams* marginPhaseParameters;
#endif
};

#endif
