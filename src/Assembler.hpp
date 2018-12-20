#ifndef CZI_SHASTA_ASSEMBLER_HPP
#define CZI_SHASTA_ASSEMBLER_HPP

// Shasta
#include "Alignment.hpp"
#include "AssemblyGraph.hpp"
#include "dset64.hpp"
#include "HttpServer.hpp"
#include "Kmer.hpp"
#include "LongBaseSequence.hpp"
#include "Marker.hpp"
#include "MarkerId.hpp"
#include "MemoryMappedObject.hpp"
#include "MultitreadedObject.hpp"
#include "OrientedReadPair.hpp"
#include "ReadId.hpp"

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
        class ConsensusCaller;
        class LocalAlignmentGraph;
        class LocalAssemblyGraph;
        class LocalMarkerGraph;
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

    // Flag for the read representation in use:
    // false: Raw reads.
    // true:  Run-length representation.
    // See comments later near Assembler::reads for more information.
    bool useRunLengthReads = false;
};



class ChanZuckerberg::shasta::Assembler :
    public MultithreadedObject<Assembler>,
    public HttpServer {
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
        size_t largeDataPageSize,
        bool useRunLengthReads
        );

    // Constructor to be called to continue an existing run.
    Assembler(
        const string& largeDataFileNamePrefix,
        size_t largeDataPageSize
        );


    // Add reads from a fasta file.
    // The reads are added to those already previously present.
    void addReadsFromFasta(
        const string& fileName,
        size_t minReadLength,
        size_t blockSize,
        size_t threadCountForReading,
        size_t threadCountForProcessing);

    // Access the reads and read names.
    void accessReadsReadOnly();
    void accessReadsReadWrite();
    void accessReadNamesReadOnly();
    void accessReadNamesReadWrite();

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
    void findAlignmentCandidates(
        size_t m,                       // Number of consecutive k-mers that define a feature.
        size_t minHashIterationCount,   // Number of minHash iterations.
        size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
        size_t maxBucketSize,           // The maximum size for a bucket to be used.
        size_t minFrequency,            // Minimum number of minHash hits for a pair to become a candidate.
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
        size_t maxVertexCountPerKmer
    );

    // Compute marker alignments of an oriented read with all reads
    // for which we have an Overlap.
    void alignOverlappingOrientedReads(
        ReadId, Strand,
        size_t maxSkip,                 // Maximum ordinal skip allowed.
        size_t maxVertexCountPerKmer,
        size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
        size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
    );



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



    // Compute an alignment for each alignment candidate.
    // Store summary information for the ones that are good enough,
    // without storing details of the alignment.
    void computeAlignments(

        // The  maximum number of vertices in the alignment graph
        // that we allow a single k-mer to generate.
        size_t alignmentMaxVertexCountPerKmer,

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

        // The  maximum number of vertices in the alignment graph
        // that we allow a single k-mer to generate.
        size_t alignmentMaxVertexCountPerKmer,

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
    GlobalMarkerGraphVertexId getGlobalMarkerGraphVertex(
        ReadId,
        Strand,
        uint32_t ordinal) const;

    // Find the markers contained in a given vertex of the global marker graph.
    // Returns the markers as tuples(read id, strand, ordinal).
    vector< tuple<ReadId, Strand, uint32_t> >
        getGlobalMarkerGraphVertexMarkers(GlobalMarkerGraphVertexId) const;

    // Find the children or parents of a vertex of the global marker graph.
    vector<GlobalMarkerGraphVertexId>
        getGlobalMarkerGraphVertexChildren(
        GlobalMarkerGraphVertexId) const;
    vector<GlobalMarkerGraphVertexId>
        getGlobalMarkerGraphVertexParents(
        GlobalMarkerGraphVertexId) const;



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
        GlobalMarkerGraphVertexId,
        GlobalMarkerGraphVertexId
        );

private:
    // Lower-level, more efficient version of the above
    // (but it returns less information).
    void getGlobalMarkerGraphEdgeInfo(
        GlobalMarkerGraphVertexId,
        GlobalMarkerGraphVertexId,
        vector<MarkerInterval>&
        );
public:


    // Create a local marker graph and return its local assembly path.
    // The local marker graph is specified by its start vertex
    // and maximum distance (number of edges) form the start vertex.
    vector<GlobalMarkerGraphVertexId> getLocalAssemblyPath(
        GlobalMarkerGraphVertexId,
        int maxDistance
        );



    // Compute connectivity of the global marker graph.
    // Vertices with more than markerCountOverflow vertices are skipped.
    void createMarkerGraphConnectivity(size_t threadCount);
    void accessMarkerGraphConnectivity(bool accessEdgesReadWrite);
    void checkMarkerGraphConnectivityIsOpen();

    // Find weak edges in the marker graph.
    void flagMarkerGraphWeakEdges(
        size_t lowCoverageThreshold,
        size_t highCoverageThreshold,
        size_t maxDistance);



    // Call this before explore to make the documentation available.
    void setDocsDirectory(const string&);

    // Call this before explore to specify the name of the fasta
    // file containing the reference to be used with Blast commands.
    void setReferenceFastaFileName(const string&);

    // Set up the ConsensusCaller used to compute the "best"
    // base and repeat count at each assembly position.
    // The aregument specifies the consensus caller to be used.
    // For now, the only supported value is "SimpleConsensusCaller".
    void setupConsensusCaller(const string&);


private:

    // Data filled in by the constructor.
    string largeDataFileNamePrefix;
    size_t largeDataPageSize;

    // Function to construct names for binary objects.
    string largeDataName(const string& name) const
    {
        return largeDataFileNamePrefix + name;
    }

    // Various pieces of assembler information stored in shared memory.
    // See class AssemblerInfo for more information.
    MemoryMapped::Object<AssemblerInfo> assemblerInfo;



    /***************************************************************************

    The reads used for this assembly.
    Indexed by ReadId.

    Depending on the setting of assemblerInfo->useRunLengthReads,
    we represent reads in one of two ways:

    - If assemblerInfo->useRunLengthReads is false, we represent
      reads as raw reads just as read from the input fasta files.
      In this case, repeat base counts are not used.

    - If assemblerInfo->useRunLengthReads is true, we use a run-length
      representation (https://en.wikipedia.org/wiki/Run-length_encoding)
      for reads: all repeated bases are removed, and
      for each base we store a repeat base count that says how many
      times that base was repeated in the original read.
      Many assembly phases use only the run-length representation
      (without using the base repeat count).
      This includes the generation of markers, the computation of
      alignments, and the creation of the marker graph.

    For example, suppose we have the following read sequence:

    TAATCATTTTGATGTAAGTCTAAAAATTTCACCTTAATACTTATTTTTCC

    If assemblerInfo->useRunLengthReads is false, the read is represented
    just as written above, this sequence is stored in reads[readId],
    and readRepeatCount is not used.

    If assemblerInfo->useRunLengthReads is true, the read is stored like this:

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
    ReadId readCount() const
    {
        return ReadId(reads.size());
    }
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
        // This should only be called when using run-length read representation.
        CZI_ASSERT(assemblerInfo->useRunLengthReads);

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



    // Alignment candidate found by the MinHash algorithm.
    // They all have readId0<readId1.
    MemoryMapped::Vector<OrientedReadPair> alignmentCandidates;
    void checkAlignmentCandidatesAreOpen() const;



    // Compute a marker alignment of two oriented reads.
    void alignOrientedReads(
        OrientedReadId,
        OrientedReadId,
        size_t maxSkip, // Maximum ordinal skip allowed.
        size_t maxVertexCountPerKmer
    );
    // This lower level version takes as input vectors of
    // markers already sorted by kmerId.
    void alignOrientedReads(
        const vector<MarkerWithOrdinal>& markers0SortedByKmerId,
        const vector<MarkerWithOrdinal>& markers1SortedByKmerId,
        size_t maxSkip,  // Maximum ordinal skip allowed.
        size_t maxVertexCountPerKmer
    );
    // This version allows reusing the AlignmentGraph and Alignment
    void alignOrientedReads(
        const vector<MarkerWithOrdinal>& markers0SortedByKmerId,
        const vector<MarkerWithOrdinal>& markers1SortedByKmerId,
        size_t maxSkip,             // Maximum ordinal skip allowed.
        size_t maxVertexCountPerKmer,
        bool debug,
        AlignmentGraph& graph,
        Alignment& alignment
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
        size_t maxVertexCountPerKmer,
        size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
        size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
    );


    // Given two oriented reads and their computed AlignmentInfo,
    // compute the left and right trim, expressed in markers.
    // This is the minimum number of markers (over the two reads)
    // that are excluded from the alignment on each side.
    // If the trim is too high, the alignment is suspicious.
    pair<uint32_t, uint32_t> computeTrim(
        OrientedReadId orientedReadIds0,
        OrientedReadId orientedReadIds1,
        const AlignmentInfo&) const;

    // The good alignments we found.
    // They are stored with readId0<readId1 and with strand0==0.
    MemoryMapped::Vector<AlignmentData> alignmentData;
    void checkAlignmentDataAreOpen();

    // The alignment table stores the AlignmentData that each oriented read is involved in.
    // Stores, for each OrientedReadId, a vector of indexes into the alignmentData vector.
    // Indexed by OrientedReadId::getValue(),
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> alignmentTable;
    void computeAlignmentTable();

#if 0
    // Thread functions for computeAllAlignments.

    // Compute the alignments and update the disjoint set data structure
    // for each good alignment.
    void computeAllAlignmentsThreadFunction1(size_t threadId);

    // Find the set that each marker belongs to.
    void computeAllAlignmentsThreadFunction2(size_t threadId);

    // Count the number of markers in each vertex.
    void computeAllAlignmentsThreadFunction3(size_t threadId);

    // Use the work area to convert raw vertex ids to final vertex ids.
    void computeAllAlignmentsThreadFunction4(size_t threadId);



    // Data for computeAllAlignments.
    class ComputeAllAlignmentsData {
    public:

        // Parameters.
        size_t maxSkip;
        size_t maxVertexCountPerKmer;
        size_t minAlignedMarkerCount;
        size_t maxTrim;

        // The total number of oriented markers.
        uint64_t orientedMarkerCount;

        // The AlignmentInfo found by each thread.
        vector< vector<AlignmentData> > threadAlignmentData;

        // Disjoint sets data structures.
        MemoryMapped::Vector< std::atomic<DisjointSets::Aint> > disjointSetsData;
        shared_ptr<DisjointSets> disjointSetsPointer;

        // Work area used for multiple purposes.
        // See computeAllAlignments for details.
        MemoryMapped::Vector<GlobalMarkerGraphVertexId> workArea;
    };
    ComputeAllAlignmentsData computeAllAlignmentsData;
#endif


    // Private functions and data used by computeAlignments.
    void computeAlignmentsThreadFunction(size_t threadId);
    class ComputeAlignmentsData {
    public:

        // Parameters.
        size_t maxSkip;
        size_t maxVertexCountPerKmer;
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

    // Given the oriented alignments computed by findOrientedAlignments,
    // compute a vector of alignment coverage for the given oriented read.
    // This is a vector that contains, for each marker of the given
    // oriented read, the number of alignments "covering" that marker.
    // An alignment "covers" as marker if it starts at or before
    // the marker and it ends at or after the marker.
    // It does not matter whether the marker is in the alignment.
    // This is used to detect portions of the read that may contain
    // unreliable alignments due to repeats (most commonly for a human genome,
    // LINE repeats).
    void computeAlignmentCoverage(
        OrientedReadId,
        const vector< pair<OrientedReadId, AlignmentInfo> >&,    // Computed by findOrientedAlignments
        vector<uint32_t>& coverage
        ) const;

    // Find out if an alignment stored in alignmentData is a containment alignment.
    // The alignment is a containment alignment if one (or both)
    // of the two reads are entirely contained in the alignment,
    // except possibly for up to maxTrim markers at each end.
    // The alignment to be checked is specified by its alignmentId, an index
    // into the alignmentData vector.
    bool isContainmentAlignment(
        uint64_t alignmentId,
        size_t maxTrim) const;



    // Read graph and related functions and data.
    // For more information, see comments at the beginning
    // of AssemblerReadGraph.cpp.
public:
    void createReadGraph(uint32_t maxTrim);
private:



    // Edges of the read graph.
    class ReadGraphEdge {
    public:
        // Index in alignmentData vector of the alignment that generated this vertex.
        uint64_t alignmentId : 62;
        // Bidirected edge information.
        // 0 = points towards vertex, 1 = points away from vertex
        // For more information, see comments at the beginning
        // of AssemblerReadGraph.cpp.
        uint64_t direction0: 1;
        uint64_t direction1: 1;
        static const int towardsVertex = 0;
        static const int awayFromVertex = 1;
    };
    MemoryMapped::Vector<ReadGraphEdge> readGraphEdges;

    // Connectivity of the read graph.
    // Stores, for each ReadId, a vector of indexes into the readGraphEdges vector.
    // Indexed by ReadId.
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> readGraphConnectivity;

public:
    // This accesses containingOrientedReadId, readGraphEdges,
    // and readGraphConnectivity.
    void accessReadGraph();
    void checkReadGraphIsOpen();



    // Use the read graph to flag chimeric reads.
    void flagChimericReads(size_t maxDistance, size_t threadCount);
    void accessChimericReadsFlags();
private:
    MemoryMapped::Vector<bool> isChimericRead;
    class FlagChimericReadsData {
    public:
        size_t maxDistance;
    };
    FlagChimericReadsData flagChimericReadsData;
    void flagChimericReadsThreadFunction(size_t threadId);



    // Create a local subgraph of the global read graph,
    // starting at a given vertex and extending out to a specified
    // distance (number of edges).
    // If the specified Readid corresponds to a contained read,
    // which does not have a corresponding vertex in the read graph,
    // the local subgraph starts instead from the containing read
    // of the specified read.
    bool createLocalReadGraph(
        ReadId& readIdStart,    // If the specified read is contained, modified to the containing read.
        uint32_t maxDistance,   // How far to go from starting oriented read.
        uint32_t maxTrim,       // To define alignment containment
        bool allowChimericReads,
        double timeout,         // Or 0 for no timeout.
        LocalReadGraph&);



    // Compute connected components of the read graph.
    // This treats chimeric reads as isolated.
public:
    void computeReadGraphConnectedComponents();



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
        size_t maxVertexCountPerKmer;

        // The total number of oriented markers.
        uint64_t orientedMarkerCount;

        // Disjoint sets data structures.
        MemoryMapped::Vector< std::atomic<DisjointSets::Aint> > disjointSetsData;
        shared_ptr<DisjointSets> disjointSetsPointer;

        // The disjoint set that each oriented marker was assigned to.
        // See createMarkerGraphVertices for details.
        MemoryMapped::Vector<GlobalMarkerGraphVertexId> disjointSetTable;

        // Work area used for multiple purposes.
        // See createMarkerGraphVertices for details.
        MemoryMapped::Vector<GlobalMarkerGraphVertexId> workArea;

        // The markers in each disjoint set with coverage in the requested range.
        MemoryMapped::VectorOfVectors<MarkerId, GlobalMarkerGraphVertexId> disjointSetMarkers;

        // Flag disjoint sets that contain more than one marker on the same oriented read.
        MemoryMapped::Vector<bool> isBadDisjointSet;

    };
    CreateMarkerGraphVerticesData createMarkerGraphVerticesData;



    // The global marker graph vertex corresponding to each marker.
    // Indexed by MarkerId.
    MemoryMapped::Vector<CompressedGlobalMarkerGraphVertexId> globalMarkerGraphVertex;

    // The oriented marker ids of the markers corresponding to
    // each vertex of the global marker graph.
    // Indexed by GlobalMarkerGraphVertexId.
    // For a given vertex, the oriented marker ids are sorted.
    MemoryMapped::VectorOfVectors<MarkerId, CompressedGlobalMarkerGraphVertexId> globalMarkerGraphVertices;
    void checkMarkerGraphVerticesAreAvailable();

    // Check for consistency of globalMarkerGraphVertex and globalMarkerGraphVertices.
    void checkMarkerGraphVertices(
        size_t minCoverage,
        size_t maxCoverage);



    // Marker graph connectivity.
    // Also contains temporary data used by createMarkerGraphConnectivity.
    class MarkerGraphConnectivity {
    public:

        // The edges of the marker graph.
        class Edge {
        public:
            Uint40 source;  // The source vertex (index into globalMarkerGraphVertices).
            Uint40 target;  // The target vertex (index into globalMarkerGraphVertices).
            uint8_t coverage;   // (255 indicates 255 or more).

            // This flag gets set by flagWeakMarkerGraphEdges.
            uint8_t isWeak : 1;

            // Set if this edge was pruned from the strong subgraph.
            uint8_t wasPruned : 1;

            // The remaining flags are currently unused.
            uint8_t flag2 : 1;
            uint8_t flag3 : 1;
            uint8_t flag4 : 1;
            uint8_t flag5 : 1;
            uint8_t flag6 : 1;
            uint8_t flag7 : 1;
            void clearFlags()
            {
                isWeak = 0;
                wasPruned = 0;
                flag2 = 0;
                flag3 = 0;
                flag4 = 0;
                flag5 = 0;
                flag6 = 0;
                flag7 = 0;
            }
            Edge() :
                source(invalidCompressedGlobalMarkerGraphVertexId),
                target(invalidCompressedGlobalMarkerGraphVertexId),
                coverage(0)
            {
                clearFlags();
            }
        };
        MemoryMapped::Vector<Edge> edges;
        const Edge* findEdge(Uint40 source, Uint40 target) const;

        // The MarkerIntervals for each of the above edges.
        MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> edgeMarkerIntervals;

        // The edges and their MarkerIntervals found by each thread.
        // This is temporary and only used inside createMarkerGraphConnectivity.
        vector< shared_ptr<MemoryMapped::Vector<Edge> > > threadEdges;
        vector< shared_ptr< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> > > threadEdgeMarkerIntervals;

        // The edges that each vertex is the source of.
        // Contains indexes into the above edges vector.
        MemoryMapped::VectorOfVectors<Uint40, uint64_t> edgesBySource;

        // The edges that each vertex is the target of.
        // Contains indexes into the above edges vector.
        MemoryMapped::VectorOfVectors<Uint40, uint64_t> edgesByTarget;

    };
    MarkerGraphConnectivity markerGraphConnectivity;
    void createMarkerGraphConnectivityThreadFunction0(size_t threadId);
    void createMarkerGraphConnectivityThreadFunction1(size_t threadId);
    void createMarkerGraphConnectivityThreadFunction2(size_t threadId);
    void createMarkerGraphConnectivityThreadFunction12(size_t threadId, size_t pass);



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
    // invalidGlobalMarkerGraphVertexId.
    GlobalMarkerGraphVertexId getGlobalMarkerGraphVertex(
        OrientedReadId,
        uint32_t ordinal) const;

    // Find the markers contained in a given vertex of the global marker graph.
    // The markers are stored as pairs(oriented read id, ordinal).
    void getGlobalMarkerGraphVertexMarkers(
        GlobalMarkerGraphVertexId,
        vector< pair<OrientedReadId, uint32_t> >&) const;

    // Find the children or parents of a vertex of the global marker graph.
    void getGlobalMarkerGraphVertexChildren(
        GlobalMarkerGraphVertexId,
        vector<GlobalMarkerGraphVertexId>&,
        bool append = false
        ) const;
    void getGlobalMarkerGraphVertexParents(
        GlobalMarkerGraphVertexId,
        vector<GlobalMarkerGraphVertexId>&,
        bool append = false
        ) const;

    // This version also returns the oriented read ids and ordinals
    // that caused a child to be marked as such.
    void getGlobalMarkerGraphVertexChildren(
        GlobalMarkerGraphVertexId,
        vector< pair<GlobalMarkerGraphVertexId, vector<MarkerInterval> > >&,
        vector< pair<GlobalMarkerGraphVertexId, MarkerInterval> >& workArea
        ) const;
    void getGlobalMarkerGraphVertexParents(
        GlobalMarkerGraphVertexId,
        vector< pair<GlobalMarkerGraphVertexId, vector<MarkerInterval> > >&,
        vector< pair<GlobalMarkerGraphVertexId, MarkerInterval> >& workArea
        ) const;

    // Return true if a vertex of the global marker graph has more than
    // one marker for at least one oriented read id.
    bool isBadMarkerGraphVertex(GlobalMarkerGraphVertexId) const;

    // Find out if a vertex is a forward or backward leaf of the pruned
    // strong subgraph of the marker graph.
    // A forward leaf is a vertex with out-degree 0.
    // A backward leaf is a vertex with in-degree 0.
    bool isForwardLeafOfMarkerGraphPrunedStrongSubgraph(GlobalMarkerGraphVertexId) const;
    bool isBackwardLeafOfMarkerGraphPrunedStrongSubgraph(GlobalMarkerGraphVertexId) const;

    // Given an edge of the pruned strong subgraph of the marker graph,
    // return the next/previous edge in the linear chain the edge belongs to.
    // If the edge is the last/first edge in its linear chain, return invalidGlobalMarkerGraphEdgeId.
    GlobalMarkerGraphEdgeId nextEdgeInMarkerGraphPrunedStrongSubgraphChain(GlobalMarkerGraphEdgeId) const;
    GlobalMarkerGraphEdgeId previousEdgeInMarkerGraphPrunedStrongSubgraphChain(GlobalMarkerGraphEdgeId) const;

    // Return the out-degree or in-degree (number of outgoing/incoming edges)
    // of a vertex of the pruned strong subgraph of the marker graph.
    size_t markerGraphPrunedStrongSubgraphOutDegree(GlobalMarkerGraphVertexId) const;
    size_t markerGraphPrunedStrongSubgraphInDegree (GlobalMarkerGraphVertexId) const;

    // Return true if an edge disconnects the local subgraph.
    bool markerGraphEdgeDisconnectsLocalStrongSubgraph(
        GlobalMarkerGraphEdgeId edgeId,
        size_t maxDistance,

        // Work areas, to reduce memory allocation activity.

        // Each of these two must be sized maxDistance+1.
        array<vector< vector<GlobalMarkerGraphEdgeId> >, 2>& verticesByDistance,

        // Each of these two must be sized globalMarkerGraphVertices.size()
        // and set to all false on entry.
        // It is left set to all false on exit, so it can be reused.
        array<vector<bool>, 2>& vertexFlags
        ) const;



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
        GlobalMarkerGraphVertexId,
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
        LocalMarkerGraph&
        );
    bool extractLocalMarkerGraphUsingStoredConnectivity(
        GlobalMarkerGraphVertexId,
        int distance,
        double timeout,                 // Or 0 for no timeout.
        bool useWeakEdges,
        bool usePrunedEdges,
        LocalMarkerGraph&
        );

    // Compute consensus sequence for a vertex of the marker graph.
    void computeMarkerGraphVertexConsensusSequence(
        GlobalMarkerGraphVertexId,
        vector<Base>& sequence,
        vector<uint32_t>& repeatCounts
        );

    // Compute consensus sequence for an edge of the marker graph.
    // This includes the k bases corresponding to the flanking markers,
    // but computed only using reads on this edge.
    void computeMarkerGraphEdgeConsensusSequenceUsingSeqan(
        GlobalMarkerGraphEdgeId,
        vector<Base>& sequence,
        vector<uint32_t>& repeatCounts
        );
    void computeMarkerGraphEdgeConsensusSequenceUsingSpoa(
        GlobalMarkerGraphEdgeId,
        vector<Base>& sequence,
        vector<uint32_t>& repeatCounts
        );



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
    void accessAssemblyGraphEdges();
private:

    // Extract a local assembly graph from the global assembly graph.
    // This returns false if the timeout was exceeded.
    bool extractLocalAssemblyGraph(
        AssemblyGraph::VertexId,
        int distance,
        double timeout,
        LocalAssemblyGraph&) const;

    // Assemble sequence for a vertex of the assembly graph.
    // Optionally outputs detailed assembly information
    // in html (skipped if the html pointer is 0).
    void assembleAssemblyGraphVertex(
        AssemblyGraph::VertexId,
        vector<Base>&,
        vector<uint32_t>& repeatCounts,
        ostream* html = 0);


    // Assemble sequence for all vertices of the assembly graph.
public:
    void assemble(size_t threadCount);
    void accessAssemblyGraphSequences();
    void computeAssemblyStatistics();
private:
    class AssembleData {
    public:

        // The results created by each thread.
        // All indexed by threadId.
        vector< vector<AssemblyGraph::VertexId> > vertices;
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



    // Data and functions used for the http server.
    // This function puts the server into an endless loop
    // of processing requests.
    void fillServerFunctionTable();
    void processRequest(
        const vector<string>& request,
        ostream&,
        const BrowserInformation&) override;
    void writeHtmlBegin(ostream&) const;
    void writeHtmlEnd(ostream&) const;
    void writeMakeAllTablesSelectable(ostream&) const;
    void writeNavigation(ostream&) const;
    void writeNavigation(
        ostream& html,
        const string& title,
        const vector<pair <string, string> >&) const;
    void exploreSummary(const vector<string>&, ostream&);
    void exploreRead(const vector<string>&, ostream&);
    void blastRead(const vector<string>&, ostream&);
    void exploreAlignments(const vector<string>&, ostream&);
    void exploreAlignment(const vector<string>&, ostream&);
    void computeAllAlignments(const vector<string>&, ostream&);
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
        GlobalMarkerGraphVertexId vertexId;
        bool vertexIdIsPresent;
        uint32_t maxDistance;
        bool maxDistanceIsPresent;
        bool detailed;
        bool useStoredConnectivity;
        bool useWeakEdges;
        bool usePrunedEdges;
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
        void writeForm(ostream&, GlobalMarkerGraphVertexId vertexCount) const;
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



    // Access all available assembly data, without thorwing an exception
    // on failures.
public:
    void accessAllSoft();



    // Functions and data used by the http server
    // for display of the local assembly graph.
private:
    void exploreAssemblyGraph(const vector<string>&, ostream&);
    class LocalAssemblyGraphRequestParameters {
    public:
        AssemblyGraph::VertexId vertexId;
        bool vertexIdIsPresent;
        uint32_t maxDistance;
        bool maxDistanceIsPresent;
        bool detailed;
        uint32_t sizePixels;
        bool sizePixelsIsPresent;
        double timeout;
        bool timeoutIsPresent;
        void writeForm(ostream&, AssemblyGraph::VertexId vertexCount) const;
        bool hasMissingRequiredParameters() const;
    };
    void getLocalAssemblyGraphRequestParameters(
        const vector<string>&,
        LocalAssemblyGraphRequestParameters&) const;
    void exploreAssemblyGraphVertex(const vector<string>&, ostream&);


    // The ConsensusCaller used to compute the "best"
    // base and repeat count at each assembly position.
    shared_ptr<ConsensusCaller> consensusCaller;
};

#endif
