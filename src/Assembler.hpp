#ifndef CZI_SHASTA_ASSEMBLER_HPP
#define CZI_SHASTA_ASSEMBLER_HPP

// shasta
#include "Alignment.hpp"
#include "dset64.hpp"
#include "HttpServer.hpp"
#include "Kmer.hpp"
#include "LongBaseSequence.hpp"
#include "Marker.hpp"
#include "MarkerId.hpp"
#include "MemoryMappedObject.hpp"
#include "MultitreadedObject.hpp"
#include "Overlap.hpp"
#include "ReadId.hpp"

// Standard library.
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
        class LocalReadGraph;
        class LocalMarkerGraph2;
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

    The constructors specify the file name prefixes for binary data files.
    There are two prefixes, one used for small data and one used for large data.
    If these are directory names, they must include the final "/".

    The constructors also specify the page size for small and large binary data files.
    Typically, small binary data files will reside in a regular
    directory on disk or on /dev/shm mapped backed by 4K pages,
    while large binary data file will reside in a huge page
    file system backed by 2MB pages.
    1GB huge pages are also supported.
    The page sizes specified here must be equal to, or be an exact multiple of,
    the actual size of the pages backing the data.
    If the system has no large pages and it is not possible to change that,
    use 4096 for both page sizes, and for performance place both
    the small and large binary data under /dev/shm (in-memory filesystem).


    ***************************************************************************/

    // Constructor to be called one to create a new run.
    Assembler(
        const string& smallDataFileNamePrefix,
        const string& largeDataFileNamePrefix,
        size_t smallDataPageSize,
        size_t largeDataPageSize,
        bool useRunLengthReads
        );

    // Constructor to be called to continue an existing run.
    Assembler(
        const string& smallDataFileNamePrefix,
        const string& largeDataFileNamePrefix,
        size_t smallDataPageSize,
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

    // Use the minHash algorithm to find pairs of overlapping oriented reads.
    // Use as features sequences of m consecutive special k-mers.
    void findOverlaps(
        size_t m,                       // Number of consecutive k-mers that define a feature.
        size_t minHashIterationCount,   // Number of minHash iterations.
        size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for minHash.
        size_t maxBucketSize,           // The maximum size for a bucket to be used.
        size_t minFrequency,            // Minimum number of minHash hits for a pair to become a candidate.
        size_t threadCount
    );
    void accessOverlaps();

    // Write the reads that overlap a given read.
    void writeOverlappingReads(ReadId, Strand, const string& fileName);

    // Create a local read graph starting from a given oriented read.
    void createLocalReadGraph(
        ReadId, Strand,
        size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
        size_t maxTrim,                 // Maximum left/right trim (in bases) to generate an edge.
        uint32_t distance               // How far to go from starting oriented read.
    );

    // Compute connected components of the global read graph.
    void computeReadGraphComponents(
        size_t minComponentSize,        // Minimum size for a connected component to be kept.
        size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
        size_t maxTrim                  // Maximum left/right trim to generate an edge
        );

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
        size_t maxTrim                  // Maximum trim allowed in an alignment.
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



    // Compute an Alignment for each Overlap, but  only store the AlignmentInfo.
    // Optionally, the alignments are used for creation of the global marker graph.
    void computeAllAlignments(

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

        // Minimum coverage (number of markers) for a vertex
        // of the marker graph to be kept.
        size_t minCoverage,

        // Number of threads. If zero, a number of threads equal to
        // the number of virtual processors is used.
        size_t threadCount
    );
    void accessAlignmentData();



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



    // Compute connectivity of the global marker graph.
    // This code is currently not in use.
    // Vertices with more than markerCountOverflow are skipped.
    void createMarkerGraphConnectivity(
        size_t threadCount,
        size_t markerCountOverflow
        );
    void accessMarkerGraphConnectivity(bool accessEdgesReadWrite);

    // Flag as not good a marker graph edge if:
    // - It has coverage<minCoverage, AND
    // - A path of length <= maxPathLength edges exists that:
    //    * Starts at the source vertex of the edge.
    //    * Ends at the target vertex of the edge.
    //    * Only uses edges with coverage>=minCoverage.
    void flagMarkerGraphEdges(
        size_t threadCount,
        size_t minCoverage,
        size_t maxPathLength);



    // Call this before explore to make the documentation available.
    void setDocsDirectory(const string&);

private:

    // Data filled in by the constructor.
    string smallDataFileNamePrefix;
    string largeDataFileNamePrefix;
    size_t smallDataPageSize;
    size_t largeDataPageSize;

    // Functions to construct names for small and large binary objects.
    string smallDataName(const string& name) const
    {
        return smallDataFileNamePrefix + name;
    }
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
    Base getOrientedReadBase(OrientedReadId orientedReadId, uint32_t position)
    {
        const auto& read = reads[orientedReadId.getReadId()];
        if(orientedReadId.getStrand() == 0) {
            return read[position];
        } else {
            return read[read.baseCount-1-position].complement();
        }
    }

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



    // Pairs of overlapping oriented reads.
    // This is a global vector that stores all the overlaps.
    // The overlap table defined below can be used to locate
    // all the overlaps that an oriented read is involved in.
    MemoryMapped::Vector<Overlap> overlaps;
    void checkOverlapsAreOpen() const;



    // The overlaps define a global overlap graph,
    // an undirected graph in which each vertex corresponds to
    // an oriented read. Two vertices are joined by an edge
    // if there is an overlap between the corresponding
    // oriented reads. We want to process each connected
    // component separately, so we compute connected components
    // of the overlap graph.

    // The connected component that each oriented read belongs to,
    // or std::numeric_limits<ReadId>::max() if the oriented read
    // belongs to a connected component that was discarded
    // because it was too small.
    // Indexed by OrientedReadId::getValue().
    MemoryMapped::Vector<ReadId> overlapGraphComponent;

    // The OrientedReadId's of each connected component.
    // Sorted by decreasing component size.
    MemoryMapped::VectorOfVectors<OrientedReadId, ReadId> overlapGraphComponents;

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

    // Create a local read graph starting from a given oriented read
    // and walking out a given distance on the global read graph.
    void createLocalReadGraph(
        OrientedReadId,
        size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
        size_t maxTrim,                 // Maximum left/right trim to generate an edge.
        uint32_t distance               // How far to go from starting oriented read.
    );
    bool createLocalReadGraph(
        OrientedReadId,
        size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
        size_t maxTrim,                 // Maximum left/right trim to generate an edge.
        uint32_t distance,              // How far to go from starting oriented read.
        double timeout,                 // Or 0 for no timeout.
        LocalReadGraph&
    );

    // Write in fasta format the sequences of the vertices of a local read graph.
    void writeLocalReadGraphToFasta(
        const LocalReadGraph&,
        const string& fileName);

    // Compute marker alignments of an oriented read with all reads
    // for which we have an Overlap.
    void alignOverlappingOrientedReads(
        OrientedReadId,
        size_t maxSkip,                 // Maximum ordinal skip allowed.
        size_t maxVertexCountPerKmer,
        size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
        size_t maxTrim                  // Maximum trim allowed in an alignment.
    );


    // Given two oriented reads and their computed AlignmentInfo,
    // compute the left and right trim.
    // This is the minimum number of bases (over the two reads)
    // that are excluded from the alignment on each size.
    pair<uint32_t, uint32_t> computeTrim(
        OrientedReadId orientedReadIds0,
        OrientedReadId orientedReadIds1,
        const AlignmentInfo&);

    // The good alignments we found.
    // These are computed with the first oriented read on strand 0.
    MemoryMapped::Vector<AlignmentData> alignmentData;
    void checkAlignmentDataAreOpen();

    // The alignment table stores the AlignmentData that each oriented read is involved in.
    // Stores, for each OrientedReadId, a vector of indexes into the alignmentData vector.
    // Indexed by OrientedReadId::getValue(),
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> alignmentTable;
    void computeAlignmentTable();


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
        std::shared_ptr<DisjointSets> disjointSetsPointer;

        // Work area used for multiple purposes.
        // See computeAllAlignments for details.
        MemoryMapped::Vector<GlobalMarkerGraphVertexId> workArea;
    };
    ComputeAllAlignmentsData computeAllAlignmentsData;



    // The global marker graph vertex corresponding to each marker.
    // Indexed by MarkerId.
    MemoryMapped::Vector<CompressedGlobalMarkerGraphVertexId> globalMarkerGraphVertex;

    // The oriented marker ids of the markers corresponding to
    // each vertex of the global marker graph.
    // Indexed by GlobalMarkerGraphVertexId.
    MemoryMapped::VectorOfVectors<MarkerId, CompressedGlobalMarkerGraphVertexId> globalMarkerGraphVertices;
    void checkMarkerGraphVerticesAreAvailable();



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
            bool isGood = true;
        };
        MemoryMapped::Vector<Edge> edges;
        const Edge* findEdge(Uint40 source, Uint40 target) const;

        // The edges found by each thread.
        // This is temporary and only used inside createMarkerGraphConnectivity.
        vector< std::shared_ptr<MemoryMapped::Vector<Edge> > > threadEdges;

        // The edges that each vertex is the source of.
        // Contains indexes into the above edges vector.
        MemoryMapped::VectorOfVectors<Uint40, uint64_t> edgesBySource;

        // The edges that each vertex is the target of.
        // Contains indexes into the above edges vector.
        MemoryMapped::VectorOfVectors<Uint40, uint64_t> edgesByTarget;

        // Vertices with more than markerCountOverflow are skipped.
        size_t markerCountOverflow;

    };
    MarkerGraphConnectivity markerGraphConnectivity;
    void createMarkerGraphConnectivityThreadFunction0(size_t threadId);
    void createMarkerGraphConnectivityThreadFunction1(size_t threadId);
    void createMarkerGraphConnectivityThreadFunction2(size_t threadId);
    void createMarkerGraphConnectivityThreadFunction12(size_t threadId, size_t pass);



    // Data used by flagMarkerGraphEdges.
    class FlagMarkerGraphEdgesData {
    public:
        size_t minCoverage;
        size_t maxPathLength;
    };
    FlagMarkerGraphEdgesData flagMarkerGraphEdgesData;
    void flagMarkerGraphEdgesThreadFunction(size_t threadId);




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
    class MarkerGraphNeighborInfo {
    public:
        OrientedReadId orientedReadId;
        uint32_t ordinal0;
        uint32_t ordinal1;
        bool operator<(const MarkerGraphNeighborInfo& that) const
        {
            return tie(orientedReadId, ordinal0, ordinal1) <
                tie(that.orientedReadId, that.ordinal0, that.ordinal1);
        }
    };
    void getGlobalMarkerGraphVertexChildren(
        GlobalMarkerGraphVertexId,
        vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > >&,
        vector< pair<GlobalMarkerGraphVertexId, MarkerGraphNeighborInfo> >& workArea
        ) const;
    void getGlobalMarkerGraphVertexParents(
        GlobalMarkerGraphVertexId,
        vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > >&,
        vector< pair<GlobalMarkerGraphVertexId, MarkerGraphNeighborInfo> >& workArea
        ) const;

    // Return true if a vertex of the global marker graph has more than
    // one marker for at least one oriented read id.
    bool isBadMarkerGraphVertex(GlobalMarkerGraphVertexId) const;



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
        LocalMarkerGraph2&
        );



    // Data and functions used for the http server.
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
    void exploreOverlappingReads(const vector<string>&, ostream&);
    void exploreAlignment(const vector<string>&, ostream&);
    void exploreReadGraph(const vector<string>&, ostream&);
    class HttpServerData {
    public:

        using ServerFunction = void (Assembler::*) (
            const vector<string>& request,
            ostream&);
        std::map<string, ServerFunction> functionTable;
        string docsDirectory;

    };
    HttpServerData httpServerData;

    // Functions and data used for display of the local marker graph.
    void exploreMarkerGraph(const vector<string>&, ostream&);
    class LocalMarkerGraphRequestParameters {
    public:
        ReadId readId;
        bool readIdIsPresent;
        Strand strand;
        bool strandIsPresent;
        uint32_t ordinal;
        bool ordinalIsPresent;
        uint32_t maxDistance;
        bool maxDistanceIsPresent;
        bool detailed;
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
        void writeForm(ostream&, size_t readCount) const;
        bool hasMissingRequiredParameters() const;
    };
    void getLocalMarkerGraphRequestParameters(
        const vector<string>&,
        LocalMarkerGraphRequestParameters&) const;
    void showLocalMarkerGraphAlignments(
        ostream& html,
        const LocalMarkerGraph2&,
        const LocalMarkerGraphRequestParameters&
        );
};

#endif
