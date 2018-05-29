#ifndef CZI_NANOPORE2_ASSEMBLER_HPP
#define CZI_NANOPORE2_ASSEMBLER_HPP

// Nanopore2
#include "Alignment.hpp"
#include "Kmer.hpp"
#include "LongBaseSequence.hpp"
#include "Marker.hpp"
#include "MemoryMappedObject.hpp"
#include "MultitreadedObject.hpp"
#include "Overlap.hpp"
#include "ReadId.hpp"

// Standard library.
#include "string.hpp"

namespace ChanZuckerberg {
    namespace Nanopore2 {

        // Forward declarations of classes defined here.
        class Assembler;
        class AssemblerInfo;

        // Forward declarations of classes defined elsewhere.
        class Alignment;
        class AlignmentGraph;
        class AlignmentInfo;
        class LocalReadGraph;
    }
}



// Class used to store various pieces of assembler information in shared memory.
class ChanZuckerberg::Nanopore2::AssemblerInfo {
public:

    // The length of k-mers used to define markers.
    size_t k;
};



class ChanZuckerberg::Nanopore2::Assembler :
    public MultithreadedObject<Assembler> {
public:



    // The constructor specifies the file name prefixes for binary data files.
    // There are two prefixes, one used for small data and one used for large data.
    // If these are directory names, they must include the final "/".
    // It also specifies the page size for small and large binary data files.
    // Typically, small binary data files will reside in a regular
    // directory on disk or on /dev/shm mapped backed by 4K pages,
    // while large binary data file will reside in a huge page
    // file system backed by 2MB pages.
    // 1GB huge pages are also supported.
    // The page sizes specified here must be equal to, or be an exact multiple of,
    // the actual size of the pages backing the data.
    // If the system has no large pages and it is not possible to change that,
    // use 4096 for both page sizes, and for performance place both
    // the small and large binary data under /dev/shm (in-memory filesystem).
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
        size_t minFrequency,            // Minimum number of minHash hits to generate an edge.
        size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
        size_t maxTrim,                 // Maximum left/right trim to generate an edge.
        size_t distance                 // How far to go from starting oriented read.
    );

    // Compute connected components of the global read graph.
    void computeReadGraphComponents(
        size_t minFrequency,            // Minimum number of minHash hits to generate an edge.
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

    // Compute a local marker graph for a set of oriented reads.
    void createLocalMarkerGraph(
        const vector< pair<ReadId, Strand> >&,
        bool alignAllPairs,
        size_t alignmentMaxSkip,
        size_t alignmentMaxVertexCountPerKmer,
        size_t minAlignedMarkerCount,
        size_t minCoverage,
        size_t minConsensus);

    // Create the local marker graph that corresponds to a local read graph
    // constructed starting at a given oriented read and extending out
    // up to a specified distance.
    void createLocalMarkerGraph(
        ReadId, Strand,
        size_t minFrequency,            // Minimum number of minHash hits to generate an edge.
        size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
        size_t maxTrim,                 // Maximum left/right trim to generate an edge.
        size_t distance,                // How far to go from starting oriented read.
        size_t alignmentMaxSkip,
        size_t alignmentMaxVertexCountPerKmer,
        size_t minCoverage,
        size_t minConsensus);

    // Compute an Alignment for each Overlap.
    // Only store the AlignmentInfo.
    void computeAllAlignments(
        size_t maxSkip, // Maximum ordinal skip allowed.
        size_t maxVertexCountPerKmer,
        size_t threadCount
    );
    void accessAlignmentInfos();

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

    // The reads used for this assembly.
    // Indexed by ReadId.
    LongBaseSequences reads;
    ReadId readCount() const
    {
        return ReadId(reads.size());
    }
    void checkReadsAreOpen() const;
    void checkReadNamesAreOpen() const;
    void checkReadId(ReadId) const;

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


    // The markers on all reads. Indexed by ReadId.
    MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t> markers;
    void getMarkers(ReadId, vector<Marker>&) const;
    void getMarkers(OrientedReadId, vector<Marker>&);
    void checkMarkersAreOpen() const;



    // Pairs of overlapping oriented reads.
    // This is a global vector that stores all the overlaps.
    // The overlap table defined below can be used to locate
    // all the overlaps that an oriented read is involved in.
    MemoryMapped::Vector<Overlap> overlaps;

    // The overlap table stores the overlaps that each oriented read is involved in.
    // Stores, for each OrientedReadId, a vector of indexes into the overlaps vector.
    // Indexed by OrientedReadId::getValue(),
    MemoryMapped::VectorOfVectors<uint64_t, uint64_t> overlapTable;
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
        const vector<Marker>& markers0SortedByKmerId,
        const vector<Marker>& markers1SortedByKmerId,
        size_t maxSkip,  // Maximum ordinal skip allowed.
        size_t maxVertexCountPerKmer
    );
    // This version allows reusing the AlignmentGraph and Alignment
    void alignOrientedReads(
        const vector<Marker>& markers0SortedByKmerId,
        const vector<Marker>& markers1SortedByKmerId,
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
        size_t minFrequency,            // Minimum number of minHash hits to generate an edge.
        size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
        size_t maxTrim,                 // Maximum left/right trim to generate an edge.
        size_t distance                 // How far to go from starting oriented read.
    );
    void createLocalReadGraph(
        OrientedReadId,
        size_t minFrequency,            // Minimum number of minHash hits to generate an edge.
        size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
        size_t maxTrim,                 // Maximum left/right trim to generate an edge.
        size_t distance,                // How far to go from starting oriented read.
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
    // This takes as input markers sorted by position.
    pair<uint32_t, uint32_t> computeTrim(
        OrientedReadId orientedReadIds0,
        OrientedReadId orientedReadIds1,
        vector<Marker>& markers0,
        vector<Marker>& markers1,
        const AlignmentInfo&);

    // The AlignmentInfo corresponding to each overlap.
    MemoryMapped::Vector<AlignmentInfo> alignmentInfos;
    void checkAlignmentInfosAreOpen();
    void computeAllAlignmentsThreadFunction(size_t threadId);
    class ComputeAllAlignmentsData {
    public:
        size_t maxSkip;
        size_t maxVertexCountPerKmer;
    };
    ComputeAllAlignmentsData computeAllAlignmentsData;


    // Compute a local marker graph for a set of oriented reads.
    void createLocalMarkerGraph(
        const vector<OrientedReadId>&,
        bool alignAllPairs,
        size_t alignmentMaxSkip,
        size_t alignmentMaxVertexCountPerKmer,
        size_t minAlignedMarkerCount,
        size_t minCoverage,
        size_t minConsensus);

};

#endif
