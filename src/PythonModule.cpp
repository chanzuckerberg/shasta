#ifdef SHASTA_PYTHON_API

// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "Base.hpp"
#include "CompactUndirectedGraph.hpp"
#include "compressAlignment.hpp"
#include "deduplicate.hpp"
#include "dset64Test.hpp"
#include "shastaLapack.hpp"
#include "LongBaseSequence.hpp"
#include "mappedCopy.hpp"
#include "Map.hpp"
#include "MedianConsensusCaller.hpp"
#include "MemoryMappedAllocator.hpp"
#include "MultithreadedObject.hpp"
#include "ShortBaseSequence.hpp"
#include "splitRange.hpp"
#include "testSpoa.hpp"
#include "shortestPathBoundedDistance.hpp"
#include "SimpleBayesianConsensusCaller.hpp"
using namespace shasta;

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
using namespace pybind11;



PYBIND11_MODULE(shasta, module)
{

    // Expose class OrientedReadPair to Python.
    class_<OrientedReadPair>(module, "OrientedReadPair")
        .def_readonly("readIds", &OrientedReadPair::readIds)
        .def_readonly("isSameStrand", &OrientedReadPair::isSameStrand)
        ;

    // Expose class Reads to Python
    class_<Reads>(module, "Reads")
        .def("readCount", &Reads::readCount, "Get the number of reads.")
        .def("writeReads",
            &Reads::writeReads,
            "Write all reads to a file in fasta format.",
            arg("fileName"))
        .def("writeRead",
            (
                void (Reads::*)
                (ReadId, const string&)
            )
            &Reads::writeRead,
            "Write one read to a file in fasta format.",
            arg("readId"),
            arg("fileName"))
        .def("writeOrientedRead",
            (
                void (Reads::*)
                (ReadId, Strand, const string&)
            )
            &Reads::writeOrientedRead,
            "Write one oriented read to a file in fasta format.",
            arg("readId"),
            arg("strand"),
            arg("fileName"))
        .def("getReadId",
            (
                ReadId (Reads::*)
                (const string&) const
            )
            &Reads::getReadId,
            "Find the ReadId corresponding to a given read name.")
        ;



    // Expose class AlignOptions to Python.
    class_<AlignOptions>(module, "AlignOptions")
        .def(pybind11::init<>())
        .def_readwrite("alignMethod", &AlignOptions::alignMethod)
        .def_readwrite("maxSkip", &AlignOptions::maxSkip)
        .def_readwrite("maxDrift", &AlignOptions::maxDrift)
        .def_readwrite("maxTrim", &AlignOptions::maxTrim)
        .def_readwrite("maxMarkerFrequency", &AlignOptions::maxMarkerFrequency)
        .def_readwrite("minAlignedMarkerCount", &AlignOptions::minAlignedMarkerCount)
        .def_readwrite("minAlignedFraction", &AlignOptions::minAlignedFraction)
        .def_readwrite("matchScore", &AlignOptions::matchScore)
        .def_readwrite("mismatchScore", &AlignOptions::mismatchScore)
        .def_readwrite("gapScore", &AlignOptions::gapScore)
        .def_readwrite("downsamplingFactor", &AlignOptions::downsamplingFactor)
        .def_readwrite("bandExtend", &AlignOptions::bandExtend)
        .def_readwrite("maxBand", &AlignOptions::maxBand)
        .def_readwrite("sameChannelReadAlignmentSuppressDeltaThreshold",
            &AlignOptions::sameChannelReadAlignmentSuppressDeltaThreshold)
        .def_readwrite("suppressContainments", &AlignOptions::suppressContainments)
        .def_readwrite("align4DeltaX", &AlignOptions::align4DeltaX)
        .def_readwrite("align4DeltaY", &AlignOptions::align4DeltaY)
        .def_readwrite("align4MinEntryCountPerCell", &AlignOptions::align4MinEntryCountPerCell)
        .def_readwrite("align4MaxDistanceFromBoundary", &AlignOptions::align4MaxDistanceFromBoundary)
        ;



    // Expose class Assembler to Python.
    class_<Assembler>(module, "Assembler")

        // Constructor.
        .def(pybind11::init<const string&, bool, size_t>(),
            "Assembler constructor.",
            arg("largeDataFileNamePrefix") = "Data/",
            arg("createNew") = false,
            arg("largeDataPageSize") = 2*1024*1024)



        // Reads
        .def("getReads", &Assembler::getReads, return_value_policy::reference)
        .def("histogramReadLength",
            &Assembler::histogramReadLength,
            "Create a histogram of read length and write it to a csv file.",
            arg("fileName") = "ReadLengthHistogram.csv")

        // K-mers.
        .def("accessKmers",
            &Assembler::accessKmers)
        .def("writeKmers",
            &Assembler::writeKmers,
            arg("fileName") = "Kmers.csv")
        .def("randomlySelectKmers",
            &Assembler::randomlySelectKmers,
            arg("k"),
            arg("probability"),
            arg("seed") = 231)
        .def("selectKmersBasedOnFrequency",
            &Assembler::selectKmersBasedOnFrequency,
            arg("k"),
            arg("markerDensity"),
            arg("seed") = 231,
            arg("enrichmentThreshold"),
            arg("threadCount") = 0)
        .def("selectKmers2",
            &Assembler::selectKmers2,
            arg("k"),
            arg("markerDensity"),
            arg("seed") = 231,
            arg("enrichmentThreshold"),
            arg("threadCount") = 0)
        .def("selectKmers4",
            &Assembler::selectKmers4,
            arg("k"),
            arg("markerDensity"),
            arg("seed") = 231,
            arg("distanceThreshold"),
            arg("threadCount") = 0)



         // Markers.
        .def("accessMarkers",
            &Assembler::accessMarkers)
        .def("findMarkers",
            &Assembler::findMarkers,
            "Find markers in reads.",
            arg("threadCount") = 0)
        .def("writeMarkers",
            (
                void (Assembler::*)
                (ReadId, Strand, const string&)
            )
            &Assembler::writeMarkers,
            "Write the markers of an oriented read.",
            arg("readId"),
            arg("strand"),
            arg("fileName"))
        .def("getMarkers",
            &Assembler::getMarkers)
        .def("writeMarkerFrequency",
            &Assembler::writeMarkerFrequency)
        .def("computeSortedMarkers",
            &Assembler::computeSortedMarkers,
            arg("threadCount") = 0)
        .def("accessSortedMarkers",
            &Assembler::accessSortedMarkers)


        // Alignment candidates.
        .def("findAlignmentCandidatesLowHash0",
            &Assembler::findAlignmentCandidatesLowHash0,
            arg("m"),
            arg("hashFraction"),
            arg("minHashIterationCount"),
            arg("alignmentCandidatesPerRead"),
            arg("log2MinHashBucketCount") = 0,
            arg("minBucketSize"),
            arg("maxBucketSize"),
            arg("minFrequency"),
            arg("threadCount") = 0)
        .def("findAlignmentCandidatesLowHash1",
            &Assembler::findAlignmentCandidatesLowHash1,
            arg("m"),
            arg("hashFraction"),
            arg("minHashIterationCount"),
            arg("log2MinHashBucketCount") = 0,
            arg("minBucketSize"),
            arg("maxBucketSize"),
            arg("minFrequency"),
            arg("threadCount") = 0)
        .def("accessAlignmentCandidates",
            &Assembler::accessAlignmentCandidates)
        .def("accessAlignmentCandidateTable",
            &Assembler::accessAlignmentCandidateTable)
        .def("getAlignmentCandidates",
            &Assembler::getAlignmentCandidates)
        .def("writeOverlappingReads",
            &Assembler::writeOverlappingReads,
            "Write in fasta format the reads that overlap a given read.",
            arg("readId"),
            arg("strand"),
            arg("fileName") = "OverlappingReads.fasta")
        .def("flagPalindromicReads",
            &Assembler::flagPalindromicReads,
            arg("maxSkip"),
            arg("maxDrift"),
            arg("maxMarkerFrequency"),
            arg("alignedFractionThreshold"),
            arg("nearDiagonalFractionThreshold"),
            arg("deltaThreshold"),
            arg("writeTocsv") = false,
            arg("threadCount") = 0)

        // Alignments.
        .def("writeAlignmentCandidates",
             &Assembler::writeAlignmentCandidates,
             arg("useReadName") = false,
             arg("verbose") = false)
        .def("writeAlignmentDetails",
             &Assembler::writeAlignmentDetails)
        .def("writeLocalAlignmentCandidateReads",
             &Assembler::writeLocalAlignmentCandidateReads,
             arg("readId"),
             arg("strand"),
             arg("maxDistance"),
             arg("useReadName"),
             arg("allowChimericReads"),
             arg("allowCrossStrandEdges"),
             arg("allowInconsistentAlignmentEdges"))
        .def("alignOrientedReads",
            (
                void (Assembler::*)
                (ReadId, Strand, ReadId, Strand, size_t, size_t, uint32_t)
            )
            &Assembler::alignOrientedReads,
            arg("readId0"),
            arg("strand0"),
            arg("readId1"),
            arg("strand1"),
            arg("maxSkip"),
            arg("maxDrift"),
            arg("maxMarkerFrequency"))
        .def("alignOrientedReads1",
            (
                void (Assembler::*)
                (ReadId, Strand, ReadId, Strand, int, int, int)
            )
            &Assembler::alignOrientedReads1,
            arg("readId0"),
            arg("strand0"),
            arg("readId1"),
            arg("strand1"),
            arg("matchCount"),
            arg("mismatchCount"),
            arg("gapCount"))
        .def("alignOrientedReads4",
            (
                void (Assembler::*)
                (ReadId, Strand, ReadId, Strand,
                    uint64_t, uint64_t, uint64_t, uint64_t,
                    uint64_t, double, uint64_t, uint64_t, uint64_t, uint64_t,
                    int64_t, int64_t, int64_t) const
            )
            &Assembler::alignOrientedReads4,
            arg("readId0"),
            arg("strand0"),
            arg("readId1"),
            arg("strand1"),
            arg("deltaX"),
            arg("deltaY"),
            arg("minEntryCountPerCell"),
            arg("maxDistanceFromBoundary"),
            arg("minAlignedMarkerCount"),
            arg("minAlignedFraction"),
            arg("maxSkip"),
            arg("maxDrift"),
            arg("maxTrim"),
            arg("maxBand"),
            arg("matchScore"),
            arg("mismatchScore"),
            arg("gapScore"))
        .def("alignOverlappingOrientedReads",
            (
                void (Assembler::*)
                (ReadId, Strand, size_t, size_t, uint32_t, size_t, size_t)
            )
            &Assembler::alignOverlappingOrientedReads,
            arg("readId"),
            arg("strand"),
            arg("maxSkip"),
            arg("maxDrift"),
            arg("maxMarkerFrequency"),
            arg("minAlignedMarkerCount"),
            arg("maxTrim")
            )

        // Compute an alignment for each alignment candidate.
        .def("computeAlignments",
            &Assembler::computeAlignments)
        .def("accessCompressedAlignments",
            &Assembler::accessCompressedAlignments)
        .def("accessAlignmentData",
            &Assembler::accessAlignmentData)
        .def("accessAlignmentDataReadWrite",
            &Assembler::accessAlignmentDataReadWrite)
        .def("analyzeAlignmentMatrix",
            &Assembler::analyzeAlignmentMatrix,
            arg("readId0"),
            arg("strand0"),
            arg("readId1"),
            arg("strand1"))



        // Undirected read graph
        .def("createReadGraph",
            &Assembler::createReadGraph,
            arg("maxAlignmentCount"),
            arg("maxTrim"))
        .def("createReadGraph2",
             &Assembler::createReadGraph2,
            arg("maxAlignmentCount"),
            arg("markerCountPercentile"),
            arg("alignedFractionPercentile"),
            arg("maxSkipPercentile"),
            arg("maxDriftPercentile"),
            arg("maxTrimPercentile"))
        .def("createReadGraphUsingPseudoPaths",
             &Assembler::createReadGraphUsingPseudoPaths,
             arg("matchScore"),
             arg("mismatchScore"),
             arg("gapScore"),
             arg("mismatchSquareFactor"),
             arg("minScore"),
             arg("maxAlignmentCount"),
             arg("threadCount") = 1
             )
        .def("accessReadGraph",
            &Assembler::accessReadGraph)
        .def("accessReadGraphReadWrite",
            &Assembler::accessReadGraphReadWrite)
        .def("flagCrossStrandReadGraphEdges",
            &Assembler::flagCrossStrandReadGraphEdges,
            arg("maxDistance"),
            arg("threadCount") = 0)
        .def("flagChimericReads",
             &Assembler::flagChimericReads,
            arg("maxChimericReadDistance"),
            arg("threadCount") = 0)
        .def("computeReadGraphConnectedComponents",
            &Assembler::computeReadGraphConnectedComponents,
            arg("minComponentSize"))
        .def("writeLocalReadGraphReads",
            &Assembler::writeLocalReadGraphReads,
            arg("readId"),
            arg("strand"),
            arg("maxDistance"),
            arg("allowChimericReads"),
            arg("allowCrossStrandEdges"),
            arg("allowInconsistentAlignmentEdges"))
        .def("removeReadGraphBridges",
             &Assembler::removeReadGraphBridges,
             arg("maxDistance"))
        .def("analyzeReadGraph",
             &Assembler::analyzeReadGraph)
        .def("readGraphClustering",
             &Assembler::readGraphClustering)
        .def("writeReadGraphEdges",
             &Assembler::writeReadGraphEdges,
             arg("useReadName") = false)
        .def("flagInconsistentAlignments",
             &Assembler::flagInconsistentAlignments,
             arg("triangleErrorThreshold"),
             arg("leastSquareErrorThreshold"),
             arg("leastSquareMaxDistance"),
             arg("threadCount") = 0)



        // Global marker graph.
        .def("createMarkerGraphVertices",
            &Assembler::createMarkerGraphVertices,
            arg("minCoverage"),
            arg("maxCoverage"),
            arg("minCoveragePerStrand"),
            arg("allowDuplicateMarkers"),
            arg("peakFinderMinAreaFraction"),
            arg("peakFinderAreaStartIndex"),
            arg("threadCount") = 0)
        .def("accessMarkerGraphVertices",
             &Assembler::accessMarkerGraphVertices,
             arg("readWriteAccess") = false)
        .def("getGlobalMarkerGraphVertex",
            (
                MarkerGraph::VertexId (Assembler::*)
                (ReadId, Strand, uint32_t) const
            )
            &Assembler::getGlobalMarkerGraphVertex,
            arg("readId"),
            arg("strand"),
            arg("ordinal"))
        .def("getGlobalMarkerGraphVertexMarkers",
            (
                vector< std::tuple<ReadId, Strand, uint32_t> > (Assembler::*)
                (MarkerGraph::VertexId) const
            )
            &Assembler::getGlobalMarkerGraphVertexMarkers,
            arg("vertexId"))
        .def("findMarkerGraphReverseComplementVertices",
            &Assembler::findMarkerGraphReverseComplementVertices,
            arg("threadCount") = 0)
        .def("accessMarkerGraphReverseComplementVertex",
            &Assembler::accessMarkerGraphReverseComplementVertex,
            arg("readWriteAccess") = false)
        .def("findMarkerGraphReverseComplementEdges",
            &Assembler::findMarkerGraphReverseComplementEdges,
            arg("threadCount") = 0)
        .def("accessMarkerGraphReverseComplementEdge",
            &Assembler::accessMarkerGraphReverseComplementEdge)
        .def("checkMarkerGraphIsStrandSymmetric",
            &Assembler::checkMarkerGraphIsStrandSymmetric,
            arg("threadCount") = 0)
        .def("computeMarkerGraphCoverageHistogram",
            &Assembler::computeMarkerGraphCoverageHistogram)
        .def("refineMarkerGraph",
            &Assembler::refineMarkerGraph,
            arg("refineThreshold"),
            arg("threadCount") = 0)
        .def("writeBadMarkerGraphVertices",
            &Assembler::writeBadMarkerGraphVertices)
        .def("cleanupDuplicateMarkers",
            &Assembler::cleanupDuplicateMarkers,
            arg("threadCount") = 0,
            arg("minCoverage"),
            arg("minCoveragePerStrand"),
            arg("duplicateMarkersPattern1Threshold"),
            arg("pattern1CreateNewVertices"),
            arg("pattern2CreateNewVertices"))
        .def("getMarkerGraphMinCoverageUsed",
            &Assembler::getMarkerGraphMinCoverageUsed)
        .def("vertexCoverageStatisticsByKmerId",
            &Assembler::vertexCoverageStatisticsByKmerId)
        .def("analyzeMarkerGraphForks",
            &Assembler::analyzeMarkerGraphForks)

        // Edges of the global marker graph.
        .def("createMarkerGraphEdges",
            &Assembler::createMarkerGraphEdges,
            arg("threadCount") = 0)
        .def("createMarkerGraphEdgesStrict",
            &Assembler::createMarkerGraphEdgesStrict,
            arg("minEdgeCoverage"),
            arg("minEdgeCoveragePerStrand"),
            arg("threadCount") = 0)
        .def("createMarkerGraphSecondaryEdges",
            &Assembler::createMarkerGraphSecondaryEdges,
            arg("minEdgeCoverage"),
            arg("minEdgeCoveragePerStrand"),
            arg("neighborhoodSize"),
            arg("threadCount") = 0)
        .def("accessMarkerGraphEdges",
            &Assembler::accessMarkerGraphEdges,
            arg("accessEdgesReadWrite") = false,
            arg("accessConnectivityReadWrite") = false)
            .def("transitiveReduction",
            &Assembler::transitiveReduction,
            arg("lowCoverageThreshold"),
            arg("highCoverageThreshold"),
            arg("maxDistance"),
            arg("edgeMarkerSkipThreshold"))
        .def("reverseTransitiveReduction",
            &Assembler::reverseTransitiveReduction,
            arg("lowCoverageThreshold"),
            arg("highCoverageThreshold"),
            arg("maxDistance"))
        .def("pruneMarkerGraphStrongSubgraph",
            &Assembler::pruneMarkerGraphStrongSubgraph,
            arg("iterationCount"))
        .def("simplifyMarkerGraph",
            &Assembler::simplifyMarkerGraph,
            arg("maxLength"),
            arg("debug") = false)
        .def("assembleMarkerGraphVertices",
            &Assembler::assembleMarkerGraphVertices,
            arg("threadCount") = 0)
        .def("accessMarkerGraphVertexRepeatCounts",
            &Assembler::accessMarkerGraphVertexRepeatCounts)
        .def("computeMarkerGraphVerticesCoverageData",
            &Assembler::computeMarkerGraphVerticesCoverageData,
            arg("threadCount") = 0)
        .def("assembleMarkerGraphEdges",
            &Assembler::assembleMarkerGraphEdges,
            arg("threadCount") = 0,
            arg("markerGraphEdgeLengthThresholdForConsensus"),
            arg("storeCoverageData"))
        .def("accessMarkerGraphConsensus",
            &Assembler::accessMarkerGraphConsensus)
        .def("accessMarkerGraphCoverageData",
            &Assembler::accessMarkerGraphCoverageData)
        .def("setMarkerGraphEdgeFlags",
            &Assembler::setMarkerGraphEdgeFlags,
            arg("wasRemovedByTransitiveReduction"),
            arg("wasPruned"),
            arg("isSuperBubbleEdge"),
            arg("isLowCoverageCrossEdge"),
            arg("wasAssembled"))
        .def("writeParallelMarkerGraphEdges",
            &Assembler::writeParallelMarkerGraphEdges)

        // Assembly graph.
        .def("createAssemblyGraphEdges",
            &Assembler::createAssemblyGraphEdges)
        .def("createAssemblyGraphVertices",
            &Assembler::createAssemblyGraphVertices)
        .def("accessAssemblyGraphEdgeLists",
            &Assembler::accessAssemblyGraphEdgeLists)
        .def("accessAssemblyGraphEdges",
            &Assembler::accessAssemblyGraphEdges)
        .def("accessAssemblyGraphOrientedReadsByEdge",
            &Assembler::accessAssemblyGraphOrientedReadsByEdge)
        .def("accessAssemblyGraphVertices",
            &Assembler::accessAssemblyGraphVertices)
        .def("writeAssemblyGraph",
            &Assembler::writeAssemblyGraph)
        .def("assemble",
            &Assembler::assemble,
            arg("threadCount") = 0,
            arg("storeCoverageDataCsvLengthThreshold") = 0)
        .def("accessAssemblyGraphSequences",
            &Assembler::accessAssemblyGraphSequences)
        .def("computeAssemblyStatistics",
            &Assembler::computeAssemblyStatistics)
        .def("writeGfa1",
            &Assembler::writeGfa1,
            arg("fileName"))
        .def("writeGfa1BothStrands",
            &Assembler::writeGfa1BothStrands,
            arg("fileName"))
        .def("writeFasta",
            &Assembler::writeFasta,
            arg("fileName"))
        .def("colorGfaWithTwoReads",
            &Assembler::colorGfaWithTwoReads,
            arg("readId0"),
            arg("strand0"),
            arg("readId1"),
            arg("strand1"),
            arg("fileName") = "Assembly-BothStrands-Color.csv")
        .def("colorGfaKeySegments",
            &Assembler::colorGfaKeySegments,
            arg("fileName") = "Assembly-BothStrands-Color.csv")
        .def("writeOrientedReadPath",
            &Assembler::writeOrientedReadPath,
            arg("readId"),
            arg("strand"),
            arg("fileName") = "OrientedReadPath.csv")
        .def("writePseudoPath",
            &Assembler::writePseudoPath,
            arg("readId"),
            arg("strand"))
        .def("colorGfaBySimilarityToSegment",
            &Assembler::colorGfaBySimilarityToSegment,
            arg("segmentId"),
            arg("minVertexCount"),
            arg("minEdgeCount"))
        .def("assembleAssemblyGraphEdge",
            (
                AssembledSegment (Assembler::*)
                (AssemblyGraph::EdgeId, bool)
            )
            &Assembler::assembleAssemblyGraphEdge,
            arg("edgeId"),
            arg("storeCoverageData") = true)
        .def("gatherOrientedReadsByAssemblyGraphEdge",
            &Assembler::gatherOrientedReadsByAssemblyGraphEdge,
            arg("threadCount") = 0)
        .def("writeOrientedReadsByAssemblyGraphEdge",
            &Assembler::writeOrientedReadsByAssemblyGraphEdge)
        .def("detangle",
            &Assembler::detangle)
        .def("detangle2",
            &Assembler::detangle2,
            arg("diagonalReadCountMin"),
            arg("offDiagonalReadCountMax"),
            arg("offDiagonalRatio"))
        .def("alignPseudoPaths",
            &Assembler::alignPseudoPaths)
        .def("analyzeAssemblyGraphBubbles",
            &Assembler::analyzeAssemblyGraphBubbles,
            arg("debug") = true)



            // Consensus caller.
        .def("setupConsensusCaller",
            &Assembler::setupConsensusCaller)

        // CompressedAssemblyGraph.
        .def("createCompressedAssemblyGraph",
            &Assembler::createCompressedAssemblyGraph)
        .def("colorCompressedAssemblyGraph",
            &Assembler::colorCompressedAssemblyGraph,
            arg("gfaId"))



        // Mode 1 assembly.
        .def("createMode1AssemblyGraph",
            &Assembler::createMode1AssemblyGraph,
            arg("minEdgeCoverage"),
            arg("minEdgeCoveragePerStrand"))



        .def("test", &Assembler::test)

        // Definition of class_Assembler ends here.
    ;



    // Expose class AssembledSegment to Python.
    class_<AssembledSegment>(module, "AssembledSegment")
        .def("size", &AssembledSegment::size)
        .def("getBase", &AssembledSegment::getBase)
        .def("getRepeatCount", &AssembledSegment::getRepeatCount)
        .def("getCoverageData", &AssembledSegment::getCoverageData)
        ;



    // Expose class CompressedCoverageData to Python.
    class_<CompressedCoverageData>(module, "CompressedCoverageData")
        .def("getBase", &CompressedCoverageData::getBase)
        .def("getStrand", &CompressedCoverageData::getStrand)
        .def("getRepeatCount", &CompressedCoverageData::getRepeatCount)
        .def("getFrequency", &CompressedCoverageData::getFrequency)
        ;



    // Constants.
    module.attr("invalidGlobalMarkerGraphVertexId") = MarkerGraph::invalidVertexId;
    module.attr("invalidCompressedGlobalMarkerGraphVertexId") =
        uint64_t(MarkerGraph::invalidCompressedVertexId);



    // Non-member functions exposed to Python.
    module.def("testMultithreadedObject",
        testMultithreadedObject
        );
    module.def("testMemoryMappedVector",
        testMemoryMappedVector
        );
    module.def("testBase",
        testBase
        );
    module.def("testShortBaseSequence",
        testShortBaseSequence
        );
    module.def("testLongBaseSequence",
        testLongBaseSequence
        );
    module.def("testSplitRange",
        testSplitRange
        );
    module.def("testCompactUndirectedGraph1",
        testCompactUndirectedGraph1
        );
    module.def("testCompactUndirectedGraph2",
        testCompactUndirectedGraph1
        );
    module.def("testSpoa",
        testSpoa
        );
    module.def("testSimpleBayesianConsensusCaller",
        testSimpleBayesianConsensusCaller
        );
    module.def("testMedianConsensusCaller",
        testMedianConsensusCaller
        );
    module.def("testDeduplicateAndCount",
        testDeduplicateAndCount
        );
    module.def("dset64Test",
        dset64Test,
        arg("n"),
        arg("m"),
        arg("threadCount"),
        arg("batchSize"),
        arg("seed")
        );
    module.def("mappedCopy",
        mappedCopy
        );
    module.def("testAlignmentCompression",
        testAlignmentCompression
        );
    module.def("testMap",
        testMap
        );
    module.def("testMemoryMappedAllocator",
        MemoryMapped::testMemoryMappedAllocator
        );
    module.def("testLapack",
        testLapack
        );
    module.def("testShortestPathBoundedDistance",
        testShortestPathBoundedDistance
        );
}

#endif
