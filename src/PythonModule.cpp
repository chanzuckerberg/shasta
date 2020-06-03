#ifdef SHASTA_PYTHON_API

// Shasta.
#include "Assembler.hpp"
#include "Base.hpp"
#include "CompactUndirectedGraph.hpp"
#include "compressAlignment.hpp"
#include "deduplicate.hpp"
#include "dset64Test.hpp"
#include "LongBaseSequence.hpp"
#include "mappedCopy.hpp"
#include "MultithreadedObject.hpp"
#include "ShortBaseSequence.hpp"
#include "splitRange.hpp"
#include "testSpoa.hpp"
#include "SimpleBayesianConsensusCaller.hpp"
#include "MedianConsensusCaller.hpp"
using namespace shasta;

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
using namespace pybind11;



PYBIND11_MODULE(shasta, module)
{

    {
        // Class used by Assembler::getGlobalMarkerGraphEdgeInformation.
        using Info = Assembler::GlobalMarkerGraphEdgeInformation;
        class_<Info>(module, "GlobalMarkerGraphEdgeInformation")
            .def_readwrite("readId", &Info::readId)
            .def_readwrite("strand", &Info::strand)
            .def_readwrite("ordinal0", &Info::ordinal0)
            .def_readwrite("ordinal1", &Info::ordinal1)
            .def_readwrite("position0", &Info::position0)
            .def_readwrite("position1", &Info::position1)
            .def_readwrite("overlappingBaseCount", &Info::overlappingBaseCount)
            .def_readwrite("sequence", &Info::sequence)
            ;
    }



    // Expose class OrientedReadPair to Python.
    class_<OrientedReadPair>(module, "OrientedReadPair")
        .def_readonly("readIds", &OrientedReadPair::readIds)
        .def_readonly("isSameStrand", &OrientedReadPair::isSameStrand)
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
        .def("readCount",
            &Assembler::readCount,
            "Get the number of reads.")
        .def("histogramReadLength",
            &Assembler::histogramReadLength,
            "Create a histogram of read length and write it to a csv file.",
            arg("fileName") = "ReadLengthHistogram.csv")
        .def("writeReads",
            &Assembler::writeReads,
            "Write all reads to a file in fasta format.",
            arg("fileName"))
        .def("writeRead",
            (
                void (Assembler::*)
                (ReadId, const string&)
            )
            &Assembler::writeRead,
            "Write one read to a file in fasta format.",
            arg("readId"),
            arg("fileName"))
        .def("writeOrientedRead",
            (
                void (Assembler::*)
                (ReadId, Strand, const string&)
            )
            &Assembler::writeOrientedRead,
            "Write one oriented read to a file in fasta format.",
            arg("readId"),
            arg("strand"),
            arg("fileName"))
        .def("initializeReadFlags",
            &Assembler::initializeReadFlags)
        .def("accessReadFlags",
            &Assembler::accessReadFlags,
            arg("readWriteAccess") = false)



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
            arg("threadCount") = 0)

        // Alignments.
        .def("writeAlignmentCandidates",
            &Assembler::writeAlignmentCandidates)
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
            &Assembler::computeAlignments,
            arg("alignmentMethod") = 0,
            arg("maxMarkerFrequency"),
            arg("maxSkip"),
            arg("maxDrift"),
            arg("minAlignedMarkerCount"),
            arg("minAlignedFraction"),
            arg("maxTrim"),
            arg("matchScore"),
            arg("mismatchScore"),
            arg("gapScore"),
            arg("downsamplingFactor"),
            arg("bandExtend"),
            arg("suppressContainments"),
            arg("storeAlignments"),
            arg("threadCount") = 0)
        .def("accessCompressedAlignments",
            &Assembler::accessCompressedAlignments)
        .def("accessAlignmentData",
            &Assembler::accessAlignmentData)
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
            arg("allowCrossStrandEdges"))



        // Directed read graph.
        .def("createDirectedReadGraph",
            &Assembler::createDirectedReadGraph,
            arg("maxTrim"),
            arg("containedNeighborCount"),
            arg("uncontainedNeighborCountPerDirection"))
        .def("accessDirectedReadGraphReadOnly",
            &Assembler::accessDirectedReadGraphReadOnly)
        .def("accessDirectedReadGraphReadWrite",
            &Assembler::accessDirectedReadGraphReadWrite)
        .def("markDirectedReadGraphConflictEdges1",
            &Assembler::markDirectedReadGraphConflictEdges1)
        .def("markDirectedReadGraphConflictEdges2",
            &Assembler::markDirectedReadGraphConflictEdges2,
            arg("radius"))
        .def("markDirectedReadGraphConflictEdges3",
            &Assembler::markDirectedReadGraphConflictEdges3,
            arg("radius"))


        // Global marker graph.
        .def("createMarkerGraphVertices",
            &Assembler::createMarkerGraphVertices,
            arg("alignMethod"),
            arg("maxMarkerFrequency"),
            arg("maxSkip"),
            arg("maxDrift"),
            arg("matchScore"),
            arg("mismatchScore"),
            arg("gapScore"),
            arg("downsamplingFactor"),
            arg("bandExtend"),
            arg("readGraphCreationMethod"),
            arg("minCoverage"),
            arg("maxCoverage"),
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
        .def("getGlobalMarkerGraphVertexChildren",
            (
                vector<MarkerGraph::VertexId> (Assembler::*)
                (MarkerGraph::VertexId) const
            )
            &Assembler::getGlobalMarkerGraphVertexChildren,
            arg("vertexId"))
        .def("getGlobalMarkerGraphVertexParents",
            (
                vector<MarkerGraph::VertexId> (Assembler::*)
                (MarkerGraph::VertexId) const
            )
            &Assembler::getGlobalMarkerGraphVertexParents,
            arg("vertexId"))
        .def("getGlobalMarkerGraphEdgeInformation",
                &Assembler::getGlobalMarkerGraphEdgeInformation,
                arg("vertexId0"),
                arg("vertexId1"))
        .def("findMarkerGraphReverseComplementVertices",
            &Assembler::findMarkerGraphReverseComplementVertices,
            arg("threadCount") = 0)
        .def("accessMarkerGraphReverseComplementVertex",
            &Assembler::accessMarkerGraphReverseComplementVertex)
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
        .def("analyzeMarkerGraphVertex",
            &Assembler::analyzeMarkerGraphVertex)
        .def("refineMarkerGraph",
            &Assembler::refineMarkerGraph,
            arg("refineThreshold"),
            arg("threadCount") = 0)



        // Edges of the global marker graph.
        .def("createMarkerGraphEdges",
            &Assembler::createMarkerGraphEdges,
            arg("threadCount") = 0)
        .def("accessMarkerGraphEdges",
            &Assembler::accessMarkerGraphEdges,
            arg("accessEdgesReadWrite") = false)
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
        .def("createConflictReadGraph",
            &Assembler::createConflictReadGraph,
            arg("threadCount") = 0,
            arg("maxOffsetSigma"),
            arg("maxTrim"),
            arg("maxSkip"),
            arg("minAlignedMarkerCount"))
        .def("accessConflictReadGraph",
            &Assembler::accessConflictReadGraph)
        .def("cleanupConflictReadGraph",
            &Assembler::cleanupConflictReadGraph)


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
        .def("findAssemblyGraphBubbles",
            &Assembler::findAssemblyGraphBubbles)
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
        .def("createSegmentGraph",
            &Assembler::createSegmentGraph)
        .def("colorGfaBySegmentGraphChain",
            &Assembler::colorGfaBySegmentGraphChain)
        .def("analyzeOrientedReadPaths",
            &Assembler::analyzeOrientedReadPaths,
            arg("readGraphCreationMethod"))
        .def("analyzeOrientedReadPathsThroughSegment",
            &Assembler::analyzeOrientedReadPathsThroughSegment,
            arg("segmentId"))

        // Consensus caller.
        .def("setupConsensusCaller",
            &Assembler::setupConsensusCaller)

        // CompressedAssemblyGraph.
        .def("createCompressedAssemblyGraph",
            &Assembler::createCompressedAssemblyGraph)
        .def("colorCompressedAssemblyGraph",
            &Assembler::colorCompressedAssemblyGraph,
            arg("gfaId"))

        // Phasing.
        .def("createPhasingData",
            &Assembler::createPhasingData,
            arg("threadCount") = 0,
            arg("phasingSimilarityThreshold"),
            arg("maxNeighborCount"))
        .def("accessPhasingData",
            &Assembler::accessPhasingData)
#if 0
        .def("computePhasingSimilarity",
            (
                double (Assembler::*)
                (ReadId, Strand, ReadId, Strand)
            )
            &Assembler::computePhasingSimilarity,
            arg("readId0"),
            arg("strand0"),
            arg("readId1"),
            arg("strand1")
            )
#endif
        .def("computePhasingSimilarity",
            (
                double (Assembler::*)
                (AssemblyGraph::EdgeId, AssemblyGraph::EdgeId)
            )
            &Assembler::computePhasingSimilarity,
            arg("edgeId0"),
            arg("edgeId1")
            )
        .def("countCommonInternalOrientedReads",
            &Assembler::countCommonInternalOrientedReads,
            arg("edgeId0"),
            arg("edgeId1")
            )

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
}

#endif
