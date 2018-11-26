
// shasta.
#include "Assembler.hpp"
#include "Base.hpp"
#include "CompactUndirectedGraph.hpp"
#include "dset64Test.hpp"
#include "LongBaseSequence.hpp"
#include "mappedCopy.hpp"
#include "MultitreadedObject.hpp"
#include "ShortBaseSequence.hpp"
#include "splitRange.hpp"
using namespace ChanZuckerberg;
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


    class_<Assembler>(module, "Assembler")

        // Constructors.
        .def(init<string, size_t, bool>(),
            "Create a new assembler.",
            arg("largeDataFileNamePrefix") = "Data/",
            arg("largeDataPageSize") = 2*1024*1024,
            arg("useRunLengthReads"))
        .def(init<string, size_t>(),
            "Access an existing Assembler.",
            arg("largeDataFileNamePrefix") = "Data/",
            arg("largeDataPageSize") = 2*1024*1024)



        // Reads
        .def("addReadsFromFasta",
            &Assembler::addReadsFromFasta,
            "Add reads from a fasta file.",
            arg("fileName"),
            arg("minReadLength"),
            arg("blockSize") = 2ULL * 1024ULL * 1024ULL * 1024ULL,
            arg("threadCountForReading") = 1,
            arg("threadCountForProcessing") = 0)
        .def("accessReadsReadOnly",
            &Assembler::accessReadsReadOnly)
        .def("accessReadsReadWrite",
            &Assembler::accessReadsReadWrite)
        .def("accessReadNamesReadOnly",
            &Assembler::accessReadNamesReadOnly)
        .def("accessReadNamesReadWrite",
            &Assembler::accessReadNamesReadWrite)
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



        // Alignment candidates.
        .def("findAlignmentCandidates",
            &Assembler::findAlignmentCandidates,
            arg("m"),
            arg("minHashIterationCount"),
            arg("log2MinHashBucketCount"),
            arg("maxBucketSize"),
            arg("minFrequency"),
            arg("threadCount") = 0)
        .def("accessCandidateAlignments",
            &Assembler::accessAlignmentCandidates)
        .def("writeOverlappingReads",
            &Assembler::writeOverlappingReads,
            "Write in fasta format the reads that overlap a given read.",
            arg("readId"),
            arg("strand"),
            arg("fileName") = "OverlappingReads.fasta")



        // Alignments.
        .def("alignOrientedReads",
            (
                void (Assembler::*)
                (ReadId, Strand, ReadId, Strand, size_t, size_t)
            )
            &Assembler::alignOrientedReads,
            arg("readId0"),
            arg("strand0"),
            arg("readId1"),
            arg("strand1"),
            arg("maxSkip"),
            arg("maxVertexCountPerKmer"))
        .def("alignOverlappingOrientedReads",
            (
                void (Assembler::*)
                (ReadId, Strand, size_t, size_t, size_t, size_t)
            )
            &Assembler::alignOverlappingOrientedReads,
            arg("readId"),
            arg("strand"),
            arg("maxSkip"),
            arg("maxVertexCountPerKmer"),
            arg("minAlignedMarkerCount"),
            arg("maxTrim")
            )

        // Compute an alignment for each alignment candidate.
        .def("computeAlignments",
            &Assembler::computeAlignments,
            arg("maxVertexCountPerKmer"),
            arg("maxSkip"),
            arg("minAlignedMarkerCount"),
            arg("maxTrim"),
            arg("threadCount") = 0)
        .def("accessAlignmentData",
            &Assembler::accessAlignmentData)



        // Read graph
        .def("createReadGraph",
            &Assembler::createReadGraph,
            arg("maxTrim"))
        .def("accessReadGraph",
            &Assembler::accessReadGraph)
        .def("flagChimericReads",
             &Assembler::flagChimericReads,
            arg("maxChimericReadDistance"),
            arg("threadCount") = 0)
        .def("accessChimericReadsFlags",
            &Assembler::accessChimericReadsFlags)



        // Global marker graph.
        .def("createMarkerGraphVertices",
            &Assembler::createMarkerGraphVertices,
            arg("maxVertexCountPerKmer"),
            arg("maxSkip"),
            arg("minCoverage"),
            arg("maxCoverage"),
            arg("threadCount") = 0)
        .def("accessMarkerGraphVertices",
             &Assembler::accessMarkerGraphVertices)
        .def("getGlobalMarkerGraphVertex",
            (
                GlobalMarkerGraphVertexId (Assembler::*)
                (ReadId, Strand, uint32_t) const
            )
            &Assembler::getGlobalMarkerGraphVertex,
            arg("readId"),
            arg("strand"),
            arg("ordinal"))
        .def("getGlobalMarkerGraphVertexMarkers",
            (
                vector< std::tuple<ReadId, Strand, uint32_t> > (Assembler::*)
                (GlobalMarkerGraphVertexId) const
            )
            &Assembler::getGlobalMarkerGraphVertexMarkers,
            arg("vertexId"))
        .def("getGlobalMarkerGraphVertexChildren",
            (
                vector<GlobalMarkerGraphVertexId> (Assembler::*)
                (GlobalMarkerGraphVertexId) const
            )
            &Assembler::getGlobalMarkerGraphVertexChildren,
            arg("vertexId"))
        .def("getGlobalMarkerGraphVertexParents",
            (
                vector<GlobalMarkerGraphVertexId> (Assembler::*)
                (GlobalMarkerGraphVertexId) const
            )
            &Assembler::getGlobalMarkerGraphVertexParents,
            arg("vertexId"))
        .def("getGlobalMarkerGraphEdgeInformation",
                &Assembler::getGlobalMarkerGraphEdgeInformation,
                arg("vertexId0"),
                arg("vertexId1"))
        .def("extractLocalMarkerGraph",
        (
            void (Assembler::*)
            (ReadId, Strand, uint32_t, int, size_t)
        )
        &Assembler::extractLocalMarkerGraph,
        arg("readId"),
        arg("strand"),
        arg("ordinal"),
        arg("distance"),
        arg("minCoverage"))
        .def("getLocalAssemblyPath",
            &Assembler::getLocalAssemblyPath,
            arg("startVertexId"),
            arg("maxDistance"))


        // Connectivity of the global marker graph.
        .def("createMarkerGraphConnectivity",
             &Assembler::createMarkerGraphConnectivity,
             arg("threadCount") = 0)
         .def("accessMarkerGraphConnectivity",
              &Assembler::accessMarkerGraphConnectivity,
              arg("accessEdgesReadWrite") = false)
         .def("flagMarkerGraphWeakEdges",
              &Assembler::flagMarkerGraphWeakEdges,
                 arg("threadCount") = 0,
                 arg("minCoverage"),
                 arg("maxPathLength"))
         .def("computeMarkerGraphSpanningSubgraph",
              &Assembler::computeMarkerGraphSpanningSubgraph,
                 arg("minCoverage"))
         .def("pruneMarkerGraphSpanningSubgraph",
              &Assembler::pruneMarkerGraphSpanningSubgraph,
                 arg("iterationCount"))



         // Assembly graph.
         .def("createAssemblyGraphVertices",
              &Assembler::createAssemblyGraphVertices)
         .def("createAssemblyGraphEdges",
             &Assembler::createAssemblyGraphEdges)
         .def("accessAssemblyGraphVertices",
              &Assembler::accessAssemblyGraphVertices)
          .def("accessAssemblyGraphEdges",
               &Assembler::accessAssemblyGraphEdges)



        // Http server.
        .def("explore",
            &Assembler::explore,
            arg("port") = 17100,
            arg("localOnly") = false)
        .def("setDocsDirectory",
            &Assembler::setDocsDirectory)
        .def("setReferenceFastaFileName",
            &Assembler::setReferenceFastaFileName)

        // Consensus caller.
        .def("setupConsensusCaller",
            &Assembler::setupConsensusCaller)



        // Definition of class_<Assembler> ends here.
    ;



    // Constants.
    module.attr("invalidGlobalMarkerGraphVertexId") = invalidGlobalMarkerGraphVertexId;
    module.attr("invalidCompressedGlobalMarkerGraphVertexId") =
        uint64_t(invalidCompressedGlobalMarkerGraphVertexId);



    // Non-member functions exposed to Python.
    module.def("testMultithreadedObject",
        testMultithreadedObject
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

}

