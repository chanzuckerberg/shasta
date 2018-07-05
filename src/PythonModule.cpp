
// Nanopore2.
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
using namespace Nanopore2;

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
using namespace pybind11;



PYBIND11_MODULE(shasta, module)
{
    class_<Assembler>(module, "Assembler")

        // Constructor.
        .def(init<string, string, size_t, size_t>(),
            "Access an existing Assembler or create a new one.",
            arg("smallDataFileNamePrefix") = "data/",
            arg("largeDataFileNamePrefix") = "Data/",
            arg("smallDataPageSize") = 4096,
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



        // Overlaps
        .def("findOverlaps",
            &Assembler::findOverlaps,
            arg("m"),
            arg("minHashIterationCount"),
            arg("log2MinHashBucketCount"),
            arg("maxBucketSize"),
            arg("minFrequency"),
            arg("threadCount") = 0)
        .def("accessOverlaps",
            &Assembler::accessOverlaps)
        .def("writeOverlappingReads",
            &Assembler::writeOverlappingReads,
            "Write in fasta format the reads that overlap a given read.",
            arg("readId"),
            arg("strand"),
            arg("fileName") = "OverlappingReads.fasta")

        // Read graph.
        .def("computeReadGraphComponents",
            &Assembler::computeReadGraphComponents,
            arg("minComponentSize"),
            arg("minAlignedMarkerCount"),
            arg("maxTrim"))
        .def("createLocalReadGraph",
            (
                void (Assembler::*)
                (ReadId, Strand, size_t, size_t, uint32_t)
            )
            &Assembler::createLocalReadGraph,
            arg("readId"),
            arg("strand"),
            arg("minAlignedMarkerCount"),
            arg("maxTrim"),
            arg("distance"))



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



        // Compute all alignments and, optionally, the global marker graph.
        .def("computeAllAlignments",
            &Assembler::computeAllAlignments,
            arg("maxVertexCountPerKmer"),
            arg("maxSkip"),
            arg("minAlignedMarkerCount"),
            arg("maxTrim"),
            arg("minCoverage"),
            arg("threadCount") = 0)
        .def("accessAlignmentData",
            &Assembler::accessAlignmentData)



        // Global marker graph.
        .def("accessGlobalMarkerGraph",
             &Assembler::accessGlobalMarkerGraph)
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
            arg("globalMarkerGraphVertexId"))
        .def("getGlobalMarkerGraphVertexChildren",
            (
                vector<GlobalMarkerGraphVertexId> (Assembler::*)
                (GlobalMarkerGraphVertexId) const
            )
            &Assembler::getGlobalMarkerGraphVertexChildren,
            arg("globalMarkerGraphVertexId"))
        .def("getGlobalMarkerGraphVertexParents",
            (
                vector<GlobalMarkerGraphVertexId> (Assembler::*)
                (GlobalMarkerGraphVertexId) const
            )
            &Assembler::getGlobalMarkerGraphVertexParents,
            arg("globalMarkerGraphVertexId"))
        .def("extractLocalMarkerGraph",
            (
                void (Assembler::*)
                (ReadId, Strand, uint32_t, int, size_t, size_t)
            )
            &Assembler::extractLocalMarkerGraph,
            arg("readId"),
            arg("strand"),
            arg("ordinal"),
            arg("distance"),
            arg("minCoverage"),
            arg("minConsensus"))



        // Http server.
        .def("explore",
            &Assembler::explore,
            arg("port") = 17100,
            arg("localOnly") = false);



        // Definition of class_<Assembler> ends here.
    ;



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

