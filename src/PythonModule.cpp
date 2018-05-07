
// Nanopore2.
#include "Assembler.hpp"
#include "Base.hpp"
#include "LongBaseSequence.hpp"
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



PYBIND11_MODULE(Nanopore2, module)
{
    class_<Assembler>(module, "Assembler")

        // Constructor.
        .def(init<string, string, size_t, size_t>(),
            "Access an existing Assembler or create a new one.",
            arg("smallDataFileNamePrefix") = "data/",
            arg("largeDataFileNamePrefix") = "Data/",
            arg("smallDataPageSize") = 4096,
            arg("largeDataPageSize") = 2*1024*1024)

        .def("addReadsFromFasta",
            &Assembler::addReadsFromFasta,
            "Add reads from a fasta file.",
            arg("fileName"),
            arg("blockSize") = 64 * 1024 * 1024,
            arg("threadCountForReading") = 1,
            arg("threadCountForProcessing") = 0)

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

}

