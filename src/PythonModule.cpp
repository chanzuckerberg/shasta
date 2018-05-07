
// Nanopore2.
#include "MultitreadedObject.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
using namespace pybind11;



PYBIND11_MODULE(Nanopore2, module)
{

    // Non-member functions exposed to Python.
    module.def("testMultithreadedObject",
        testMultithreadedObject
        );

}

