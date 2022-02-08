#include "Assembler.hpp"
#include "mode3.hpp"
#include "Reads.hpp"
using namespace shasta;
using namespace mode3;



void Assembler::mode3Assembly(
    size_t threadCount)
{
    DynamicAssemblyGraph g(reads->getFlags(), markers, markerGraph, threadCount);

}
