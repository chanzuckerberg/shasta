// Alternative alignment functions with 1 suffix.
#include "Assembler.hpp"
using namespace shasta;



void Assembler::alignOrientedReads1(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1)
{
    alignOrientedReads1(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1));
}



void Assembler::alignOrientedReads1(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1)
{
    SHASTA_ASSERT(0);
}


