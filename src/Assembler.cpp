#include "Assembler.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



// Constructor to be called one to create a new run.
Assembler::Assembler(
    const string& smallDataFileNamePrefix,
    const string& largeDataFileNamePrefix,
    size_t smallDataPageSize,
    size_t largeDataPageSize,
    bool useRunLengthReads) :
    MultithreadedObject(*this),
    smallDataFileNamePrefix(smallDataFileNamePrefix),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    smallDataPageSize(smallDataPageSize),
    largeDataPageSize(largeDataPageSize)
{
    assemblerInfo.createNew(smallDataName("Info"));
    assemblerInfo->useRunLengthReads = useRunLengthReads;

    reads.createNew(largeDataName("Reads"), largeDataPageSize);
    reads.close();

    readNames.createNew(largeDataName("ReadNames"), largeDataPageSize);
    readNames.close();

    if(useRunLengthReads) {
        readRepeatCounts.createNew(largeDataName("ReadRepeatCounts"), largeDataPageSize);
        readRepeatCounts.close();

    }

    // assemblerInfo is the only open object
    // when the constructor finishes.

    fillServerFunctionTable();

}



// Constructor to be called to continue an existing run.
Assembler::Assembler(
    const string& smallDataFileNamePrefix,
    const string& largeDataFileNamePrefix,
    size_t smallDataPageSize,
    size_t largeDataPageSize) :
    MultithreadedObject(*this),
    smallDataFileNamePrefix(smallDataFileNamePrefix),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    smallDataPageSize(smallDataPageSize),
    largeDataPageSize(largeDataPageSize)
{

    assemblerInfo.accessExistingReadWrite(smallDataName("Info"));

    // assemblerInfo is the only open object
    // when the constructor finishes.

    fillServerFunctionTable();

}

