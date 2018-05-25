#include "Assembler.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;



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
    try {
        assemblerInfo.accessExistingReadWrite(smallDataName("Info"));
    } catch(runtime_error) {
        assemblerInfo.createNew(smallDataName("Info"));
        reads.createNew(largeDataName("Reads"), largeDataPageSize);
        readNames.createNew(largeDataName("ReadNames"), largeDataPageSize);
        reads.close();
        readNames.close();
    }

    // Either way, assemblerInfo is the only open object
    // when the constructor finishes.
}

