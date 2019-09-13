#include "ReadRepeatCodes.hpp"
using namespace shasta;



void ReadRepeatCodes::createNew(const string& name, uint64_t pageSize)
{
    readRepeatCodes.createNew(name, pageSize);
}
void ReadRepeatCodes::accessExisting(const string& name)
{
    readRepeatCodes.accessExistingReadWrite(name);
}



// Function called to set the probabilistic mode of operation.
// and specify the nominal repeat count corresponding to each
// repeat code.
// If this is not called, the standard mode of operation is used.
void ReadRepeatCodes::setProbabilisticMode(
    const array<uint8_t, 256>& nominalRepeatCountArgument)
{
    nominalRepeatCount = nominalRepeatCountArgument;
}
