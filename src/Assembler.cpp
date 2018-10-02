#include "Assembler.hpp"
#include "SimpleConsensusCaller.hpp"
#include "SimpleBayesianConsensusCaller.hpp"
#include "TrainedBayesianConsensusCaller.hpp"
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
    assemblerInfo.createNew(largeDataName("Info"), largeDataPageSize);
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

    assemblerInfo.accessExistingReadWrite(largeDataName("Info"));

    // assemblerInfo is the only open object
    // when the constructor finishes.

    fillServerFunctionTable();

}


// Set up the ConsensusCaller used to compute the "best"
// base and repeat count at each assembly position.
void Assembler::setupConsensusCaller(const string& s)
{
    if(s == "SimpleConsensusCaller") {
        consensusCaller = std::make_shared<SimpleConsensusCaller>();
        return;
    }

    if(s == "SimpleBayesianConsensusCaller") {
        consensusCaller = std::make_shared<SimpleBayesianConsensusCaller>();
        return;
    }

    if(s == "TrainedBayesianConsensusCaller") {
        consensusCaller = std::make_shared<TrainedBayesianConsensusCaller>();
        return;
    }

    // If getting here, the argument does not specify a supported
    // consensus caller.
    throw runtime_error("Unsupported consensus caller " + s);
}
