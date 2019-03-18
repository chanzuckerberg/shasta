#include "Assembler.hpp"
#include "BiasedGaussianConsensusCaller.hpp"
#include "SimpleConsensusCaller.hpp"
#include "SimpleBayesianConsensusCaller.hpp"
#include "TrainedBayesianConsensusCaller.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



// Constructor to be called one to create a new run.
Assembler::Assembler(
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize,
    bool useRunLengthReads) :
    MultithreadedObject(*this),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize),
    marginPhaseParameters(0)
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
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize) :
    MultithreadedObject(*this),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize),
    marginPhaseParameters(0)
{

    assemblerInfo.accessExistingReadWrite(largeDataName("Info"));

    // assemblerInfo is the only open object
    // when the constructor finishes.

    fillServerFunctionTable();

}



// Destructor.
Assembler::~Assembler()
{
    if(marginPhaseParameters) {
        destroyConsensusParameters(marginPhaseParameters);
        marginPhaseParameters = 0;
    }
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

    if(s == "BiasedGaussianConsensusCaller") {
        consensusCaller = std::make_shared<BiasedGaussianConsensusCaller>();
        return;
    }

    // If getting here, the argument does not specify a supported
    // consensus caller.
    throw runtime_error("Unsupported consensus caller " + s);
}



// Read marginPhase parameters from file MarginPhase.json in the run directory.
void Assembler::setupMarginPhase()
{
    const string fileName = "MarginPhase.json";
    marginPhaseParameters = getConsensusParameters(const_cast<char*>(fileName.c_str()));
    if(!marginPhaseParameters) {
        throw runtime_error("Error reading marginPhase parameters from " + fileName);
    }
}
void Assembler::checkMarginPhaseWasSetup()
{
    if(!marginPhaseParameters) {
        throw runtime_error("MarginPhase was not set up.");
    }
}

