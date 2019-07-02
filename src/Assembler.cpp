#include "Assembler.hpp"
#include "buildId.hpp"
#include "BiasedGaussianConsensusCaller.hpp"
#include "SimpleConsensusCaller.hpp"
#include "SimpleBayesianConsensusCaller.hpp"
#include "MedianConsensusCaller.hpp"
using namespace ::shasta;



// Constructor to be called one to create a new run.
Assembler::Assembler(
    const string& largeDataFileNamePrefix,
    bool createNew,
    size_t largeDataPageSize) :

    MultithreadedObject(*this),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize)
#ifndef SHASTA_STATIC_EXECUTABLE
    , marginPhaseParameters(0)
#endif
{
    cout << buildId() << endl;

    if(createNew) {

        // Create a new assembly.
        assemblerInfo.createNew(largeDataName("Info"), largeDataPageSize);
        assemblerInfo->largeDataPageSize = largeDataPageSize;

        reads.createNew(largeDataName("Reads"), largeDataPageSize);
        readNames.createNew(largeDataName("ReadNames"), largeDataPageSize);
        readRepeatCounts.createNew(largeDataName("ReadRepeatCounts"), largeDataPageSize);

    } else {

        // Access an existing assembly.
        assemblerInfo.accessExistingReadWrite(largeDataName("Info"));
        largeDataPageSize = assemblerInfo->largeDataPageSize;

        reads.accessExistingReadWrite(largeDataName("Reads"));
        readNames.accessExistingReadWrite(largeDataName("ReadNames"));
        readRepeatCounts.accessExistingReadWrite(largeDataName("ReadRepeatCounts"));

    }

    // In both cases, assemblerInfo, reads, readNames, readRepeatCounts are all open for write.

#ifndef SHASTA_STATIC_EXECUTABLE
    fillServerFunctionTable();
#endif
}



// Destructor.
Assembler::~Assembler()
{
#ifndef SHASTA_STATIC_EXECUTABLE
    if(marginPhaseParameters) {
        destroyConsensusParameters(marginPhaseParameters);
        marginPhaseParameters = 0;
    }
#endif
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

    if(s == "BiasedGaussianConsensusCaller") {
        consensusCaller = std::make_shared<BiasedGaussianConsensusCaller>();
        return;
    }

    if(s == "MedianConsensusCaller") {
        consensusCaller = std::make_shared<MedianConsensusCaller>();
        return;
    }


    // If getting here, the argument does not specify a supported
    // consensus caller.
    throw runtime_error("Unsupported consensus caller " + s);
}



// Read marginPhase parameters from file MarginPhase.json in the run directory.
void Assembler::setupMarginPhase()
{
#ifndef SHASTA_STATIC_EXECUTABLE
    const string fileName = "MarginPhase.json";
    marginPhaseParameters = getConsensusParameters(const_cast<char*>(fileName.c_str()));
    if(!marginPhaseParameters) {
        throw runtime_error("Error reading marginPhase parameters from " + fileName);
    }
#else
    // The static executable does not support MarginPhase.
    SHASTA_ASSERT(0);
#endif
}



void Assembler::checkMarginPhaseWasSetup()
{
#ifndef SHASTA_STATIC_EXECUTABLE
    if(!marginPhaseParameters) {
        throw runtime_error("MarginPhase was not set up.");
    }
#else
    // The static executable does not support MarginPhase.
    SHASTA_ASSERT(0);
#endif
}


// Store assembly time.
void Assembler::storeAssemblyTime(
    double elapsedTimeSeconds,
    double averageCpuUtilization)
{
    assemblerInfo->assemblyElapsedTimeSeconds = elapsedTimeSeconds;
    assemblerInfo->averageCpuUtilization = averageCpuUtilization;
}

