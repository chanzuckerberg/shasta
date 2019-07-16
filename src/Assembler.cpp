#include "Assembler.hpp"
#include "buildId.hpp"
#include "BiasedGaussianConsensusCaller.hpp"
#include "SimpleConsensusCaller.hpp"
#include "SimpleBayesianConsensusCaller.hpp"
#include "MedianConsensusCaller.hpp"
using namespace shasta;



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



// The ConsensusCaller used to compute the "best"
// base and repeat count at each assembly position.
// The argument to setupConsensusCaller specifies
// the consensus caller to be used. It can be one of the following:
// - The type of a consensus caller that
//   does not require any configuration information. Possibilities are:
//   * "SimpleConsensusCaller".
//   * "BiasedGaussianConsensusCaller".
//   * "MedianConsensusCaller".
// - The absolute path to a configuration file for the
//   consensus caller to be used. A relative path is not accepted.
//   The file name portion of this path
//   must begin with the type of the consensus caller followed by a dash.
//   Currently, the only such type of consensus caller is
//   SimpleBayesianConsensusCaller, so this requires
//   an absolute path of the form "/*/SimpleBayesianConsensusCaller-*",
//   where the two "*" can be replaced by anything.
void Assembler::setupConsensusCaller(const string& s)
{
    // Types that don't require a configuration file.
    if(s == "SimpleConsensusCaller") {
        consensusCaller = std::make_shared<SimpleConsensusCaller>();
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



    // If getting here, the argument  must be an absolute path.
    if(s.size() == 0 || s[0] != '/') {
        throw runtime_error("Invalid absolute path for consensus caller "
            "configuration file " + s);
    }

    // Extract the file name portion of the path.
    const size_t lastSlashPosition = s.find_last_of('/');
    const size_t fileNameBegin = lastSlashPosition + 1;
    const string fileName = s.substr(fileNameBegin);

    // The consensus caller type is the portion up to the dash.
    const size_t dashPosition = fileName.find_first_of('-');
    if(dashPosition == string::npos) {
        throw runtime_error("Invalid configuration file name for consensus caller. "
            "The file name portion must begin with the consensus caller type "
            "followed by a dash. " +  s);
    }
    const string consensusCallerType = fileName.substr(0, dashPosition);



    // Now that we know the type we can create it.
    if(consensusCallerType == "SimpleBayesianConsensusCaller") {
        consensusCaller = std::make_shared<SimpleBayesianConsensusCaller>(s);
        return;
    }



    // If getting here, the argument does not specify a supported
    // consensus caller.
    throw runtime_error("Unsupported consensus caller type " + consensusCallerType);
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

