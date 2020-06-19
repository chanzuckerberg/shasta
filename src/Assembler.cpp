#include "Assembler.hpp"
#include "buildId.hpp"
#include "SimpleConsensusCaller.hpp"
#include "SimpleBayesianConsensusCaller.hpp"
#include "MedianConsensusCaller.hpp"
using namespace shasta;



// Constructor to be called one to create a new run.
Assembler::Assembler(
    const string& largeDataFileNamePrefix,
    bool createNew,
    size_t largeDataPageSizeArgument) :

    MultithreadedObject(*this),
    largeDataFileNamePrefix(largeDataFileNamePrefix)
{

    if(createNew) {

        // Create a new assembly.
        assemblerInfo.createNew(largeDataName("Info"), largeDataPageSizeArgument);
        assemblerInfo->largeDataPageSize = largeDataPageSizeArgument;
        largeDataPageSize = largeDataPageSizeArgument;

        reads.createNew(largeDataName("Reads"), largeDataPageSize);
        readNames.createNew(largeDataName("ReadNames"), largeDataPageSize);
        readMetaData.createNew(largeDataName("ReadMetaData"), largeDataPageSize);
        readRepeatCounts.createNew(largeDataName("ReadRepeatCounts"), largeDataPageSize);
        // cout << "Created a new assembly with page size " << largeDataPageSize << endl;

    } else {

        // Access an existing assembly.
        assemblerInfo.accessExistingReadWrite(largeDataName("Info"));
        largeDataPageSize = assemblerInfo->largeDataPageSize;

        reads.accessExistingReadWrite(largeDataName("Reads"));
        readNames.accessExistingReadWrite(largeDataName("ReadNames"));
        readMetaData.accessExistingReadWrite(largeDataName("ReadMetaData"));
        readRepeatCounts.accessExistingReadWrite(largeDataName("ReadRepeatCounts"));
        // cout << "Accessed an existing assembly with page size " << largeDataPageSize << endl;

    }
    SHASTA_ASSERT(largeDataPageSize == assemblerInfo->largeDataPageSize);

    // In both cases, assemblerInfo, reads, readNames, readRepeatCounts are all open for write.

#ifdef SHASTA_HTTP_SERVER
    fillServerFunctionTable();
#endif
}



// Set up the ConsensusCaller used to compute the "best"
// base and repeat count at each assembly position.
// The argument to setupConsensusCaller specifies
// the consensus caller to be used.
// It can be one of the following:
// - Modal
//   Selects the SimpleConsensusCaller.
// - Median
//   Selects the MedianConsensusCaller.
// - Bayesian:fileName
//   Selects the SimpleBayesianConsensusCaller,
//   using fileName as the configuration file.
//   Filename must be an absolute path (it must begin with "/").
void Assembler::setupConsensusCaller(const string& s)
{

    // Parse the argument.
    // The portion up to the first colon is the type string (Modal, Median, or Bayesian).
    // The portion following the first colon is the constructor string.
    const size_t colonPosition = s.find_first_of(':');
    string typeString, constructorString;
    if(colonPosition==string::npos) {
        // There is no colon.
        typeString = s;
        constructorString = "";
    } else {
        // The colon is present.
        typeString =  s.substr(0, colonPosition);
        constructorString =  s.substr(colonPosition+1);
    }



    // The Modal consensus caller (class SimpleConsensusCaller)
    // takes no constructor string (nothing after the colon).
    if(typeString == "Modal") {

        //  The constructor string must be empty.
        if(constructorString.size() > 0) {
            throw runtime_error("Invalid consensus caller " + s);
        }

        consensusCaller = std::make_shared<SimpleConsensusCaller>();
        return;
    }



    // The Median consensus caller (class MedianConsensusCaller)
    // takes no  constructor string (nothing after the colon).
    if(typeString == "Median") {

        //  The constructor string must be empty.
        if(constructorString.size() > 0) {
            throw runtime_error("Invalid consensus caller " + s);
        }

        consensusCaller = std::make_shared<MedianConsensusCaller>();
        return;
    }



    // The Bayesian consensus caller (class SimpleBayesianConsensusCaller)
    // requires a constructor string (portion following the first colon).
    // The constructor string is either the name of a built-in Bayesian model
    // or the absolute path to the configuration file.
    if(typeString == "Bayesian") {
        if(constructorString.size() == 0) {
            throw runtime_error("Invalid consensus caller " + s +
                ". Bayesian:builtinName required or Bayesian::absolutePath");
        }

        consensusCaller = std::make_shared<SimpleBayesianConsensusCaller>(constructorString);
        return;
    }



    // If getting here, the argument does not specify a supported
    // consensus caller.
    throw runtime_error("Invalid consensus caller " + s +
        ". Valid choices are: Modal, Median, Bayesian:absolutePathToConfigFile.");
}



// Store assembly time.
void Assembler::storeAssemblyTime(
    double elapsedTimeSeconds,
    double averageCpuUtilization)
{
    assemblerInfo->assemblyElapsedTimeSeconds = elapsedTimeSeconds;
    assemblerInfo->averageCpuUtilization = averageCpuUtilization;
}

void Assembler::storePeakMemoryUsage(const uint64_t peakMemoryUsage) {
    assemblerInfo->peakMemoryUsage = peakMemoryUsage;
}

