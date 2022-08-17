// Main program for the Shasta static executable.
// The static executable provides
// basic functionality and reduced performance.
// For full functionality use the shared library built
// under directory src.

// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "AssemblyGraph.hpp"
#include "buildId.hpp"
#include "ConfigurationTable.hpp"
#include "Coverage.hpp"
#include "filesystem.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "Tee.hpp"
#include "timestamp.hpp"
#include "platformDependent.hpp"
#include "SimpleBayesianConsensusCaller.hpp"

// Standard library.
#include <filesystem>

namespace shasta {
    namespace main {

        void main(int argumentCount, const char** arguments);

        void assemble(
            Assembler&,
            const AssemblerOptions&,
            vector<string> inputNames);

        void mode0Assembly(
            Assembler&,
            const AssemblerOptions&,
            uint32_t threadCount);
        void mode2Assembly(
            Assembler&,
            const AssemblerOptions&,
            uint32_t threadCount);
        void mode3Assembly(
            Assembler&,
            const AssemblerOptions&,
            uint32_t threadCount);

        void setupRunDirectory(
            const string& memoryMode,
            const string& memoryBacking,
            size_t& pageSize,
            string& dataDirectory
            );

        void setupHugePages();
        void segmentFaultHandler(int);

        // Functions that implement --command keywords
        void assemble(const AssemblerOptions&, int argumentCount, const char** arguments);
        void saveBinaryData(const AssemblerOptions&);
        void cleanupBinaryData(const AssemblerOptions&);
        void createBashCompletionScript(const AssemblerOptions&);
        void listCommands();
        void listConfigurations();
        void listConfiguration(const AssemblerOptions&);
        void explore(const AssemblerOptions&);

        const std::set<string> commands = {
            "assemble",
            "cleanupBinaryData",
            "createBashCompletionScript",
            "explore",
            "listCommands",
            "listConfiguration",
            "listConfigurations",
            "saveBinaryData"};

    }

    Tee tee;
    ofstream shastaLog;
}
using namespace shasta;

// Boost libraries.
#include <boost/program_options.hpp>
#include  <boost/chrono/process_cpu_clocks.hpp>

//  Linux.
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>


// Standard library.
#include "chrono.hpp"
#include "iostream.hpp"
#include "iterator.hpp"
#include "stdexcept.hpp"



int main(int argumentCount, const char** arguments)
{
    try {

        shasta::main::main(argumentCount, arguments);

    } catch(const boost::program_options::error_with_option_name& e) {
        cout << "Invalid option: " << e.what() << endl;
        return 1;
    } catch (const runtime_error& e) {
        cout << timestamp << e.what() << endl;
        return 2;
    } catch (const std::bad_alloc& e) {
        cout << timestamp << e.what() << endl;
        cout << "Memory allocation failure." << endl;
        cout << "This assembly requires more memory than available." << endl;
        cout << "Rerun on a larger machine." << endl;
        return 2;
    } catch (const exception& e) {
        cout << timestamp << e.what() << endl;
        return 3;
    } catch (...) {
        cout << timestamp << "Terminated after catching a non-standard exception." << endl;
        return 4;
    }
    return 0;
}



void shasta::main::segmentFaultHandler(int)
{
    char message[] = "\nA segment fault occurred. Please report it by filing an "
        "issue on the Shasta repository and attaching the entire log output. "
        "To file an issue, point your browser to https://github.com/chanzuckerberg/shasta/issues\n";
    ::write(fileno(stderr), message, sizeof(message));
    ::_exit(1);
}


void shasta::main::main(int argumentCount, const char** arguments)
{
    struct sigaction action;
    ::memset(&action, 0, sizeof(action));
    action.sa_handler = &segmentFaultHandler;
    sigaction(SIGSEGV, &action, 0);

    // Parse command line options and the configuration file, if one was specified.
    AssemblerOptions assemblerOptions(argumentCount, arguments);
    cout << buildId() << endl;

    // Check that we have a valid command.
    auto it = commands.find(assemblerOptions.commandLineOnlyOptions.command);
    if(it ==commands.end()) {
        const string message = "Invalid command " + assemblerOptions.commandLineOnlyOptions.command;
        listCommands();
        throw runtime_error(message);
    }



    // Execute the requested command.
    if(assemblerOptions.commandLineOnlyOptions.command == "assemble") {
        assemble(assemblerOptions, argumentCount, arguments);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "cleanupBinaryData") {
        cleanupBinaryData(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "saveBinaryData") {
        saveBinaryData(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "explore") {
        explore(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "createBashCompletionScript") {
        createBashCompletionScript(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "listCommands") {
        listCommands();
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "listConfigurations") {
        listConfigurations();
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "listConfiguration") {
        listConfiguration(assemblerOptions);
        return;
    }

    // We already checked for a valid command above, so if we get here
    // the above logic is missing code for one of the valid commands.
    SHASTA_ASSERT(0);

}



// Implementation of --command assemble.
void shasta::main::assemble(
    const AssemblerOptions& assemblerOptions,
    int argumentCount, const char** arguments)
{
    SHASTA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "assemble");


    // Various checks for option validity.

    if(assemblerOptions.commandLineOnlyOptions.configName.empty()) {
        cout <<
            "Option \"--config\" is missing and is now required to "
            "run an assembly.\n"
            "It must specify either a configuration file\n"
            "or one of the following built-in configurations:\n";
        for(const auto& p: configurationTable) {
            cout << p.first << endl;
        }
        throw runtime_error(
            "Option \"--config\" is missing "
            "and is now required to run an assembly.");
    }

    // Check that we have at least one input file.
    if(assemblerOptions.commandLineOnlyOptions.inputFileNames.empty()) {
        cout << assemblerOptions.allOptionsDescription << endl;
        throw runtime_error("Specify at least one input file "
            "using command line option \"--input\".");
    }

    // Check assemblerOptions.minHashOptions.version.
    if( assemblerOptions.minHashOptions.version!=0 and
        assemblerOptions.minHashOptions.version!=1) {
        throw runtime_error("Invalid value " +
            to_string(assemblerOptions.minHashOptions.version) +
            " specified for --MinHash.version. Must be 0 or 1.");
    }

    // Check assemblerOptions.minHashOptions minimum/maximum bucket size.
    if( assemblerOptions.minHashOptions.maxBucketSize <=
        assemblerOptions.minHashOptions.minBucketSize) {
        throw runtime_error("MinHash maximum bucket size must be greater than minimum bucket size. "
            "The following values were specified:"
            " minimum bucket size " +
            to_string(assemblerOptions.minHashOptions.minBucketSize) +
            ", maximum bucket size " +
            to_string(assemblerOptions.minHashOptions.maxBucketSize) + "."
            );
    }

    // If coverage data was requested, memoryMode should be filesystem,
    // otherwise the coverage data cannot be accessed.
    if(assemblerOptions.assemblyOptions.storeCoverageData) {
        if(assemblerOptions.commandLineOnlyOptions.memoryMode != "filesystem") {
            throw runtime_error("To obtain usable coverage data, "
                "you must use --memoryMode filesystem.");
        }
    }

    if( assemblerOptions.alignOptions.alignMethod <  0 or
        assemblerOptions.alignOptions.alignMethod == 2 or
        assemblerOptions.alignOptions.alignMethod >  4) {
        throw runtime_error("Align method " + to_string(assemblerOptions.alignOptions.alignMethod) +
            " is not valid. Valid options are 0, 1, 3, and 4.");
    }

    if(assemblerOptions.readGraphOptions.creationMethod != 0 and
        assemblerOptions.readGraphOptions.creationMethod != 2) {
        throw runtime_error("--ReadGraph.creationMethod " +
            to_string(assemblerOptions.readGraphOptions.creationMethod) +
            " is not valid. Valid values are 0 and 2.");
    }

    // Check assemblerOptions.assemblyOptions.detangleMethod.
    if( assemblerOptions.assemblyOptions.detangleMethod!=0 and
        assemblerOptions.assemblyOptions.detangleMethod!=1 and
        assemblerOptions.assemblyOptions.detangleMethod!=2) {
        throw runtime_error("Invalid value " +
            to_string(assemblerOptions.assemblyOptions.detangleMethod) +
            " specified for --AssemblyOptions.detangleMethod. Must be 0, 1, or 2.");
    }

    // Check readGraphOptions.strandSeparationMethod.
    if(assemblerOptions.assemblyOptions.mode == 2 and
        assemblerOptions.readGraphOptions.strandSeparationMethod != 2) {
        throw runtime_error("--Assembly.mode 2 requires --ReadGraph.strandSeparationMethod 2.");
    }

    // Find absolute paths of the input files.
    // We will use them below after changing directory to the output directory.
    vector<string> inputFileAbsolutePaths;
    for(const string& inputFileName: assemblerOptions.commandLineOnlyOptions.inputFileNames) {
        if(!std::filesystem::exists(inputFileName)) {
            throw runtime_error("Input file not found: " + inputFileName);
        }
        if(!std::filesystem::is_regular_file(inputFileName)) {
            throw runtime_error("Input file is not a regular file: " + inputFileName);
        }
        inputFileAbsolutePaths.push_back(filesystem::getAbsolutePath(inputFileName));
    }



    // Create the assembly directory. If it exists and is not empty then stop.
    bool exists = std::filesystem::exists(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);
    bool isDir = std::filesystem::is_directory(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);
    if (exists) {
        if (!isDir) {
            throw runtime_error(
                assemblerOptions.commandLineOnlyOptions.assemblyDirectory +
                " already exists and is not a directory.\n"
                "Use --assemblyDirectory to specify a different assembly directory."
            );
        }
        bool isEmpty = filesystem::directoryContents(assemblerOptions.commandLineOnlyOptions.assemblyDirectory).empty();
        if (!isEmpty) {
            throw runtime_error(
                "Assembly directory " +
                assemblerOptions.commandLineOnlyOptions.assemblyDirectory +
                " exists and is not empty.\n"
                "Empty it for reuse or use --assemblyDirectory to specify a different assembly directory.");
        }
    } else {
        SHASTA_ASSERT(std::filesystem::create_directory(assemblerOptions.commandLineOnlyOptions.assemblyDirectory));
    }

    // Make the assembly directory current.
    std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);

    // Open the performance log.
    openPerformanceLog("performance.log");
    performanceLog << timestamp << "Assembly begins." << endl;

    // Open stdout.log and "tee" (duplicate) stdout to it.
    if(not assemblerOptions.commandLineOnlyOptions.suppressStdoutLog) {
        shastaLog.open("stdout.log");
        tee.duplicate(cout, shastaLog);
    }

    // Echo out the command line options.
    cout << timestamp << "Assembly begins.\nCommand line:" << endl;
    for(int i=0; i<argumentCount; i++) {
        cout << arguments[i] << " ";
    }
    cout << endl;



    // Set up the run directory as required by the memoryMode and memoryBacking options.
    size_t pageSize = 0;
    string dataDirectory;
    setupRunDirectory(
        assemblerOptions.commandLineOnlyOptions.memoryMode,
        assemblerOptions.commandLineOnlyOptions.memoryBacking,
        pageSize,
        dataDirectory);



    // Write out the option in effect to shasta.conf.
    {
        ofstream configurationFile("shasta.conf");
        assemblerOptions.write(configurationFile);
    }
    cout << "For options in use for this assembly, see shasta.conf in the assembly directory." << endl;



    // Initial disclaimer message.
    if(assemblerOptions.commandLineOnlyOptions.memoryBacking != "2M" &&
        assemblerOptions.commandLineOnlyOptions.memoryMode != "filesystem") {
        cout << "This run uses options \"--memoryBacking " << assemblerOptions.commandLineOnlyOptions.memoryBacking <<
            " --memoryMode " << assemblerOptions.commandLineOnlyOptions.memoryMode << "\".\n"
            "This could result in performance degradation.\n"
            "For full performance, use \"--memoryBacking 2M --memoryMode filesystem\"\n"
            "(root privilege via sudo required).\n"
            "Therefore the results of this run should not be used\n"
            "for benchmarking purposes." << endl;
    }

    // Create the Assembler.
    Assembler assembler(dataDirectory, true, assemblerOptions.readsOptions.representation, pageSize);
    assembler.assemblerInfo->readGraphCreationMethod = assemblerOptions.readGraphOptions.creationMethod;
    assembler.assemblerInfo->assemblyMode = assemblerOptions.assemblyOptions.mode;


    // Run the assembly.
    assemble(assembler, assemblerOptions, inputFileAbsolutePaths);

    // Final disclaimer message.
    if(assemblerOptions.commandLineOnlyOptions.memoryBacking != "2M" &&
        assemblerOptions.commandLineOnlyOptions.memoryMode != "filesystem") {
        cout << "This run used options \"--memoryBacking " << assemblerOptions.commandLineOnlyOptions.memoryBacking <<
            " --memoryMode " << assemblerOptions.commandLineOnlyOptions.memoryMode << "\".\n"
            "This could have resulted in performance degradation.\n"
            "For full performance, use \"--memoryBacking 2M --memoryMode filesystem\"\n"
            "(root privilege via sudo required).\n"
            "Therefore the results of this run should not be used\n"
            "for benchmarking purposes." << endl;
    }

    // Write out the build id again.
    cout << buildId() << endl;

    performanceLog << timestamp << "Assembly ends." << endl;
    cout << timestamp << "Assembly ends." << endl;
}



// Set up the run directory as required by the memoryMode and memoryBacking options.
void shasta::main::setupRunDirectory(
    const string& memoryMode,
    const string& memoryBacking,
    size_t& pageSize,
    string& dataDirectory
    )
{

    if(memoryMode == "anonymous") {

        if(memoryBacking == "disk") {

            // This combination is meaningless.
            throw runtime_error("\"--memoryMode anonymous\" is not allowed in combination "
                "with \"--memoryBacking disk\".");

        } else if(memoryBacking == "4K") {

            // Anonymous memory on 4KB pages.
            // This combination is the default.
            // It does not require root privilege.
            dataDirectory = "";
            pageSize = 4096;

        } else if(memoryBacking == "2M") {

            // Anonymous memory on 2MB pages.
            // This may require root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            // Root privilege is not required if 2M pages have already
            // been set up as required.
            setupHugePages();
            pageSize = 2 * 1024 * 1024;

        } else {
            throw runtime_error("Invalid value specified for --memoryBacking: " + memoryBacking +
                "\nValid values are: disk, 4K, 2M.");
        }

    } else if(memoryMode == "filesystem") {

        if(memoryBacking == "disk") {

            // Binary files on disk.
            // This does not require root privilege.
            SHASTA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 4096;

        } else if(memoryBacking == "4K") {

            // Binary files on the tmpfs filesystem
            // (filesystem in memory backed by 4K pages).
            // This requires root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            SHASTA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 4096;
            const string command = "sudo mount -t tmpfs -o size=0 tmpfs Data";
            const int errorCode = ::system(command.c_str());
            if(errorCode != 0) {
                throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
                    " running command: " + command);
            }

        } else if(memoryBacking == "2M") {

            // Binary files on the hugetlbfs filesystem
            // (filesystem in memory backed by 2M pages).
            // This requires root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            setupHugePages();
            SHASTA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 2 * 1024 * 1024;
            const uid_t userId = ::getuid();
            const gid_t groupId = ::getgid();
            const string command = "sudo mount -t hugetlbfs -o pagesize=2M"
                ",uid=" + to_string(userId) +
                ",gid=" + to_string(groupId) +
                " none Data";
            const int errorCode = ::system(command.c_str());
            if(errorCode != 0) {
                throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
                    " running command: " + command);
            }

        } else {
            throw runtime_error("Invalid value specified for --memoryBacking: " + memoryBacking +
                "\nValid values are: disk, 4K, 2M.");
        }

    } else {
        throw runtime_error("Invalid value specified for --memoryMode: " + memoryMode +
            "\nValid values are: anonymous, filesystem.");
    }
}



// This runs the entire assembly, under the following assumptions:
// - The current directory is the run directory.
// - The Data directory has already been created and set up, if necessary.
// - The input file names are either absolute,
//   or relative to the run directory, which is the current directory.
void shasta::main::assemble(
    Assembler& assembler,
    const AssemblerOptions& assemblerOptions,
    vector<string> inputFileNames)
{
    const auto steadyClock0 = std::chrono::steady_clock::now();
    const auto userClock0 = boost::chrono::process_user_cpu_clock::now();
    const auto systemClock0 = boost::chrono::process_system_cpu_clock::now();

    // Adjust the number of threads, if necessary.
    uint32_t threadCount = assemblerOptions.commandLineOnlyOptions.threadCount;
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "This assembly will use " << threadCount << " threads." << endl;

    // Set up the consensus caller.
    cout << "Setting up consensus caller " <<
        assemblerOptions.assemblyOptions.consensusCaller << endl;
    assembler.setupConsensusCaller(assemblerOptions.assemblyOptions.consensusCaller);



    // Add reads from the specified input files.
    performanceLog << timestamp << "Begin loading reads from " << inputFileNames.size() << " files." << endl;
    const auto t0 = steady_clock::now();
    for(const string& inputFileName: inputFileNames) {

        assembler.addReads(
            inputFileName,
            assemblerOptions.readsOptions.minReadLength,
            assemblerOptions.readsOptions.noCache,
            threadCount);
    }

    if(assembler.getReads().readCount() == 0) {
        throw runtime_error("There are no input reads.");
    }



    // If requested, increase the read length cutoff
    // to reduce coverage to the specified amount.
    if (assemblerOptions.readsOptions.desiredCoverage > 0) {
        // Write out the read length histogram using provided minReadLength.
        assembler.histogramReadLength("ExtendedReadLengthHistogram.csv");

        const auto newMinReadLength = assembler.adjustCoverageAndGetNewMinReadLength(
            assemblerOptions.readsOptions.desiredCoverage);

        const auto oldMinReadLength = uint64_t(assemblerOptions.readsOptions.minReadLength);

        if (newMinReadLength == 0ULL) {
            throw runtime_error(
                "With Reads.minReadLength " +
                to_string(assemblerOptions.readsOptions.minReadLength) +
                ", total available coverage is " +
                to_string(assembler.getReads().getTotalBaseCount()) +
                ", less than desired coverage " +
                to_string(assemblerOptions.readsOptions.desiredCoverage) +
                ". Try reducing Reads.minReadLength if appropriate or get more coverage."
            );
        }

        // Adjusting coverage should only ever reduce coverage if necessary.
        SHASTA_ASSERT(newMinReadLength >= oldMinReadLength);
    }

    assembler.computeReadIdsSortedByName();
    assembler.histogramReadLength("ReadLengthHistogram.csv");

    const auto t1 = steady_clock::now();
    performanceLog << timestamp << "Done loading reads from " << inputFileNames.size() << " files." << endl;
    performanceLog << "Read loading took " << seconds(t1-t0) << "s." << endl;



    // Select the k-mers that will be used as markers.
    switch(assemblerOptions.kmersOptions.generationMethod) {
    case 0:
        assembler.randomlySelectKmers(
            assemblerOptions.kmersOptions.k,
            assemblerOptions.kmersOptions.probability, 231);
        break;

    case 1:
        // Randomly select the k-mers to be used as markers, but
        // excluding those that are globally overenriched in the input reads,
        // as measured by total frequency in all reads.
        assembler.selectKmersBasedOnFrequency(
            assemblerOptions.kmersOptions.k,
            assemblerOptions.kmersOptions.probability, 231,
            assemblerOptions.kmersOptions.enrichmentThreshold, threadCount);
        break;

    case 2:
        // Randomly select the k-mers to be used as markers, but
        // excluding those that are overenriched even in a single oriented read.
        assembler.selectKmers2(
            assemblerOptions.kmersOptions.k,
            assemblerOptions.kmersOptions.probability, 231,
            assemblerOptions.kmersOptions.enrichmentThreshold, threadCount);
        break;

    case 3:
        // Read the k-mers to be used as markers from a file.
        if(assemblerOptions.kmersOptions.file.empty() or
            assemblerOptions.kmersOptions.file[0] != '/') {
            throw runtime_error("Option --Kmers.file must specify an absolute path. "
                "A relative path is not accepted.");
        }
        assembler.readKmersFromFile(
            assemblerOptions.kmersOptions.k,
            assemblerOptions.kmersOptions.file);
        break;

    case 4:
        // Randomly select the k-mers to be used as markers, but
        // excluding those that appear in two copies close to each other
        // even in a single oriented read.
        assembler.selectKmers4(
            assemblerOptions.kmersOptions.k,
            assemblerOptions.kmersOptions.probability, 231,
            assemblerOptions.kmersOptions.distanceThreshold, threadCount);
        break;

    default:
        throw runtime_error("Invalid --Kmers generationMethod. "
            "Specify a value between 0 and 4, inclusive.");
    }



    // Find the markers in the reads.
    assembler.findMarkers(0);

    if(!assemblerOptions.readsOptions.palindromicReads.skipFlagging) {

        // Flag palindromic reads.
        // These will be excluded from further processing.
        assembler.flagPalindromicReads(
            assemblerOptions.readsOptions.palindromicReads.maxSkip,
            assemblerOptions.readsOptions.palindromicReads.maxDrift,
            assemblerOptions.readsOptions.palindromicReads.maxMarkerFrequency,
            assemblerOptions.readsOptions.palindromicReads.alignedFractionThreshold,
            assemblerOptions.readsOptions.palindromicReads.nearDiagonalFractionThreshold,
            assemblerOptions.readsOptions.palindromicReads.deltaThreshold,
            threadCount);
    }



    // Find alignment candidates.
    if(assemblerOptions.minHashOptions.allPairs) {
        assembler.markAlignmentCandidatesAllPairs();
    } else if(assemblerOptions.minHashOptions.version == 0) {
        assembler.findAlignmentCandidatesLowHash0(
            assemblerOptions.minHashOptions.m,
            assemblerOptions.minHashOptions.hashFraction,
            assemblerOptions.minHashOptions.minHashIterationCount,
            assemblerOptions.minHashOptions.alignmentCandidatesPerRead,
            0,
            assemblerOptions.minHashOptions.minBucketSize,
            assemblerOptions.minHashOptions.maxBucketSize,
            assemblerOptions.minHashOptions.minFrequency,
            threadCount);
    } else {
        SHASTA_ASSERT(assemblerOptions.minHashOptions.version == 1);    // Already checked for that.
        assembler.findAlignmentCandidatesLowHash1(
            assemblerOptions.minHashOptions.m,
            assemblerOptions.minHashOptions.hashFraction,
            assemblerOptions.minHashOptions.minHashIterationCount,
            0,
            assemblerOptions.minHashOptions.minBucketSize,
            assemblerOptions.minHashOptions.maxBucketSize,
            assemblerOptions.minHashOptions.minFrequency,
            threadCount);
    }



    // Suppress alignment candidates where reads are close on the same channel.
    if(assemblerOptions.alignOptions.sameChannelReadAlignmentSuppressDeltaThreshold > 0) {
        assembler.suppressAlignmentCandidates(
            assemblerOptions.alignOptions.sameChannelReadAlignmentSuppressDeltaThreshold,
            threadCount);
    }


    // For http server and debugging/development purposes, generate an exhaustive table of candidates
    assembler.computeCandidateTable();


    // Compute alignments.
    assembler.computeAlignments(
        assemblerOptions.alignOptions,
        threadCount);



    // Create the read graph.
    if(assemblerOptions.readGraphOptions.creationMethod == 0) {
        assembler.createReadGraph(
            assemblerOptions.readGraphOptions.maxAlignmentCount,
            assemblerOptions.alignOptions.maxTrim);

        // Actual alignment criteria are as specified in the command line options
        // and/or configuration.
        assembler.assemblerInfo->actualMinAlignedFraction = assemblerOptions.alignOptions.minAlignedFraction;
        assembler.assemblerInfo->actualMinAlignedMarkerCount = assemblerOptions.alignOptions.minAlignedMarkerCount;
        assembler.assemblerInfo->actualMaxDrift = assemblerOptions.alignOptions.maxDrift;
        assembler.assemblerInfo->actualMaxSkip = assemblerOptions.alignOptions.maxSkip;
        assembler.assemblerInfo->actualMaxTrim = assemblerOptions.alignOptions.maxTrim;


    } else if(assemblerOptions.readGraphOptions.creationMethod == 2) {
        assembler.createReadGraph2(
            assemblerOptions.readGraphOptions.maxAlignmentCount,
            assemblerOptions.readGraphOptions.markerCountPercentile,
            assemblerOptions.readGraphOptions.alignedFractionPercentile,
            assemblerOptions.readGraphOptions.maxSkipPercentile,
            assemblerOptions.readGraphOptions.maxDriftPercentile,
            assemblerOptions.readGraphOptions.maxTrimPercentile);
    } else {
        throw runtime_error("Invalid value for --ReadGraph.creationMethod.");
    }

    // Limited strand separation.
    // If strict strand separation is requested, it is done later,
    // after chimera detection.
    if(assemblerOptions.readGraphOptions.strandSeparationMethod == 1) {
        assembler.flagCrossStrandReadGraphEdges1(
            assemblerOptions.readGraphOptions.crossStrandMaxDistance,
            threadCount);
    }

    // Flag chimeric reads.
    assembler.flagChimericReads(assemblerOptions.readGraphOptions.maxChimericReadDistance, threadCount);

    // Flag inconsistent alignments, if requested.
    if(assemblerOptions.readGraphOptions.flagInconsistentAlignments) {
        assembler.flagInconsistentAlignments(
            assemblerOptions.readGraphOptions.flagInconsistentAlignmentsTriangleErrorThreshold,
            assemblerOptions.readGraphOptions.flagInconsistentAlignmentsLeastSquareErrorThreshold,
            assemblerOptions.readGraphOptions.flagInconsistentAlignmentsLeastSquareMaxDistance,
            threadCount);
    }

    // Strict strand separation.
    if(assemblerOptions.readGraphOptions.strandSeparationMethod == 2) {
        assembler.flagCrossStrandReadGraphEdges2();
    }

    // Compute connected components of the read graph.
    // These are currently not used.
    // For strand separation method 2 this was already done
    // in flagCrossStrandReadGraphEdges2.
    if(assemblerOptions.readGraphOptions.strandSeparationMethod != 2) {
        assembler.computeReadGraphConnectedComponents();
    }



    // Do the rest of the assembly using the selected assembly mode.
    switch(assemblerOptions.assemblyOptions.mode) {
    case 0:
        mode0Assembly(assembler, assemblerOptions, threadCount);
        break;
    case 2:
        mode2Assembly(assembler, assemblerOptions, threadCount);
        break;
    case 3:
        mode3Assembly(assembler, assemblerOptions, threadCount);
        break;
    default:
        throw runtime_error("Invalid value specified for --Assembly.mode. "
            "Valid values are 0 (haploid assembly) and 2 (phased diploid assembly), but " +
            to_string(assemblerOptions.assemblyOptions.mode) +
            " was specified.");
    }


    // Store elapsed time for assembly.
    const auto steadyClock1 = std::chrono::steady_clock::now();
    const auto userClock1 = boost::chrono::process_user_cpu_clock::now();
    const auto systemClock1 = boost::chrono::process_system_cpu_clock::now();
    const double elapsedTime = 1.e-9 * double((
        std::chrono::duration_cast<std::chrono::nanoseconds>(steadyClock1 - steadyClock0)).count());
    const double userTime = 1.e-9 * double((
        boost::chrono::duration_cast<boost::chrono::nanoseconds>(userClock1 - userClock0)).count());
    const double systemTime = 1.e-9 * double((
        boost::chrono::duration_cast<boost::chrono::nanoseconds>(systemClock1 - systemClock0)).count());
    const double averageCpuUtilization =
        (userTime + systemTime) / (double(std::thread::hardware_concurrency()) * elapsedTime);
    assembler.storeAssemblyTime(elapsedTime, averageCpuUtilization);

    // Store peak memory usage.
    uint64_t peakMemoryUsage = getPeakMemoryUsage();
    assembler.storePeakMemoryUsage(peakMemoryUsage);

    // Store other performance information.
    assembler.assemblerInfo->threadCount = threadCount;
    assembler.assemblerInfo->virtualCpuCount = std::thread::hardware_concurrency();
    assembler.assemblerInfo->totalAvailableMemory = getTotalPhysicalMemory();

    // Write a summary of read information.
    assembler.writeReadsSummary();

    // Write the assembly summary.
    ofstream html("AssemblySummary.html");
    assembler.writeAssemblySummary(html);
    ofstream json("AssemblySummary.json");
    assembler.writeAssemblySummaryJson(json);
    ofstream htmlIndex("index.html");
    assembler.writeAssemblyIndex(htmlIndex);

    performanceLog << timestamp << endl;
    performanceLog << "Assembly time statistics:\n"
        "    Elapsed seconds: " << elapsedTime << "\n"
        "    Elapsed minutes: " << elapsedTime/60. << "\n"
        "    Elapsed hours:   " << elapsedTime/3600. << "\n";
    performanceLog << "Average CPU utilization: " << averageCpuUtilization << endl;
    performanceLog << "Peak Memory usage: " << peakMemoryUsage << " bytes = " <<
        int(std::round(double(peakMemoryUsage) / (1024. * 1024. * 1024.)) ) << " GiB" << endl;

}



void shasta::main::mode0Assembly(
    Assembler& assembler,
    const AssemblerOptions& assemblerOptions,
    uint32_t threadCount)
{

    // Iterative assembly, if requested (experimental).
    if(assemblerOptions.assemblyOptions.iterative) {
        for(uint64_t iteration=0;
            iteration<assemblerOptions.assemblyOptions.iterativeIterationCount;
            iteration++) {
            cout << timestamp << "Iterative assembly iteration " << iteration << " begins." << endl;

            // Do an assembly with the current read graph, without marker graph
            // simplification or detangling.
            assembler.createMarkerGraphVertices(
                assemblerOptions.markerGraphOptions.minCoverage,
                assemblerOptions.markerGraphOptions.maxCoverage,
                assemblerOptions.markerGraphOptions.minCoveragePerStrand,
                assemblerOptions.markerGraphOptions.allowDuplicateMarkers,
                assemblerOptions.markerGraphOptions.peakFinderMinAreaFraction,
                assemblerOptions.markerGraphOptions.peakFinderAreaStartIndex,
                threadCount);
            assembler.findMarkerGraphReverseComplementVertices(threadCount);
            assembler.createMarkerGraphEdges(threadCount);
            assembler.findMarkerGraphReverseComplementEdges(threadCount);
            assembler.transitiveReduction(
                assemblerOptions.markerGraphOptions.lowCoverageThreshold,
                assemblerOptions.markerGraphOptions.highCoverageThreshold,
                assemblerOptions.markerGraphOptions.maxDistance,
                assemblerOptions.markerGraphOptions.edgeMarkerSkipThreshold);
            assembler.pruneMarkerGraphStrongSubgraph(
                assemblerOptions.markerGraphOptions.pruneIterationCount);
            assembler.createAssemblyGraphEdges();
            assembler.createAssemblyGraphVertices();

            // Recreate the read graph using pseudo-paths from this assembly.
            assembler.createReadGraphUsingPseudoPaths(
                assemblerOptions.assemblyOptions.iterativePseudoPathAlignMatchScore,
                assemblerOptions.assemblyOptions.iterativePseudoPathAlignMismatchScore,
                assemblerOptions.assemblyOptions.iterativePseudoPathAlignGapScore,
                assemblerOptions.assemblyOptions.iterativeMismatchSquareFactor,
                assemblerOptions.assemblyOptions.iterativeMinScore,
                assemblerOptions.assemblyOptions.iterativeMaxAlignmentCount,
                threadCount);
            for(uint64_t bridgeRemovalIteration=0;
                bridgeRemovalIteration<assemblerOptions.assemblyOptions.iterativeBridgeRemovalIterationCount;
                bridgeRemovalIteration++) {
                assembler.removeReadGraphBridges(
                    assemblerOptions.assemblyOptions.iterativeBridgeRemovalMaxDistance);
            }

            // Remove the marker graph and assembly graph we created in the process.
            assembler.markerGraph.remove();
            assembler.assemblyGraphPointer.reset();

        }

        // Now we have a new read graph with some amount of separation
        // between copies of long repeats and/or haplotypes.
        // The rest of the assembly continues normally.
    }



    // Create marker graph vertices.
    // This uses a disjoint sets data structure to merge markers
    // that are aligned based on an alignment present in the read graph.
    assembler.createMarkerGraphVertices(
        assemblerOptions.markerGraphOptions.minCoverage,
        assemblerOptions.markerGraphOptions.maxCoverage,
        assemblerOptions.markerGraphOptions.minCoveragePerStrand,
        assemblerOptions.markerGraphOptions.allowDuplicateMarkers,
        assemblerOptions.markerGraphOptions.peakFinderMinAreaFraction,
        assemblerOptions.markerGraphOptions.peakFinderAreaStartIndex,
        threadCount);

    // Find the reverse complement of each marker graph vertex.
    assembler.findMarkerGraphReverseComplementVertices(threadCount);

    // Clean up of duplicate markers, if requested and necessary.
    if(assemblerOptions.markerGraphOptions.allowDuplicateMarkers and
        assemblerOptions.markerGraphOptions.cleanupDuplicateMarkers) {
        assembler.cleanupDuplicateMarkers(
            threadCount,
            assembler.getMarkerGraphMinCoverageUsed(),    // Stored by createMarkerGraphVertices.
            assemblerOptions.markerGraphOptions.minCoveragePerStrand,
            assemblerOptions.markerGraphOptions.duplicateMarkersPattern1Threshold,
            false, false);
    }

    // Create edges of the marker graph.
    assembler.createMarkerGraphEdges(threadCount);
    assembler.findMarkerGraphReverseComplementEdges(threadCount);

    // Approximate transitive reduction.
    assembler.transitiveReduction(
        assemblerOptions.markerGraphOptions.lowCoverageThreshold,
        assemblerOptions.markerGraphOptions.highCoverageThreshold,
        assemblerOptions.markerGraphOptions.maxDistance,
        assemblerOptions.markerGraphOptions.edgeMarkerSkipThreshold);
    if(assemblerOptions.markerGraphOptions.reverseTransitiveReduction) {
        assembler.reverseTransitiveReduction(
            assemblerOptions.markerGraphOptions.lowCoverageThreshold,
            assemblerOptions.markerGraphOptions.highCoverageThreshold,
            assemblerOptions.markerGraphOptions.maxDistance);
    }



    // Prune the marker graph.
    assembler.pruneMarkerGraphStrongSubgraph(
        assemblerOptions.markerGraphOptions.pruneIterationCount);

    // Compute marker graph coverage histogram.
    assembler.computeMarkerGraphCoverageHistogram();

    // Simplify the marker graph to remove bubbles and superbubbles.
    // The maxLength parameter controls the maximum number of markers
    // for a branch to be collapsed during each iteration.
    assembler.simplifyMarkerGraph(assemblerOptions.markerGraphOptions.simplifyMaxLengthVector, false);

    // Create the assembly graph.
    assembler.createAssemblyGraphEdges();
    assembler.createAssemblyGraphVertices();

    // Remove low-coverage cross-edges from the assembly graph and
    // the corresponding marker graph edges.
    if(assemblerOptions.markerGraphOptions.crossEdgeCoverageThreshold > 0.) {
        assembler.removeLowCoverageCrossEdges(
            uint32_t(assemblerOptions.markerGraphOptions.crossEdgeCoverageThreshold));
        assembler.assemblyGraphPointer->remove();
        assembler.createAssemblyGraphEdges();
        assembler.createAssemblyGraphVertices();
    }

    // Prune the assembly graph, if requested.
    if(assemblerOptions.assemblyOptions.pruneLength > 0) {
        assembler.pruneAssemblyGraph(assemblerOptions.assemblyOptions.pruneLength);
    }

    // Detangle, if requested.
    if(assemblerOptions.assemblyOptions.detangleMethod == 1) {
        assembler.detangle();
    } else if(assemblerOptions.assemblyOptions.detangleMethod == 2) {
        assembler.detangle2(
            assemblerOptions.assemblyOptions.detangleDiagonalReadCountMin,
            assemblerOptions.assemblyOptions.detangleOffDiagonalReadCountMax,
            assemblerOptions.assemblyOptions.detangleOffDiagonalRatio
            );
    }

    // If any detangling was done, remove low-coverage cross-edges again.
    if(assemblerOptions.assemblyOptions.detangleMethod != 0 and
        assemblerOptions.markerGraphOptions.crossEdgeCoverageThreshold > 0.) {
        assembler.removeLowCoverageCrossEdges(
            uint32_t(assemblerOptions.markerGraphOptions.crossEdgeCoverageThreshold));
        assembler.assemblyGraphPointer->remove();
        assembler.createAssemblyGraphEdges();
        assembler.createAssemblyGraphVertices();
    }

    // Compute optimal repeat counts for each vertex of the marker graph.
    if(assemblerOptions.readsOptions.representation == 1) {
        assembler.assembleMarkerGraphVertices(threadCount);
    }

    // If coverage data was requested, compute and store coverage data for the vertices.
    if(assemblerOptions.assemblyOptions.storeCoverageData or
        assemblerOptions.assemblyOptions.storeCoverageDataCsvLengthThreshold>0) {
        assembler.computeMarkerGraphVerticesCoverageData(threadCount);
    }

    // Compute consensus sequence for marker graph edges to be used for assembly.
    assembler.assembleMarkerGraphEdges(
        threadCount,
        assemblerOptions.assemblyOptions.markerGraphEdgeLengthThresholdForConsensus,
        assemblerOptions.assemblyOptions.storeCoverageData or
        assemblerOptions.assemblyOptions.storeCoverageDataCsvLengthThreshold>0,
        false
        );

    // Use the assembly graph for global assembly.
    assembler.assemble(
        threadCount,
        assemblerOptions.assemblyOptions.storeCoverageDataCsvLengthThreshold);
    // assembler.findAssemblyGraphBubbles();
    assembler.computeAssemblyStatistics();
    assembler.writeGfa1("Assembly.gfa");
    assembler.writeGfa1BothStrands("Assembly-BothStrands.gfa");
    assembler.writeGfa1BothStrandsNoSequence("Assembly-BothStrands-NoSequence.gfa");
    assembler.writeFasta("Assembly.fasta");

    // If requested, write out the oriented reads that were used to assemble
    // each assembled segment.
    if(assemblerOptions.assemblyOptions.writeReadsByAssembledSegment) {
        cout << timestamp << " Writing the oriented reads that were used to assemble each segment." << endl;
        assembler.gatherOrientedReadsByAssemblyGraphEdge(threadCount);
        assembler.writeOrientedReadsByAssemblyGraphEdge();
    }
}



void shasta::main::mode2Assembly(
    Assembler& assembler,
    const AssemblerOptions& assemblerOptions,
    uint32_t threadCount)
{
    // Create marker graph vertices.
    assembler.createMarkerGraphVertices(
        assemblerOptions.markerGraphOptions.minCoverage,
        assemblerOptions.markerGraphOptions.maxCoverage,
        assemblerOptions.markerGraphOptions.minCoveragePerStrand,
        assemblerOptions.markerGraphOptions.allowDuplicateMarkers,
        assemblerOptions.markerGraphOptions.peakFinderMinAreaFraction,
        assemblerOptions.markerGraphOptions.peakFinderAreaStartIndex,
        threadCount);
    assembler.findMarkerGraphReverseComplementVertices(threadCount);

    // Create marker graph edges.
    // For assembly mode 1 we use createMarkerGraphEdgesStrict
    // with minimum edge coverage (total and per strand).
    assembler.createMarkerGraphEdgesStrict(
        assemblerOptions.markerGraphOptions.minEdgeCoverage,
        assemblerOptions.markerGraphOptions.minEdgeCoveragePerStrand, threadCount);
    assembler.findMarkerGraphReverseComplementEdges(threadCount);

    // Coverage histograms for vertices and edges of the marker graph.
    assembler.computeMarkerGraphCoverageHistogram();

    // To recover contiguity, add secondary edges.
    assembler.createMarkerGraphSecondaryEdges(
        uint32_t(assemblerOptions.markerGraphOptions.secondaryEdgesMaxSkip),
        threadCount);
    assembler.splitMarkerGraphSecondaryEdges(
        assemblerOptions.markerGraphOptions.secondaryEdgesSplitErrorRateThreshold,
        assemblerOptions.markerGraphOptions.secondaryEdgesSplitMinCoverage,
        threadCount);

    // Coverage histograms for vertices and edges of the marker graph.
    assembler.computeMarkerGraphCoverageHistogram();

    // Compute optimal repeat counts for each vertex of the marker graph.
    if(assemblerOptions.readsOptions.representation == 1) {
        assembler.assembleMarkerGraphVertices(threadCount);
    }

    // Compute consensus sequence for all marker graph edges.
    assembler.assembleMarkerGraphEdges(
        threadCount,
        assemblerOptions.assemblyOptions.markerGraphEdgeLengthThresholdForConsensus,
        assemblerOptions.assemblyOptions.storeCoverageData or
        assemblerOptions.assemblyOptions.storeCoverageDataCsvLengthThreshold>0,
        true
        );

    // Create the mode 2 assembly graph.
    assembler.createAssemblyGraph2(
        assemblerOptions.assemblyOptions.pruneLength,
        assemblerOptions.assemblyOptions.mode2Options,
        threadCount, false);


}






void shasta::main::mode3Assembly(
    Assembler& assembler,
    const AssemblerOptions& assemblerOptions,
    uint32_t threadCount)
{
    // Create marker graph vertices.
    assembler.createMarkerGraphVertices(
        assemblerOptions.markerGraphOptions.minCoverage,
        assemblerOptions.markerGraphOptions.maxCoverage,
        assemblerOptions.markerGraphOptions.minCoveragePerStrand,
        assemblerOptions.markerGraphOptions.allowDuplicateMarkers,
        assemblerOptions.markerGraphOptions.peakFinderMinAreaFraction,
        assemblerOptions.markerGraphOptions.peakFinderAreaStartIndex,
        threadCount);
    assembler.findMarkerGraphReverseComplementVertices(threadCount);

    // Create marker graph edges.
    // For assembly mode 1 we use createMarkerGraphEdgesStrict
    // with minimum edge coverage (total and per strand).
    assembler.createMarkerGraphEdgesStrict(
        assemblerOptions.markerGraphOptions.minEdgeCoverage,
        assemblerOptions.markerGraphOptions.minEdgeCoveragePerStrand, threadCount);
    assembler.findMarkerGraphReverseComplementEdges(threadCount);

    // Coverage histograms for vertices and edges of the marker graph.
    assembler.computeMarkerGraphCoverageHistogram();

    // In mode 3 assembly, we don't add secondary edges.

    // Coverage histograms for vertices and edges of the marker graph.
    assembler.computeMarkerGraphCoverageHistogram();

    // Compute optimal repeat counts for each vertex of the marker graph.
    if(assemblerOptions.readsOptions.representation == 1) {
        assembler.assembleMarkerGraphVertices(threadCount);
    }

    // Compute consensus sequence for all marker graph edges.
    assembler.assembleMarkerGraphEdges(
        threadCount,
        assemblerOptions.assemblyOptions.markerGraphEdgeLengthThresholdForConsensus,
        assemblerOptions.assemblyOptions.storeCoverageData or
        assemblerOptions.assemblyOptions.storeCoverageDataCsvLengthThreshold>0,
        true
        );

    // Run mode 3 assembly.
    assembler.mode3Assembly(
        threadCount);


}



// This function sets nr_overcommit_hugepages for 2MB pages
// to a little below total memory.
// If the setting needs to be modified, it acquires
// root privilege via sudo. This may result in the
// user having to enter a password.
void shasta::main::setupHugePages()
{

    // Get the total memory size.
    const uint64_t totalMemoryBytes = sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES);

    // Figure out how much memory we want to allow for 2MB pages.
    const uint64_t MB = 1024 * 1024;
    const uint64_t GB = MB * 1024;
    const uint64_t maximumHugePageMemoryBytes = totalMemoryBytes - 8 * GB;
    const uint64_t maximumHugePageMemoryHugePages = maximumHugePageMemoryBytes / (2 * MB);

    // Check what we have it set to.
    const string fileName = "/sys/kernel/mm/hugepages/hugepages-2048kB/nr_overcommit_hugepages";
    ifstream file(fileName);
    if(!file) {
        throw runtime_error("Error opening " + fileName + " for read.");
    }
    uint64_t currentValue = 0;
    file >> currentValue;
    file.close();

    // If it's set to at least what we want, don't do anything.
    // When this happens, root access is not required.
    if(currentValue >= maximumHugePageMemoryHugePages) {
        return;
    }

    // Use sudo to set.
    const string command =
        "sudo sh -c \"echo " +
        to_string(maximumHugePageMemoryHugePages) +
        " > " + fileName + "\"";
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " running command: " + command);
    }

}



// Implementation of --command saveBinaryData.
// This copies Data to DataOnDisk.
void shasta::main::saveBinaryData(
    const AssemblerOptions& assemblerOptions)
{
    SHASTA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "saveBinaryData");

    // Locate the Data directory.
    const string dataDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/Data";
    if(!std::filesystem::exists(dataDirectory)) {
        throw runtime_error(dataDirectory + " does not exist, nothing done.");
    }

    // Check that the DataOnDisk directory does not exist.
    const string dataOnDiskDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/DataOnDisk";
    if(std::filesystem::exists(dataOnDiskDirectory)) {
        throw runtime_error(dataOnDiskDirectory + " already exists, nothing done.");
    }

    // Copy Data to DataOnDisk.
    const string command = "cp -rp " + dataDirectory + " " + dataOnDiskDirectory;
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " running command:\n" + command);
    }
    cout << "Binary data successfully saved." << endl;
}



// Implementation of --command cleanupBinaryData.
void shasta::main::cleanupBinaryData(
    const AssemblerOptions& assemblerOptions)
{
    SHASTA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "cleanupBinaryData");

    // Locate the Data directory.
    const string dataDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/Data";
    if(!std::filesystem::exists(dataDirectory)) {
        cout << dataDirectory << " does not exist, nothing done." << endl;
        return;
    }

    // Unmount it and remove it.
    ::system(("sudo umount " + dataDirectory).c_str());
    const int errorCode = ::system(string("rm -rf " + dataDirectory).c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " removing " + dataDirectory);
    }
    cout << "Cleanup of " << dataDirectory << " successful." << endl;

    // If the DataOnDisk directory exists, create a symbolic link
    // Data->DataOnDisk.
    const string dataOnDiskDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/DataOnDisk";
    if(std::filesystem::exists(dataOnDiskDirectory)) {
        std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);
        const string command = "ln -s DataOnDisk Data";
        ::system(command.c_str());
    }

}

// Implementation of --command explore.
void shasta::main::explore(
    const AssemblerOptions& assemblerOptions)
{
    // If a paf file was specified, find its absolute path
    // before we switch to the assembly directory.
    string alignmentsPafFileAbsolutePath;
    if(not assemblerOptions.commandLineOnlyOptions.alignmentsPafFile.empty()) {
        if(!std::filesystem::exists(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile)) {
            throw runtime_error(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile + " not found.");
        }
        if(!std::filesystem::is_regular_file(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile)) {
            throw runtime_error(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile + " is not a regular file.");
        }
        alignmentsPafFileAbsolutePath = filesystem::getAbsolutePath(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile);
    }

    // Go to the assembly directory.
    std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);

    // Check that we have the binary data.
    if(!std::filesystem::exists("Data")) {
        throw runtime_error("Binary directory \"Data\" not available "
        " in assembly directory " +
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory +
        ". Use \"--memoryMode filesystem\", possibly followed by "
        "\"--command saveBinaryData\" and \"--command cleanupBinaryData\" "
        "if you want to make sure the binary data are persistently available on disk. "
        "See the documentations are some of these options require root access."
        );
        return;
    }

    // Create the Assembler.
    Assembler assembler("Data/", false, 1, 0);

    // Set up the consensus caller.
    cout << "Setting up consensus caller " <<
        assemblerOptions.assemblyOptions.consensusCaller << endl;
    assembler.setupConsensusCaller(assemblerOptions.assemblyOptions.consensusCaller);

    // Access all available binary data.
    assembler.accessAllSoft();

    string executablePath = filesystem::executablePath();
    // On Linux it will be something like - `/path/to/install_root/bin/shasta`

    string executableBinPath = executablePath.substr(0, executablePath.find_last_of('/'));
    string installRootPath = executableBinPath.substr(0, executableBinPath.find_last_of('/'));
    string docsPath = installRootPath + "/docs";

    if (std::filesystem::is_directory(docsPath)) {
        assembler.httpServerData.docsDirectory = docsPath;
    } else {
        cout << "Documentation is not available." << endl;
        assembler.httpServerData.docsDirectory = "";
    }

    // Load the paf file, if one was specified.
    if(not alignmentsPafFileAbsolutePath.empty()) {
        assembler.loadAlignmentsPafFile(alignmentsPafFileAbsolutePath);
    }

    // Start the http server.
    assembler.httpServerData.assemblerOptions = &assemblerOptions;
    bool localOnly;
    bool sameUserOnly;
    if(assemblerOptions.commandLineOnlyOptions.exploreAccess == "user") {
        localOnly = true;
        sameUserOnly = true;
    } else if(assemblerOptions.commandLineOnlyOptions.exploreAccess == "local") {
        localOnly = true;
        sameUserOnly = false;
    } else if (assemblerOptions.commandLineOnlyOptions.exploreAccess == "unrestricted"){
        localOnly = false;
        sameUserOnly = false;
    } else {
        throw runtime_error("Invalid value specified for --exploreAccess. "
            "Only use this option if you understand its security implications."
        );
    }
    assembler.explore(
        assemblerOptions.commandLineOnlyOptions.port,
        localOnly,
        sameUserOnly);
}



// This creates a bash completion script for the Shasta executable,
// which makes it easier to type long option names.
// To use it:
// shasta --command createBashCompletionScript; source shastaCompletion.sh
// Then, press TAB once or twice while editing a Shasta command line
// to get the Bash shell to suggest or fill in possibilities.
// You can put the "source" command in your .bashrc or other
// appropriate location.
// THIS IS AN INITIAL CUT AND LACKS MANY DESIRABLE FEATURES,
// LIKE FOR EXAMPLE COMPLETION OF FILE NAMES (AFTER --input),
// AND THE ABILITY TO COMPLETE KEYWORDS ONLY AFTER THE OPTION THEY
// SHOULD BE PRECEDED BY.
// IF SOMEBODY WITH A GOOD UNDERSTANDING OF BASH COMPLETION SEES THIS,
// PLEASE MAKE IT BETTER AND SUBMIT A PULL REQUEST!
void shasta::main::createBashCompletionScript(const AssemblerOptions& assemblerOptions)
{
    const string fileName = "shastaCompletion.sh";
    ofstream file(fileName);

    file << "#!/bin/bash\n";
    file << "complete -o default -W \"\\\n";

    // Options.
    for(const auto& option: assemblerOptions.allOptionsDescription.options()) {
        file << "--" << option->long_name() << " \\\n";
    }

    // Commands.
    for(const auto& command: commands) {
        file << command << " \\\n";
    }

    // Built-in configurations.
    for(const auto& p: configurationTable) {
        file << p.first << " \\\n";
    }

    // Bayesian models.
    for(const string& name: SimpleBayesianConsensusCaller::builtIns) {
        file << name << " \\\n";
    }

    // Other keywords. This should be modified to only accept them after the appropriate option.
    file << "filesystem anonymous \\\n";
    file << "disk 4K 2M \\\n";
    file << "user local unrestricted \\\n";
    file << "Bayesian Modal Median \\\n";

    // Finish the "complete" command.
    file << "\" shasta\n";

    cout << "Created shastaCompletion.sh. "
        "In the bash shell, use the following command to "
        "get shell command completion when invoking Shasta:\n"
        "source shastaCompletion.sh\n"
        "This makes it easier to type when running Shasta." << endl;

}



void shasta::main::listCommands()
{
    cout << "Valid commands are:" << endl;
    for(const string& command: commands) {
        cout << command << endl;
    }
}



void shasta::main::listConfigurations()
{
    cout << "Valid Shasta built-in configurations, in chronological order, are:\n" << endl;
    for(const auto& p: configurationTable) {
        cout << p.first << endl;
    }
    cout <<
        "\nUse \"shasta --command listConfiguration --config configurationName\" "
        "to list the details of one of the above configurations.\n\n"
        "When running an assembly, you can use option \"--config\" "
        "to specify any of the above configuration names, "
        "or the name of a configuration file. "
        "See shasta/conf for examples of configuration files. "
        "Each of the above configurations has a corresponding "
        "configuration file in shasta/conf." << endl;
}



void shasta::main::listConfiguration(const AssemblerOptions& options)
{
    const string& configName = options.commandLineOnlyOptions.configName;

    if(configName.empty()) {
        throw runtime_error("Specify --config with a valid configuration name.");
    }

    const string* configuration = getConfiguration(configName);
    if(configuration == 0) {
        const string message = configName + " is not a valid configuration name.";
        cout << message << endl;
        listConfigurations();
        throw runtime_error(configName);
    }

    cout << *configuration << flush;
}


