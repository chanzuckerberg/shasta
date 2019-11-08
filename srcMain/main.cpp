// Main program for the Shasta static executable.
// The static executable provides
// basic functionality and reduced performance. 
// For full functionality use the shared library built
// under directory src.

// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "buildId.hpp"
#include "filesystem.hpp"
#include "timestamp.hpp"
namespace shasta {
    namespace main {

        void main(int argumentCount, const char** arguments);

        void assemble(
            Assembler&,
            const AssemblerOptions&,
            vector<string> inputNames);

        void setupRunDirectory(
            const string& memoryMode,
            const string& memoryBacking,
            size_t& pageSize,
            string& dataDirectory
            );

        void setupHugePages();

        // Functions that implement --command keywords
        void assemble(const AssemblerOptions&);
        void saveBinaryData(const AssemblerOptions&);
        void cleanupBinaryData(const AssemblerOptions&);
        void explore(const AssemblerOptions&);
        void createBashCompletionScript(const AssemblerOptions&);

    }
}
using namespace shasta;

// Boost libraries.
#include <boost/program_options.hpp>
#include  <boost/chrono/process_cpu_clocks.hpp>

//  Linux.
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
        cout << timestamp << "Terminated after catching a runtime error exception:" << endl;
        cout << e.what() << endl;
        return 2;
    } catch (const std::bad_alloc& e) {
        cout << timestamp << e.what() << endl;
        cout << "Memory allocation failure." << endl;
        cout << "This assembly requires more memory than available." << endl;
        cout << "Rerun on a larger machine." << endl;
        return 2;
    } catch (const exception& e) {
        cout << timestamp << "Terminated after catching a standard exception:" << endl;
        cout << e.what() << endl;
        return 3;
    } catch (...) {
        cout << timestamp << "Terminated after catching a non-standard exception." << endl;
        return 4;
    }
    return 0;
}



void shasta::main::main(int argumentCount, const char** arguments)
{


    // Parse command line options and the configuration file, if one was specified.
    AssemblerOptions assemblerOptions(argumentCount, arguments);
    cout << buildId() << endl;



    // Execute the requested command.
    if(assemblerOptions.commandLineOnlyOptions.command == "assemble") {
        assemble(assemblerOptions);
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
    }

    // If getting here, the requested command is invalid.
    throw runtime_error("Invalid command " + assemblerOptions.commandLineOnlyOptions.command +
        ". Valid commands are: assemble, saveBinaryData, cleanupBinaryData, createBashCompletionScript.");

}



// Implementation of --command assemble.
void shasta::main::assemble(
    const AssemblerOptions& assemblerOptions)
{
    SHASTA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "assemble");

    const string executableDescription =
        "\nThis is the static executable for the Shasta assembler. "
        "It has no dependencies and requires no installation.\n\n"
        "To run an assembly, use the \"--input\" option to specify the input files. "
        "See below for a description of the other options and parameters.\n\n"
        "Default values of assembly parameters are optimized for an assembly "
        "at coverage 60x. If your data have significantly different coverage, "
        "some changes in assembly parameters may be necessary to get good results.\n\n"
        "For more information about the Shasta assembler, see\n"
        "https://github.com/chanzuckerberg/shasta\n\n"
        "Complete documentation for the latest version of Shasta is available here:\n"
        "https://chanzuckerberg.github.io/shasta\n";



    // Various checks for option validity.

    // Check that a valid assembly strategy was specified.
    switch(assemblerOptions.assemblyOptions.strategy) {
    case 0:
    case 1:
        break;
    default:
        throw runtime_error("Invalid Assembly.strategy " +
            to_string(assemblerOptions.assemblyOptions.strategy));
    }

    // Check that we have at least one input file.
    if(assemblerOptions.commandLineOnlyOptions.inputFileNames.empty()) {
        cout << executableDescription << assemblerOptions.allOptionsDescription << endl;
        throw runtime_error("Specify at least one input file "
            "using command line option \"--input\".");
    }

    // Assembly.useMarginPhase is not supported.
    if(assemblerOptions.assemblyOptions.useMarginPhase) {
        throw runtime_error("Assembly.useMarginPhase is not supported.");
    }

    // If the build does not support GPU acceleration, reject the --gpu option.
#ifndef SHASTA_BUILD_FOR_GPU
    if(assemblerOptions.commandLineOnlyOptions.useGpu) {
        throw runtime_error("This Shasta build does not provide GPU acceleration.");
    }
#endif

    // Check assemblerOptions.minHashOptions.version.
    if( assemblerOptions.minHashOptions.version!=0 and
        assemblerOptions.minHashOptions.version!=1) {
        throw runtime_error("Invalid value " +
            to_string(assemblerOptions.minHashOptions.version) +
            " specified for --MinHash.version. Must be 0 or 1.");
    }

    // If coverage data was requested, memoryMode should be filesystem,
    // otherwise the coverage data cannot be accessed.
    if(assemblerOptions.assemblyOptions.storeCoverageData) {
        if(assemblerOptions.commandLineOnlyOptions.memoryMode != "filesystem") {
            throw runtime_error("To obtain usable coverage data, "
                "you must use --memoryMode filesystem.");
        }
    }

    // MacOS does not support alignment method 1.
#ifndef __linux__
    if( (assemblerOptions.alignOptions.alignMethodForReadGraph == 1) or
        (assemblerOptions.alignOptions.alignMethodForMarkerGraph == 1)) {
        throw runtime_error("Align method 1 is not supported on macOS.");
    }

#endif



    // Write a startup message.
    cout << timestamp <<
        "\nThis is the static executable for the Shasta assembler. "
        "It has no dependencies and requires no installation.\n\n"
        "Default values of assembly parameters are optimized for an assembly "
        "at coverage 60x. If your data have significantly different coverage, "
        "some changes in assembly parameters may be necessary to get good results.\n\n"
        "For more information about the Shasta assembler, see\n"
        "https://github.com/chanzuckerberg/shasta\n\n"
        "Complete documentation for the latest version of Shasta is available here:\n"
        "https://chanzuckerberg.github.io/shasta\n\n";

    // Find absolute paths of the input files.
    // We will use them below after changing directory to the output directory.
    vector<string> inputFileAbsolutePaths;
    for(const string& inputFileName: assemblerOptions.commandLineOnlyOptions.inputFileNames) {
        if(!filesystem::exists(inputFileName)) {
            throw runtime_error("Input file not found: " + inputFileName);
        }
        if(!filesystem::isRegularFile(inputFileName)) {
            throw runtime_error("Input file is not a regular file: " + inputFileName);
        }
        inputFileAbsolutePaths.push_back(filesystem::getAbsolutePath(inputFileName));
    }



    // Create the run the output directory. If it exists, stop.
    if(filesystem::exists(assemblerOptions.commandLineOnlyOptions.assemblyDirectory)) {
        throw runtime_error(
            "Assembly directory " +
            assemblerOptions.commandLineOnlyOptions.assemblyDirectory +
            " already exists.\n"
            "Remove it or use --assemblyDirectory to specify a different assembly directory.");
    }
    filesystem::createDirectory(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);

    // Make the output directory current.
    filesystem::changeDirectory(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);



    // Set up the run directory as required by the memoryMode and memoryBacking options.
    size_t pageSize = 0;
    string dataDirectory;
    setupRunDirectory(
        assemblerOptions.commandLineOnlyOptions.memoryMode,
        assemblerOptions.commandLineOnlyOptions.memoryBacking,
        pageSize,
        dataDirectory);



    // Write out the option values we are using.
    // This code should should be moved to AssemblerOptions::write.
    cout << "Options in use:" << endl;
    cout << "Input files: ";
    copy(
        assemblerOptions.commandLineOnlyOptions.inputFileNames.begin(),
        assemblerOptions.commandLineOnlyOptions.inputFileNames.end(),
        ostream_iterator<string>(cout, " "));
    cout << endl;
    cout << "assemblyDirectory = " <<
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory << endl;
#ifdef __linux__
    cout << "memoryMode = " << assemblerOptions.commandLineOnlyOptions.memoryMode << endl;
    cout << "memoryBacking = " << assemblerOptions.commandLineOnlyOptions.memoryBacking << endl;
    cout << "threadCount = " << assemblerOptions.commandLineOnlyOptions.threadCount << "\n" << endl;
    cout << "useGpu = " <<
        AssemblerOptions::convertBoolToPythonString(assemblerOptions.commandLineOnlyOptions.useGpu) <<
        "\n" << endl;
#endif
    assemblerOptions.write(cout);
    {
        ofstream configurationFile("shasta.conf");
        assemblerOptions.write(configurationFile);
    }



    // Initial disclaimer message.
#ifdef __linux
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
#else
    cout << "The macOS version of the Shasta assembler runs at degraded performance.\n";
    cout << "Use Linux for full performance.\n";
    cout << "Therefore the results of this run should not be used\n"
        "for benchmarking purposes." << endl;
#endif

    // Create the Assembler.
    Assembler assembler(dataDirectory, true, pageSize);

    // Run the assembly.
    assemble(assembler, assemblerOptions, inputFileAbsolutePaths);

    // Final disclaimer message.
#ifdef __linux
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
#else
    cout << "The macOS version of the Shasta assembler runs at degraded performance.\n";
    cout << "Use Linux for full performance.\n";
    cout << "Therefore the results of this run should not be used\n"
        "for benchmarking purposes." << endl;
#endif

    // Write out the build id again.
    cout << buildId() << endl;

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
            filesystem::createDirectory("Data");
            dataDirectory = "Data/";
            pageSize = 4096;

        } else if(memoryBacking == "4K") {

            // Binary files on the tmpfs filesystem
            // (filesystem in memory backed by 4K pages).
            // This requires root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            filesystem::createDirectory("Data");
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
            filesystem::createDirectory("Data");
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
    cout << timestamp << "Begin loading reads from " << inputFileNames.size() << " files." << endl;
    const auto t0 = steady_clock::now();
    for(const string& inputFileName: inputFileNames) {
        assembler.addReads(
            inputFileName,
            assemblerOptions.readsOptions.minReadLength,
            threadCount);
    }
    if(assembler.readCount() == 0) {
        throw runtime_error("There are no input reads.");
    }
    const auto t1 = steady_clock::now();
    cout << timestamp << "Done loading reads from " << inputFileNames.size() << " files." << endl;
    cout << "Read loading took " << seconds(t1-t0) << "s." << endl;



    // Initialize read flags.
    assembler.initializeReadFlags();

    // Create a histogram of read lengths.
    assembler.histogramReadLength("ReadLengthHistogram.csv");

    // Select the k-mers that will be used as markers.
    if(assemblerOptions.kmersOptions.suppressHighFrequencyMarkers) {
        assembler.selectKmersBasedOnFrequency(
            assemblerOptions.kmersOptions.k,
            assemblerOptions.kmersOptions.probability, 231,
            assemblerOptions.kmersOptions.enrichmentThreshold, threadCount);
    } else {
        assembler.randomlySelectKmers(
            assemblerOptions.kmersOptions.k,
            assemblerOptions.kmersOptions.probability, 231);
    }

    // Find the markers in the reads.
    assembler.findMarkers(0);

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



    // Find alignment candidates.
    if(assemblerOptions.minHashOptions.allPairs) {
        assembler.markAlignmentCandidatesAllPairs();
    } else if(assemblerOptions.minHashOptions.version == 0) {
        assembler.findAlignmentCandidatesLowHash0(
            assemblerOptions.minHashOptions.m,
            assemblerOptions.minHashOptions.hashFraction,
            assemblerOptions.minHashOptions.minHashIterationCount,
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



    // Compute alignments.
    if(assemblerOptions.commandLineOnlyOptions.useGpu) {
#ifdef SHASTA_BUILD_FOR_GPU
        cout << "Using GPU acceleration for alignment computation." << endl;
        cout << "This is under development and is not ready to be used." << endl;
        assembler.computeAlignmentsGpu(
            assemblerOptions.alignOptions.maxMarkerFrequency,
            assemblerOptions.alignOptions.maxSkip,
            assemblerOptions.alignOptions.maxDrift,
            assemblerOptions.alignOptions.minAlignedMarkerCount,
            assemblerOptions.alignOptions.maxTrim,
            threadCount);
#else
        throw runtime_error("This Shasta build does not provide GPU acceleration.");
#endif
    } else {
        assembler.computeAlignments(
            assemblerOptions.alignOptions.alignMethodForReadGraph,
            assemblerOptions.alignOptions.maxMarkerFrequency,
            assemblerOptions.alignOptions.maxSkip,
            assemblerOptions.alignOptions.maxDrift,
            assemblerOptions.alignOptions.minAlignedMarkerCount,
            assemblerOptions.alignOptions.maxTrim,
            assemblerOptions.alignOptions.matchScore,
            assemblerOptions.alignOptions.mismatchScore,
            assemblerOptions.alignOptions.gapScore,
            threadCount);
    }



    // Create the read graph.
    assembler.createReadGraph(
        assemblerOptions.readGraphOptions.maxAlignmentCount,
        assemblerOptions.alignOptions.maxTrim);

    // Flag read graph edges that cross strands.
    assembler.flagCrossStrandReadGraphEdges(
        assemblerOptions.readGraphOptions.crossStrandMaxDistance,
        threadCount);

    // Flag chimeric reads.
    assembler.flagChimericReads(assemblerOptions.readGraphOptions.maxChimericReadDistance, threadCount);
    assembler.computeReadGraphConnectedComponents(assemblerOptions.readGraphOptions.minComponentSize);



    // Create marker graph vertices.
    // This uses a disjoint sets data structure to merge markers
    // that are aligned based on an alignment presentin the read graph.
    if(assemblerOptions.commandLineOnlyOptions.useGpu) {

#ifdef SHASTA_BUILD_FOR_GPU

        // Create marker graph vertices: do it on the GPU.
        cout << "Using GPU acceleration for creating marker graph vertices.." << endl;
        cout << "This is under development and is not ready to be used." << endl;
        assembler.createMarkerGraphVerticesGpu(
            assemblerOptions.alignOptions.maxMarkerFrequency,
            assemblerOptions.alignOptions.maxSkip,
            assemblerOptions.alignOptions.maxDrift,
            assemblerOptions.markerGraphOptions.minCoverage,
            assemblerOptions.markerGraphOptions.maxCoverage,
            threadCount);
#else

        // The build does not have GPU support.
        throw runtime_error("This Shasta build does not provide GPU acceleration.");

#endif
    } else {

        // Create marker graph vertices: mainstream code.
        assembler.createMarkerGraphVertices(
            assemblerOptions.alignOptions.alignMethodForMarkerGraph,
            assemblerOptions.alignOptions.maxMarkerFrequency,
            assemblerOptions.alignOptions.maxSkip,
            assemblerOptions.alignOptions.maxDrift,
            assemblerOptions.alignOptions.matchScore,
            assemblerOptions.alignOptions.mismatchScore,
            assemblerOptions.alignOptions.gapScore,
            assemblerOptions.markerGraphOptions.minCoverage,
            assemblerOptions.markerGraphOptions.maxCoverage,
            threadCount);
    }



    // Find the reverse complement of each marker graph vertex.
    assembler.findMarkerGraphReverseComplementVertices(threadCount);

    // Create edges of the marker graph.
    assembler.createMarkerGraphEdges(threadCount);
    assembler.findMarkerGraphReverseComplementEdges(threadCount);

    // Approximate transitive reduction.
    assembler.transitiveReduction(
        assemblerOptions.markerGraphOptions.lowCoverageThreshold,
        assemblerOptions.markerGraphOptions.highCoverageThreshold,
        assemblerOptions.markerGraphOptions.maxDistance,
        assemblerOptions.markerGraphOptions.edgeMarkerSkipThreshold);



    // Marker graph processing that varies according to
    // --Assembly.strategy.
    if(assemblerOptions.assemblyOptions.strategy == 1) {

        // Do a reverse transitive reduction to remove
        // the short reverse bubbles, but don't do
        // bubble/superbubble removal!
        assembler.reverseTransitiveReduction(
            assemblerOptions.markerGraphOptions.lowCoverageThreshold,
            assemblerOptions.markerGraphOptions.highCoverageThreshold,
            assemblerOptions.markerGraphOptions.maxDistance);

        // Prune the marker graph.
        assembler.pruneMarkerGraphStrongSubgraph(
            assemblerOptions.markerGraphOptions.pruneIterationCount);

        // Create the assembly graph.
        assembler.createAssemblyGraphEdges();
        assembler.createAssemblyGraphVertices();

        // Remove low coverage cross edges of the assembly graph
        // and the corresponding marker graph edges.
        assembler.removeLowCoverageCrossEdges(
            assemblerOptions.assemblyOptions.crossEdgeCoverageThreshold);

        // Compute marker graph coverage histogram.
        assembler.computeMarkerGraphCoverageHistogram();

        // Recreate the assembly graph, to
        // allow edges to be combined after cross edge removal.
        assembler.assemblyGraph.remove();
        assembler.createAssemblyGraphEdges();
        assembler.createAssemblyGraphVertices();

    } else {
        SHASTA_ASSERT(assemblerOptions.assemblyOptions.strategy == 0);

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
    }



    assembler.writeAssemblyGraph("AssemblyGraph-Final.dot");

    // Compute optimal repeat counts for each vertex of the marker graph.
    assembler.assembleMarkerGraphVertices(threadCount);

    // If coverage data was requested, compute and store coverage data for the vertices.
    if(assemblerOptions.assemblyOptions.storeCoverageData or
        assemblerOptions.assemblyOptions.storeCoverageDataCsvLengthThreshold>0) {
        assembler.computeMarkerGraphVerticesCoverageData(threadCount);
    }

    // Compute consensus sequence for marker graph edges to be used for assembly.
    assembler.assembleMarkerGraphEdges(
        threadCount,
        assemblerOptions.assemblyOptions.markerGraphEdgeLengthThresholdForConsensus,
        false,
        assemblerOptions.assemblyOptions.storeCoverageData or
        assemblerOptions.assemblyOptions.storeCoverageDataCsvLengthThreshold>0
        );

    // Use the assembly graph for global assembly.
    assembler.assemble(
        threadCount,
        assemblerOptions.assemblyOptions.storeCoverageDataCsvLengthThreshold);
    assembler.computeAssemblyStatistics();
    assembler.writeGfa1("Assembly.gfa");
    assembler.writeGfa1BothStrands("Assembly-BothStrands.gfa");
    assembler.writeFasta("Assembly.fasta");

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

    cout << "Assembly time statistics:\n"
        "    Elapsed seconds: " << elapsedTime << "\n"
        "    Elapsed minutes: " << elapsedTime/60. << "\n"
        "    Elapsed hours:   " << elapsedTime/3600. << "\n";
    cout << "Average CPU utilization: " << averageCpuUtilization << endl;

    // Write the assembly summary.
    ofstream html("AssemblySummary.html");
    assembler.writeAssemblySummary(html);
    ofstream json("AssemblySummary.json");
    assembler.writeAssemblySummaryJson(json);

    // Also write a summary of read information.
    assembler.writeReadsSummary();

    // For Assembly strategy 1,write a warning message
    if(assemblerOptions.assemblyOptions.strategy == 1) {
        cout << "This assembly was created using --Assembly.strategy 1. "
            "This option is under development and "
            " does not produce usable assembly results." << endl;
        return;
    }
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
    if(!filesystem::exists(dataDirectory)) {
        cout << dataDirectory << " does not exist, nothing done." << endl;
        return;
    }

    // Check that the DataOnDisk directory does not exist.
    const string dataOnDiskDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/DataOnDisk";
    if(filesystem::exists(dataOnDiskDirectory)) {
        cout << dataOnDiskDirectory << " already exists, nothing done." << endl;
        return;
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
    if(!filesystem::exists(dataDirectory)) {
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
    if(filesystem::exists(dataOnDiskDirectory)) {
        filesystem::changeDirectory(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);
        const string command = "ln -s DataOnDisk Data";
        ::system(command.c_str());
    }

}


// Implementation of --command explore.
void shasta::main:: explore(
    const AssemblerOptions& assemblerOptions)
{
    // Go to the assembly directory.
    filesystem::changeDirectory(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);
    
    // Check that we have the binary data. 
    if(!filesystem::exists("Data")) {
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
    Assembler assembler("Data/", false, 0);
    
    // Access all available binary data.
    assembler.accessAllSoft();
    
    // Set up the consensus caller.
    cout << "Setting up consensus caller " <<
        assemblerOptions.assemblyOptions.consensusCaller << endl;
    assembler.setupConsensusCaller(assemblerOptions.assemblyOptions.consensusCaller);

 
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
// which makes it esaier to type long option names.
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
// PLESE MAKE IT BETTER AND SUBMIT A PULL REQUEST!
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

    // Other keywords. This should be modified to only accept them after the appropriate option.
    file << "assemble saveBinaryData cleanupBinaryData explore createBashCompletionScript \\\n";
    file << "filesystem anonymous \\\n";
    file << "disk 4K 2M \\\n";
    file << "user local unrestricted \\\n";

    // Finish the "complete" command.
    file << "\" shasta\n";
}
