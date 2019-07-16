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
        void runAssembly(
            Assembler&,
            const AssemblerOptions&,
            vector<string> inputFastaFileNames,
            const string& consensusCallerType);
        void setupHugePages();
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
#include "fstream.hpp"
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
    cout << buildId() << endl;

    // Some names in the boost program_options library.
    using boost::program_options::command_line_parser;
    using boost::program_options::options_description;
    using boost::program_options::value;
    using boost::program_options::variables_map;

    const string executableDescription =
        "\nThis is the static executable for the Shasta assembler. "
        "It provides limited Shasta functionality, "
        "at reduced performance when using the default options, "
        "but has no dependencies and requires no installation.\n\n"
        "To run an assembly, use the \"--input\" option to specify the input Fasta files. "
        "See below for a description of the other options and parameters.\n\n"
        "Default values of assembly parameters are optimized for an assembly "
        "at coverage 60x. If your data have significantly different coverage, "
        "some changes in assembly parameters may be necessary to get good results.\n\n"
        "For more information about the Shasta assembler, see\n"
        "https://github.com/chanzuckerberg/shasta\n\n"
        "Complete documentation for the latest version of Shasta is available here:\n"
        "https://chanzuckerberg.github.io/shasta\n";



    // Options that are only allowed on the command line.
    options_description commandLineOnlyOptions(
        "Options allowed only on the command line");
    string configFileName;
    vector < string > inputFastaFileNames;
    string outputDirectory;
    string command;
    string memoryMode;
    string memoryBacking;
    commandLineOnlyOptions.add_options()

        ("help", 
        "Write a help message.")

        ("config", 
        boost::program_options::value<string>(&configFileName),
        "Configuration file name.")

        ("input",
        value< vector<string> >(&inputFastaFileNames)->multitoken(),
        "Names of input FASTA files. Specify at least one.")

        ("output",
        value<string>(&outputDirectory)->
        default_value("ShastaRun"),
        "Name of the output directory. Must not exist.")

        ("command",
        value<string>(&command)->
        default_value("assemble"),
        "Command to run. Must be one of:\n"
        "assemble (default): run an assembly\n"
        "cleanup: cleanup the Data directory that was created during assembly\n"
        "    if --memoryMode filesystem.\n")

#ifdef __linux__
        ("memoryMode",
        value<string>(&memoryMode)->
        default_value("anonymous"),
        "Specify whether allocated memory is anonymous or backed by a filesystem. "
        "Allowed values: anonymous (default), filesystem.")

        ("memoryBacking",
        value<string>(&memoryBacking)->
        default_value("4K"),
        "Specify the type of pages used to back memory.\n"
        "Allowed values: disk, 4K (default), 2M (for best performance, Linux only). "
        "All combinations (memoryMode, memoryBacking) are allowed "
        "except for (anonymous, disk).\n"
        "Some combinations require root privilege, which is obtained using sudo "
        "and may result in a password prompting depending on your sudo set up.")
#endif
        ;


    // For reasons not completely understood and that there was no time to investigate,
    // the only combination that works on MacOS is
    // "--memoryMode filesystem --memoryBacking disk".
    // This incurs a performance price but this is not too much of a big deal
    // as macOS  is only to be used for small test runs.
#ifndef __linux__
    memoryMode = "filesystem";
    memoryBacking = "disk";
#endif



    // Options allowed on the command line and in the config file.
    // Values specified in the command line take precedence.
    options_description options(
        "Options allowed on the command line and in the config file");
    AssemblerOptions assemblerOptions;
    assemblerOptions.add(options);
    

        
    // Get options from the command line.
    // These take precedence over values entered in the config file.
    options_description commandLineOptions;
    commandLineOptions.add(commandLineOnlyOptions);
    commandLineOptions.add(options);
    variables_map variablesMap;
    store(command_line_parser(argumentCount, arguments).
          options(commandLineOptions).run(), variablesMap);
    notify(variablesMap);

    if (variablesMap.count("help")) {
        cout << executableDescription << commandLineOptions << endl;
        return;
    }

    // Get options from the config file, if one was specified.
    if(!configFileName.empty()) {
        ifstream configFile(configFileName);
        if (!configFile) {
            throw runtime_error("Unable to open open config file " + configFileName);
        }
        store(parse_config_file(configFile, options), variablesMap);
        notify(variablesMap);
    }



    // If command is "cleanup", just do it and exit.
    if(command == "cleanup") {
        const string dataDirectory = outputDirectory + "/Data";
        if(!filesystem::exists(dataDirectory)) {
            cout << dataDirectory << " does not exist, nothing done." << endl;
        }
        ::system(("sudo umount " + dataDirectory).c_str());
        const int errorCode = ::system(string("rm -rf " + dataDirectory).c_str());
        if(errorCode != 0) {
            throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
                " removing " + dataDirectory);
        }
        cout << "Cleanup of " << dataDirectory << " successful." << endl;
        return;
    }



    // Check that we have at least one input FASTA file.     
    if (inputFastaFileNames.empty()) {
        cout << executableDescription << commandLineOptions << endl;
        throw runtime_error("Specify at least one input FASTA file "
            "using command line option \"--input\".");
    }

    // Parse MarkerGraph.simplifyMaxLength.
    assemblerOptions.markerGraphOptions.parseSimplifyMaxLength();

    // Check for options unsupported by the static executable.
    if(assemblerOptions.Assembly.useMarginPhase != "False") {
        throw runtime_error("Assembly.useMarginPhase is not supported by the Shasta static executable.");
    }

    // Write a startup message.
    cout << timestamp <<
        "\nThis is the static executable for the Shasta assembler. "
        "It provides limited Shasta functionality, "
        "at reduced performance when using the default options, "
        "but has no dependencies and requires no installation.\n\n"
        "Default values of assembly parameters are optimized for an assembly "
        "at coverage 60x. If your data have significantly different coverage, "
        "some changes in assembly parameters may be necessary to get good results.\n\n"
        "For more information about the Shasta assembler, see\n"
        "https://github.com/chanzuckerberg/shasta\n\n"
        "Complete documentation for the latest version of Shasta is available here:\n"
        "https://chanzuckerberg.github.io/shasta\n\n";

    // Find absolute paths of the input fasta files.
    // We will use them below after changing directory to the output directory.
    vector<string> inputFastaFileAbsolutePaths;
    for(const string& inputFastaFileName: inputFastaFileNames) {
        inputFastaFileAbsolutePaths.push_back(filesystem::getAbsolutePath(inputFastaFileName));
    }



    // If necessary, find the absolute path of the configuration file
    // for the consensus caller.
    string consensusCallerType;
    string consensusCallerConfigurationFileName;
    if(assemblerOptions.Assembly.consensusCaller == "SimpleConsensusCaller") {
        consensusCallerType = "SimpleConsensusCaller";
    } else {

        // In this case, the option specifies the name
        // of the configuration file to use for the consensus caller.
        // The type is deduced from the file name.
        consensusCallerConfigurationFileName = assemblerOptions.Assembly.consensusCaller;

        // Deduce the type from the file name.
        const size_t lastSlashPosition = consensusCallerConfigurationFileName.find_last_of('/');
        const size_t fileNameBegin = (lastSlashPosition == string::npos ? 0 : lastSlashPosition+1);
        const string fileName = consensusCallerConfigurationFileName.substr(fileNameBegin);
        const size_t dashPosition = fileName.find_first_of('-');
        consensusCallerType = fileName.substr(0, dashPosition);

        // Check that it is one of the supported types.
        if(consensusCallerType != "SimpleBayesianConsensusCaller") {
            throw runtime_error("Invalid consensus caller " +
                assemblerOptions.Assembly.consensusCaller);
        }

    }



    // Create the run the output directory. If it exists, stop.
    if(filesystem::exists(outputDirectory)) {
        throw runtime_error("Output directory " + outputDirectory + " already exists.\n"
            "Remove it or use --output to specify a different output directory.");
    }
    filesystem::createDirectory(outputDirectory);

    // If necessary, copy the configuration file for the consensus caller
    // to the output directory.
    if(consensusCallerConfigurationFileName.size() > 0) {
        filesystem::copy(
            consensusCallerConfigurationFileName,
            outputDirectory + "/" + consensusCallerType + ".csv");
    }

    // Make the output directory current.
    filesystem::changeDirectory(outputDirectory);



    // Set up the run directory as required by the memoryMode and memoryBacking option.
    size_t pageSize = 0;
    string dataDirectory;
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



    // Write out the option values we are using.
    cout << "Options in use:" << endl;
    cout << "Input FASTA files: ";
    copy(inputFastaFileNames.begin(), inputFastaFileNames.end(), ostream_iterator<string>(cout, " "));
    cout << endl;
    cout << "outputDirectory = " << outputDirectory << endl;
#ifdef __linux__
    cout << "memoryMode = " << memoryMode << endl;
    cout << "memoryBacking = " << memoryBacking << "\n" << endl;
#endif
    assemblerOptions.write(cout);
    {
        ofstream configurationFile("shasta.conf");
        assemblerOptions.write(configurationFile);
    }

    // Initial disclaimer message.
#ifdef __linux
    if(memoryBacking != "2M" && memoryMode != "filesystem") {
        cout << "This run uses options \"--memoryBacking " << memoryBacking <<
            " --memoryMode " << memoryMode << "\".\n"
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

    // Write the number of threads the run will use.
    cout << "This assembly will use " << std::thread::hardware_concurrency() << " threads." << endl;

    // Create the Assembler.
    Assembler assembler(dataDirectory, true, pageSize);

    // Run the assembly.
    runAssembly(assembler, assemblerOptions, inputFastaFileAbsolutePaths, consensusCallerType);

    // Final disclaimer message.
#ifdef __linux
    if(memoryBacking != "2M" && memoryMode != "filesystem") {
        cout << "This run used options \"--memoryBacking " << memoryBacking <<
            " --memoryMode " << memoryMode << "\".\n"
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



// This runs the entire assembly, under the following assumptions:
// - The current directory is the run directory.
// - The Data directory has already been created and set up, if necessary.
// - The input Fasta file names are either absolute,
//   or relative to the run directory, which is the current directory.
void shasta::main::runAssembly(
    Assembler& assembler,
    const AssemblerOptions& assemblerOptions,
    vector<string> inputFastaFileNames,
    const string& consensusCallerType)
{
    const auto steadyClock0 = std::chrono::steady_clock::now();
    const auto userClock0 = boost::chrono::process_user_cpu_clock::now();
    const auto systemClock0 = boost::chrono::process_system_cpu_clock::now();

    // The executable only supports SimpleConsensusCaller,
    // at least for now.
    cout << "Setting up " << consensusCallerType << endl;
    assembler.setupConsensusCaller(consensusCallerType);

    // Add reads from the specified FASTA files.
    for(const string& inputFastaFileName: inputFastaFileNames) {
        assembler.addReadsFromFasta(
            inputFastaFileName,
            assemblerOptions.readsOptions.minReadLength,
            2ULL * 1024ULL * 1024ULL * 1024ULL,
            1,
            0);
    }
    if(assembler.readCount() == 0) {
        throw runtime_error("There are no input reads.");
    }


    // Initialize read flags.
    assembler.initializeReadFlags();

    // Create a histogram of read lengths.
    assembler.histogramReadLength("ReadLengthHistogram.csv");

    // Randomly select the k-mers that will be used as markers.
    assembler.randomlySelectKmers(
        assemblerOptions.kmersOptions.k,
        assemblerOptions.kmersOptions.probability, 231);

    // Find the markers in the reads.
    assembler.findMarkers(0);

    // Flag palindromic reads.
    // These wil be excluded from further processing.
    assembler.flagPalindromicReads(
        assemblerOptions.readsOptions.palindromicReads.maxSkip,
        assemblerOptions.readsOptions.palindromicReads.maxMarkerFrequency,
        assemblerOptions.readsOptions.palindromicReads.alignedFractionThreshold,
        assemblerOptions.readsOptions.palindromicReads.nearDiagonalFractionThreshold,
        assemblerOptions.readsOptions.palindromicReads.deltaThreshold,
        0);

    // Find alignment candidates.
    assembler.findAlignmentCandidatesLowHash(
        assemblerOptions.minHashOptions.m,
        assemblerOptions.minHashOptions.hashFraction,
        assemblerOptions.minHashOptions.minHashIterationCount,
        0,
        assemblerOptions.minHashOptions.maxBucketSize,
        assemblerOptions.minHashOptions.minFrequency,
        0);


    // Compute alignments.
    assembler.computeAlignments(
        assemblerOptions.alignOptions.maxMarkerFrequency,
        assemblerOptions.alignOptions.maxSkip,
        assemblerOptions.alignOptions.minAlignedMarkerCount,
        assemblerOptions.alignOptions.maxTrim,
        0);

    // Create the read graph.
    assembler.createReadGraph(
        assemblerOptions.readGraphOptions.maxAlignmentCount,
        assemblerOptions.alignOptions.maxTrim);

    // Flag read graph edges that cross strands.
    assembler.flagCrossStrandReadGraphEdges();

    // Flag chimeric reads.
    assembler.flagChimericReads(assemblerOptions.readGraphOptions.maxChimericReadDistance, 0);
    assembler.computeReadGraphConnectedComponents(assemblerOptions.readGraphOptions.minComponentSize);

    // Create vertices of the marker graph.
    assembler.createMarkerGraphVertices(
        assemblerOptions.alignOptions.maxMarkerFrequency,
        assemblerOptions.alignOptions.maxSkip,
        assemblerOptions.markerGraphOptions.minCoverage,
        assemblerOptions.markerGraphOptions.maxCoverage,
        0);
    assembler.findMarkerGraphReverseComplementVertices(0);

    // Create edges of the marker graph.
    assembler.createMarkerGraphEdges(0);
    assembler.findMarkerGraphReverseComplementEdges(0);

    // Approximate transitive reduction.
    assembler.flagMarkerGraphWeakEdges(
        assemblerOptions.markerGraphOptions.lowCoverageThreshold,
        assemblerOptions.markerGraphOptions.highCoverageThreshold,
        assemblerOptions.markerGraphOptions.maxDistance,
        assemblerOptions.markerGraphOptions.edgeMarkerSkipThreshold);

    // Prune the strong subgraph of the marker graph.
    assembler.pruneMarkerGraphStrongSubgraph(
        assemblerOptions.markerGraphOptions.pruneIterationCount);

    // Simplify the marker graph to remove bubbles and superbubbles.
    // The maxLength parameter controls the maximum number of markers
    // for a branch to be collapsed during each iteration.
    assembler.simplifyMarkerGraph(assemblerOptions.markerGraphOptions.simplifyMaxLengthVector, false);

    // Create the assembly graph.
    assembler.createAssemblyGraphEdges();
    assembler.createAssemblyGraphVertices();
    assembler.writeAssemblyGraph("AssemblyGraph-Final.dot");

    // Compute optimal repeat counts for each vertex of the marker graph.
    assembler.assembleMarkerGraphVertices(0);

    // If coverage data was requested, compute and store coverage data for the vertices.
    if(assemblerOptions.Assembly.storeCoverageData != "False") {
        assembler.computeMarkerGraphVerticesCoverageData(0);
    }

    // Compute consensus sequence for marker graph edges to be used for assembly.
    assembler.assembleMarkerGraphEdges(
        0,
        assemblerOptions.Assembly.markerGraphEdgeLengthThresholdForConsensus,
        false,
        assemblerOptions.Assembly.storeCoverageData != "False");

    // Use the assembly graph for global assembly.
    assembler.assemble(0);
    assembler.computeAssemblyStatistics();
    assembler.writeGfa1("Assembly.gfa");
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
