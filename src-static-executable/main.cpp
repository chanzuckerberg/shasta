// Main program for the Shasta static executable.
// The static executable provides
// basic functionality and reduced performance. 
// For full functionality use the shared library built
// under directory src.

// Shasta.
#include "Assembler.hpp"
#include "AssemblyOptions.hpp"
#include "buildId.hpp"
#include "filesystem.hpp"
#include "timestamp.hpp"
namespace ChanZuckerberg {
    namespace shasta {
        namespace main {
            void main(int argumentCount, const char** arguments);
            void runAssembly(
                Assembler&,
                const AssemblyOptions&,
                vector<string> inputFastaFileNames);
            void setupHugePages();
        }
        class AssemblyOptions;
    }
}
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/program_options.hpp>

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

    } catch(boost::program_options::error_with_option_name e) {
        cout << "Invalid option: " << e.what() << endl;
        return 1;
    } catch (runtime_error e) {
        cout << timestamp << "Terminated after catching a runtime error exception:" << endl;
        cout << e.what() << endl;
        return 2;
    } catch (exception e) {
        cout << timestamp << "Terminated after catching a standard exception:" << endl;
        cout << e.what() << endl;
        return 3;
    } catch (...) {
        cout << timestamp << "Terminated after catching a non-standard exception." << endl;
        return 4;
    }
    return 0;
}



void ChanZuckerberg::shasta::main::main(int argumentCount, const char** arguments)
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
        "at reduced performance when using the default options,"
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
        ;



    // Options allowed on the command line and in the config file.
    // Values specified in the command line take precedence.
    options_description options(
        "Options allowed on the command line and in the config file");
    AssemblyOptions assemblyOptions;
    assemblyOptions.add(options);
    

        
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
        throw runtime_error("Specify at least one input FASTA file.");
    }

    // Parse MarkerGraph.simplifyMaxLength.
    assemblyOptions.MarkerGraph.parseSimplifyMaxLength();

    // Check for options unsupported by the static executable.
    if(assemblyOptions.Assembly.consensusCaller != "SimpleConsensusCaller") {
        throw runtime_error("Assembly.consensusCaller value " + assemblyOptions.Assembly.consensusCaller +
            " is not supported.\n"
            "Only value supported by the Shasta static executable is SimpleConsensusCaller.");
    }
    if(assemblyOptions.Assembly.useMarginPhase != "False") {
        throw runtime_error("Assembly.useMarginPhase is not supported by the Shasta static executable.");
    }
    if(assemblyOptions.Assembly.storeCoverageData != "False") {
        throw runtime_error("Assembly.storeCoverageData is not supported by the Shasta static executable.");
    }

    // Write a startup message.
    cout << timestamp <<
        "\nThis is the static executable for the Shasta assembler. "
        "It provides limited Shasta functionality "
        "at but has no dependencies and requires no installation.\n\n"
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

    // If the output directory exists, stop.
    // Otherwise, create it and make it current.
    if(filesystem::exists(outputDirectory)) {
        throw runtime_error("Output directory " + outputDirectory + " already exists.\n"
            "Remove it or use --output to specify a different output directory.");
    }
    filesystem::createDirectory(outputDirectory);
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
    cout << "memoryMode = " << memoryMode << endl;
    cout << "memoryBacking = " << memoryBacking << "\n" << endl;
    assemblyOptions.write(cout);
    {
        ofstream configurationFile("shasta.conf");
        assemblyOptions.write(configurationFile);
    }

    // Create the Assembler.
    Assembler assembler(dataDirectory, true, pageSize);

    // Run the assembly.
    runAssembly(assembler, assemblyOptions, inputFastaFileAbsolutePaths);

    // Final disclaimer message.
#ifdef __linux
    if(memoryBacking != "2M" && memoryMode != "filesystem") {
        cout << "This run was done with \"--memoryBacking " << memoryBacking <<
            " --memoryMode " << memoryMode << "\".\n"
            "This could have resulted in performance degradation.\n"
            "For full performance, use \"--memoryBacking 2M --memoryMode filesystem\"\n"
            "(root privilege via sudo required).\n"
            "Therefore the results of this run should not be used\n"
            "for benchmarking purposes." << endl;
    }
#endif

    // Write out the build id again.
    cout << buildId() << endl;

}



// This runs the entire assembly, under the following assumptions:
// - The current directory is the run directory.
// - The Data directory has already been created and set up, if necessary.
// - The input Fasta file names are either absolute,
//   or relative to the run directory, which is the current directory.
void ChanZuckerberg::shasta::main::runAssembly(
    Assembler& assembler,
    const AssemblyOptions& assemblyOptions,
    vector<string> inputFastaFileNames)
{
    // The executable only supports SimpleConsensusCaller,
    // at least for now.
    assembler.setupConsensusCaller("SimpleConsensusCaller");

    // Add reads from the specified FASTA files.
    for(const string& inputFastaFileName: inputFastaFileNames) {
        assembler.addReadsFromFasta(
            inputFastaFileName,
            assemblyOptions.Reads.minReadLength,
            2ULL * 1024ULL * 1024ULL * 1024ULL,
            1,
            0);
    }


    // Initialize read flags.
    assembler.initializeReadFlags();

    // Create a histogram of read lengths.
    assembler.histogramReadLength("ReadLengthHistogram.csv");

    // Randomly select the k-mers that will be used as markers.
    assembler.randomlySelectKmers(
        assemblyOptions.Kmers.k,
        assemblyOptions.Kmers.probability, 231);

    // Find the markers in the reads.
    assembler.findMarkers(0);

    // Flag palindromic reads.
    // These wil be excluded from further processing.
    assembler.flagPalindromicReads(
        assemblyOptions.Reads.palindromicReads.maxSkip,
        assemblyOptions.Reads.palindromicReads.maxMarkerFrequency,
        assemblyOptions.Reads.palindromicReads.alignedFractionThreshold,
        assemblyOptions.Reads.palindromicReads.nearDiagonalFractionThreshold,
        assemblyOptions.Reads.palindromicReads.deltaThreshold,
        0);

    // Find alignment candidates.
    assembler.findAlignmentCandidatesLowHash(
        assemblyOptions.MinHash.m,
        assemblyOptions.MinHash.hashFraction,
        assemblyOptions.MinHash.minHashIterationCount,
        0,
        assemblyOptions.MinHash.maxBucketSize,
        assemblyOptions.MinHash.minFrequency,
        0);


    // Compute alignments.
    assembler.computeAlignments(
        assemblyOptions.Align.maxMarkerFrequency,
        assemblyOptions.Align.maxSkip,
        assemblyOptions.Align.minAlignedMarkerCount,
        assemblyOptions.Align.maxTrim,
        0);

    // Create the read graph.
    assembler.createReadGraph(
        assemblyOptions.ReadGraph.maxAlignmentCount,
        assemblyOptions.Align.maxTrim);

    // Flag read graph edges that cross strands.
    assembler.flagCrossStrandReadGraphEdges();

    // Flag chimeric reads.
    assembler.flagChimericReads(assemblyOptions.ReadGraph.maxChimericReadDistance, 0);
    assembler.computeReadGraphConnectedComponents(assemblyOptions.ReadGraph.minComponentSize);

    // Create vertices of the marker graph.
    assembler.createMarkerGraphVertices(
        assemblyOptions.Align.maxMarkerFrequency,
        assemblyOptions.Align.maxSkip,
        assemblyOptions.MarkerGraph.minCoverage,
        assemblyOptions.MarkerGraph.maxCoverage,
        0);
    assembler.findMarkerGraphReverseComplementVertices(0);

    // Create edges of the marker graph.
    assembler.createMarkerGraphEdges(0);
    assembler.findMarkerGraphReverseComplementEdges(0);

    // Approximate transitive reduction.
    assembler.flagMarkerGraphWeakEdges(
        assemblyOptions.MarkerGraph.lowCoverageThreshold,
        assemblyOptions.MarkerGraph.highCoverageThreshold,
        assemblyOptions.MarkerGraph.maxDistance,
        assemblyOptions.MarkerGraph.edgeMarkerSkipThreshold);

    // Prune the strong subgraph of the marker graph.
    assembler.pruneMarkerGraphStrongSubgraph(
        assemblyOptions.MarkerGraph.pruneIterationCount);

    // Simplify the marker graph to remove bubbles and superbubbles.
    // The maxLength parameter controls the maximum number of markers
    // for a branch to be collapsed during each iteration.
    assembler.simplifyMarkerGraph(assemblyOptions.MarkerGraph.simplifyMaxLengthVector, false);

    // Create the assembly graph.
    assembler.createAssemblyGraphEdges();
    assembler.createAssemblyGraphVertices();
    assembler.writeAssemblyGraph("AssemblyGraph-Final.dot");

    // Compute optimal repeat counts for each vertex of the marker graph.
    assembler.assembleMarkerGraphVertices(0);

    // Compute consensus sequence for marker graph edges to be used for assembly.
    assembler.assembleMarkerGraphEdges(
        0,
        assemblyOptions.Assembly.markerGraphEdgeLengthThresholdForConsensus,
        false,
        false);

    // Use the assembly graph for global assembly.
    assembler.assemble(0);
    assembler.computeAssemblyStatistics();
    assembler.writeGfa1("Assembly.gfa");
    assembler.writeFasta("Assembly.fasta");

}


// This function sets nr_overcommit_hugepages for 2MB pages
// to a little below total memory.
// If the setting needs to be modified, it acquires
// root privilege via sudo. This may result in the
// user having to enter a password.
void ChanZuckerberg::shasta::main::setupHugePages()
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
