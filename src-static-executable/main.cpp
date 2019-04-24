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
        void shastaMain(int argumentCount, const char** arguments);
        class AssemblyOptions;
    }
}
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/program_options.hpp>

// Standard library.
#include "fstream.hpp"
#include "iostream.hpp"
#include "iterator.hpp"
#include "stdexcept.hpp"



void ChanZuckerberg::shasta::shastaMain(int argumentCount, const char** arguments)
{
    cout << buildId() << endl;

    // Some names in the boost program_options library.
    using boost::program_options::command_line_parser;
    using boost::program_options::options_description;
    using boost::program_options::value;
    using boost::program_options::variables_map;

    const string executableDescription =
        buildId() +
        "\n\nThis is the static executable for the Shasta assembler. "
        "It provides limited Shasta functionality "
        "at reduced performance but has no dependencies and requires no installation.\n\n"
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

    // Check that we have at least one input FASTA file.     
    if (inputFastaFileNames.empty()) {
        cout << executableDescription << commandLineOptions << endl;
        throw runtime_error("Specify at least one input FASTA file.");
    }

    // Parse MarkerGraph.simplifyMaxLength
    vector<size_t> simplifyMaxLength;
    {
        boost::tokenizer< boost::char_separator<char> > tokenizer(
            assemblyOptions.MarkerGraph.simplifyMaxLength, boost::char_separator<char>(","));
        for(const string token: tokenizer) {
            try {
                size_t numberEndsHere;
                const size_t value = std::stoi(token, &numberEndsHere);
                if(numberEndsHere != token.size()) {
                    throw runtime_error("Error parsing MarkerGraph.simplifyMaxLength " +
                        assemblyOptions.MarkerGraph.simplifyMaxLength);
                }
                simplifyMaxLength.push_back(value);
            } catch(std::invalid_argument e) {
                throw runtime_error("Error parsing MarkerGraph,simplifyMaxLength " +
                    assemblyOptions.MarkerGraph.simplifyMaxLength);
            }
        }
    }



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
        "\n\nThis is the static executable for the Shasta assembler. "
        "It provides limited Shasta functionality "
        "at reduced performance but has no dependencies and requires no installation.\n\n"
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

    // Create the Data directory.
    filesystem::createDirectory("Data");

    // Write out the option values we are using.
    cout << "Options in use:" << endl;
    cout << "Input FASTA files: ";
    copy(inputFastaFileNames.begin(), inputFastaFileNames.end(), ostream_iterator<string>(cout, " "));
    cout << endl;
    cout << "outputDirectory = " << outputDirectory << endl << endl;
    assemblyOptions.write(cout);
    {
        ofstream configurationFile("shasta.conf");
        assemblyOptions.write(configurationFile);
    }

    // Create the Assembler.
    Assembler assembler("Data/", true, 2*1024*1024);
    assembler.setupConsensusCaller("SimpleConsensusCaller");

    // Add reads from the specified FASTA files.
    assembler.accessReadsReadWrite();
    assembler.accessReadNamesReadWrite();
    for(const string& inputFastaFileName: inputFastaFileAbsolutePaths) {
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
    assembler.simplifyMarkerGraph(simplifyMaxLength, false);

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

    // Final disclaimer message.
    cout << "\n" << buildId() << endl;
    cout << "This run was done using the Shasta static executable,\n"
        "which provides limited functionality at reduced performance.\n"
        "Depending on run and machine characteristics, the performance penalty\n"
        "can be large (up to a slow down of a factor of 3).\n"
        "The Shasta static executable should not be used for\n"
        "large runs or for benchmarking." << endl;

}



int main(int argumentCount, const char** arguments)
{
    try {
    
        shastaMain(argumentCount, arguments);      
          
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
