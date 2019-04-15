// Main program for the Shasta static executable.
// The static executable provides
// basic functionality and reduced performance. 
// For full functionality use the shared library built
// under directory src.

// Shasta.
#include "Assembler.hpp"
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



// In class AssemblyOptions, we choose names that are
// consistent with options names in shasta.conf,
// even though we are violating naming conventions used
// in the rest of the Shasta code.
class ChanZuckerberg::shasta::AssemblyOptions {
public:



    class ReadsOptions {
    public:
        // useRunLengthReads; This will be phased out and always be true.
        int minReadLength;
        class PalindromicReadOptions {
        public:
            int maxSkip;
            int maxMarkerFrequency;
            double alignedFractionThreshold;
            double nearDiagonalFractionThreshold;
            int deltaThreshold;
            void write(ostream& s) const
            {
                s << "Reads.palindromicReads.maxSkip = " << maxSkip << "\n";
                s << "Reads.palindromicReads.maxMarkerFrequency = " << maxMarkerFrequency << "\n";
                s << "Reads.palindromicReads.alignedFractionThreshold = " << alignedFractionThreshold << "\n";
                s << "Reads.palindromicReads.nearDiagonalFractionThreshold = " << nearDiagonalFractionThreshold << "\n";
                s << "Reads.palindromicReads.deltaThreshold = " << deltaThreshold << "\n";
            }
        };
        PalindromicReadOptions palindromicReads;

        void write(ostream& s) const
        {
            s << "Reads.minReadLength = " << minReadLength << "\n";
            palindromicReads.write(s);
        }
    };
    ReadsOptions Reads;



    class KmersOptions {
    public:
        int k;
        double probability;
        void write(ostream& s) const
        {
            s << "Kmers.k = " << k << "\n";
            s << "Kmers.probability = " << probability << "\n";
        }
    };
    KmersOptions Kmers;



    class MinHashOptions {
    public:
        int m;
        double hashFraction;
        int minHashIterationCount;
        int maxBucketSize;
        int minFrequency;
        void write(ostream& s) const
        {
            s << "MinHash.m = " << m << "\n";
            s << "MinHash.hashFraction = " << hashFraction << "\n";
            s << "MinHash.minHashIterationCount = " << minHashIterationCount << "\n";
            s << "MinHash.maxBucketSize = " << maxBucketSize << "\n";
            s << "MinHash.minFrequency = " << minFrequency << "\n";
        }
    };
    MinHashOptions MinHash;



    class AlignOptions {
    public:
        int maxSkip;
        int maxMarkerFrequency;
        int minAlignedMarkerCount;
        int maxTrim;
        void write(ostream& s) const
        {
            s << "Align.m = " << maxSkip << "\n";
            s << "Align.maxMarkerFrequency = " << maxMarkerFrequency << "\n";
            s << "Align.minAlignedMarkerCount = " << minAlignedMarkerCount << "\n";
            s << "Align.maxTrim = " << maxTrim << "\n";
        }
    };
    AlignOptions Align;



    class ReadGraphOptions {
    public:
        int maxAlignmentCount;
        int minComponentSize;
        int maxChimericReadDistance;
        void write(ostream& s) const
        {
            s << "ReadGraph.maxAlignmentCount = " << maxAlignmentCount << "\n";
            s << "ReadGraph.minComponentSize = " << minComponentSize << "\n";
            s << "ReadGraph.maxChimericReadDistance = " << maxChimericReadDistance << "\n";
        }
    };
    ReadGraphOptions ReadGraph;



    class MarkerGraphOptions {
    public:
        int minCoverage;
        int maxCoverage;
        int lowCoverageThreshold;
        int highCoverageThreshold;
        int maxDistance;
        int edgeMarkerSkipThreshold;
        int pruneIterationCount;
        string simplifyMaxLength;
        void write(ostream& s) const
        {
            s << "MarkerGraph.minCoverage = " << minCoverage << "\n";
            s << "MarkerGraph.maxCoverage = " << maxCoverage << "\n";
            s << "MarkerGraph.lowCoverageThreshold = " << lowCoverageThreshold << "\n";
            s << "MarkerGraph.highCoverageThreshold = " << highCoverageThreshold << "\n";
            s << "MarkerGraph.maxDistance = " << maxDistance << "\n";
            s << "MarkerGraph.edgeMarkerSkipThreshold = " << edgeMarkerSkipThreshold << "\n";
            s << "MarkerGraph.pruneIterationCount = " << pruneIterationCount << "\n";
            s << "MarkerGraph.simplifyMaxLength = " << simplifyMaxLength << "\n";
        }
    };
    MarkerGraphOptions MarkerGraph;



    class AssemblyOptionsInner {
    public:
        int markerGraphEdgeLengthThresholdForConsensus;
        void write(ostream& s) const
        {
            s << "Assembly.markerGraphEdgeLengthThresholdForConsensus = " <<
                    markerGraphEdgeLengthThresholdForConsensus << "\n";
        }
    };
    AssemblyOptionsInner Assembly;


    void write(ostream& s) const
    {
        Reads.write(s);
        Kmers.write(s);
        MinHash.write(s);
        Align.write(s);
        ReadGraph.write(s);
        MarkerGraph.write(s);
        Assembly.write(s);
        s << endl;
    }
};



void ChanZuckerberg::shasta::shastaMain(int argumentCount, const char** arguments)
{
    // Some names in the boost program_options library.
    using boost::program_options::command_line_parser;
    using boost::program_options::options_description;
    using boost::program_options::value;
    using boost::program_options::variables_map;


    const string executableDescription =
        "\nThis is the static executable for the Shasta assembler. "
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

    // Options that are only allowed only on the command line.
    options_description commandLineOnlyOptions(
        "Options allowed only on the command line");
    string configFileName;
    commandLineOnlyOptions.add_options()
        ("help", 
        "Write a help message.")
        ("config", 
        boost::program_options::value<string>(&configFileName),
        "Configuration file name.")
        ;



    // Options allowed on the command line and in the config file.
    // Values specified in the command line take precedence.
    options_description options(
        "Options allowed on the command line and in the config file");
    AssemblyOptions assemblyOptions;
    vector < string > inputFastaFileNames;
    string outputDirectory;

    options.add_options()

        ("input", 
        value< vector<string> >(&inputFastaFileNames)->multitoken(),
        "Names of input FASTA files. Specify at least one.")
        
        ("output",
        value<string>(&outputDirectory)->
        default_value("ShastaRun"),
        "Name of the output directory. Must not exist.")

        ("Reads.minReadLength", 
        value<int>(&assemblyOptions.Reads.minReadLength)->
        default_value(10000),
        "Read length cutoff.")
    
        ("Reads.palindromicReads.maxSkip", 
        value<int>(&assemblyOptions.Reads.palindromicReads.maxSkip)->
        default_value(100),
        "Used for palindromic read detection.")
    
        ("Reads.palindromicReads.maxMarkerFrequency", 
        value<int>(&assemblyOptions.Reads.palindromicReads.maxMarkerFrequency)->
        default_value(10),
        "Used for palindromic read detection.")
    
        ("Reads.palindromicReads.alignedFractionThreshold", 
        value<double>(&assemblyOptions.Reads.palindromicReads.alignedFractionThreshold)->
        default_value(0.1, "0.1"),
        "Used for palindromic read detection.")
    
        ("Reads.palindromicReads.nearDiagonalFractionThreshold", 
        value<double>(&assemblyOptions.Reads.palindromicReads.nearDiagonalFractionThreshold)->
        default_value(0.1, "0.1"),
        "Used for palindromic read detection.")

        ("Reads.palindromicReads.deltaThreshold",
         value<int>(&assemblyOptions.Reads.palindromicReads.deltaThreshold)->
         default_value(100),
         "Used for palindromic read detection.")

        ("Kmers.k",
        value<int>(&assemblyOptions.Kmers.k)->
        default_value(10),
        "Length of marker k-mers (in run-length space).")

        ("Kmers.probability",
        value<double>(&assemblyOptions.Kmers.probability)->
        default_value(0.1, "0.1"),
        "Probability that a k-mer is used as a marker.")

        ("MinHash.m",
        value<int>(&assemblyOptions.MinHash.m)->
        default_value(4),
        "The number of consecutive markers that define a MinHash/LowHash feature.")

        ("MinHash.hashFraction",
        value<double>(&assemblyOptions.MinHash.hashFraction)->
        default_value(0.01, "0.01"),
        "Defines how low a hash has to be to be used with the LowHash algorithm.")

        ("MinHash.minHashIterationCount",
        value<int>(&assemblyOptions.MinHash.minHashIterationCount)->
        default_value(10),
        "The number of MinHash/LowHash iterations.")

        ("MinHash.maxBucketSize",
        value<int>(&assemblyOptions.MinHash.maxBucketSize)->
        default_value(10),
        "The maximum bucket size to be used by the MinHash/LowHash algoritm.")

        ("MinHash.minFrequency",
        value<int>(&assemblyOptions.MinHash.minFrequency)->
        default_value(2),
        "The minimum number of times a pair of reads must be found by the MinHash/LowHash algorithm "
        "in order to be considered a candidate alignment.")

        ("Align.maxSkip",
        value<int>(&assemblyOptions.Align.maxSkip)->
        default_value(30),
        "The maximum number of markers that an alignment is allowed to skip.")

        ("Align.maxMarkerFrequency",
        value<int>(&assemblyOptions.Align.maxMarkerFrequency)->
        default_value(10),
        "Marker frequency threshold.")

        ("Align.minAlignedMarkerCount",
        value<int>(&assemblyOptions.Align.minAlignedMarkerCount)->
        default_value(100),
        "The minimum number of aligned markers for an alignment to be used.")

        ("Align.maxTrim",
        value<int>(&assemblyOptions.Align.maxTrim)->
        default_value(30),
        "The maximum number of trim markers tolerated at the beginning and end of an alignment.")

        ("ReadGraph.maxAlignmentCount",
        value<int>(&assemblyOptions.ReadGraph.maxAlignmentCount)->
        default_value(6),
        "The maximum alignments to be kept for each read.")

        ("ReadGraph.minComponentSize",
        value<int>(&assemblyOptions.ReadGraph.minComponentSize)->
        default_value(100),
        "The minimum size (number of oriented reads) of a connected component "
        "of the read graph to be kept.")

        ("ReadGraph.maxChimericReadDistance",
        value<int>(&assemblyOptions.ReadGraph.maxChimericReadDistance)->
        default_value(2),
        "Used for chimeric read detection.")

        ("MarkerGraph.minCoverage",
        value<int>(&assemblyOptions.MarkerGraph.minCoverage)->
        default_value(10),
        "Minimum number of markers for a marker graph vertex.")

        ("MarkerGraph.maxCoverage",
        value<int>(&assemblyOptions.MarkerGraph.maxCoverage)->
        default_value(100),
        "Maximum number of markers for a marker graph vertex.")

        ("MarkerGraph.lowCoverageThreshold",
        value<int>(&assemblyOptions.MarkerGraph.lowCoverageThreshold)->
        default_value(0),
        "Used during approximate transitive reduction.")

        ("MarkerGraph.highCoverageThreshold",
        value<int>(&assemblyOptions.MarkerGraph.highCoverageThreshold)->
        default_value(256),
        "Used during approximate transitive reduction.")

        ("MarkerGraph.maxDistance",
        value<int>(&assemblyOptions.MarkerGraph.maxDistance)->
        default_value(30),
        "Used during approximate transitive reduction.")

        ("MarkerGraph.edgeMarkerSkipThreshold",
        value<int>(&assemblyOptions.MarkerGraph.edgeMarkerSkipThreshold)->
        default_value(100),
        "Used during approximate transitive reduction.")

        ("MarkerGraph.pruneIterationCount",
        value<int>(&assemblyOptions.MarkerGraph.pruneIterationCount)->
        default_value(6),
        "Number of prune iterations.")

        ("MarkerGraph.simplifyMaxLength",
        value<string>(&assemblyOptions.MarkerGraph.simplifyMaxLength)->
        default_value("10,100,1000"),
        "Maximum lengths (in markers) used at each iteration of simplifyMarkerGraph.")

        ("AssemblyOptions.markerGraphEdgeLengthThresholdForConsensus",
        value<int>(&assemblyOptions.Assembly.markerGraphEdgeLengthThresholdForConsensus)->
        default_value(1000),
        "Controls assembly of long marker graph edges.")

        ;
    

        
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


    // Write a startup message
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
    // We will use them below after changing directory to the outpiut directory.
    vector<string> inputFastaFileAbsolutePaths;
    for(const string& inputFastaFileName: inputFastaFileNames) {
        inputFastaFileAbsolutePaths.push_back(filesystem::getAbsolutePath(inputFastaFileName));
    }


    // If the output directory exists, stop.
    // Otherwise, create it and make it current..
    if(filesystem::exists(outputDirectory)) {
        throw runtime_error("Output directory " + outputDirectory + " already exists.");
    }
    filesystem::createDirectory(outputDirectory);
    filesystem::changeDirectory(outputDirectory);

    // Create the Data and threadLogs directories.
    filesystem::createDirectory("Data");
    filesystem::createDirectory("threadLogs");

    // Write out the option values we are using.
    cout << "Options in use:" << endl;
    cout << "Input FASTA files: ";
    copy(inputFastaFileNames.begin(), inputFastaFileNames.end(), ostream_iterator<string>(cout, " "));
    cout << endl;
    cout << "outputDirectory = " << outputDirectory << endl;
    assemblyOptions.write(cout);

    // Create the Assembler.
    Assembler assembler("Data/", 2*1024*1024, true);
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

#if 0
# Flag chimeric reads.
a.flagChimericReads(
    maxChimericReadDistance = int(config['ReadGraph']['maxChimericReadDistance']))
a.computeReadGraphConnectedComponents(
    minComponentSize = int(config['ReadGraph']['minComponentSize']))

# Create vertices of the marker graph.
a.createMarkerGraphVertices(
    maxMarkerFrequency = int(config['Align']['maxMarkerFrequency']),
    maxSkip = int(config['Align']['maxSkip']),
    minCoverage = int(config['MarkerGraph']['minCoverage']),
    maxCoverage = int(config['MarkerGraph']['maxCoverage']))
a.findMarkerGraphReverseComplementVertices()

# Create edges of the marker graph.
a.createMarkerGraphEdges()
a.findMarkerGraphReverseComplementEdges()

# Approximate transitive reduction.
a.flagMarkerGraphWeakEdges(
    lowCoverageThreshold = int(config['MarkerGraph']['lowCoverageThreshold']),
    highCoverageThreshold = int(config['MarkerGraph']['highCoverageThreshold']),
    maxDistance = int(config['MarkerGraph']['maxDistance']),
    edgeMarkerSkipThreshold = int(config['MarkerGraph']['edgeMarkerSkipThreshold'])
    )

# Prune the strong subgraph of the marker graph.
a.pruneMarkerGraphStrongSubgraph(
    iterationCount = int(config['MarkerGraph']['pruneIterationCount']))

# Simplify the marker graph to remove bubbles and superbubbles.
# The maxLength parameter controls the maximum number of markers
# for a branch to be collapsed during each iteration.
a.simplifyMarkerGraph(
    maxLength = [int(s) for s in config['MarkerGraph']['simplifyMaxLength'].split(',')],
    debug = True)

# Create the assembly graph.
a.createAssemblyGraphEdges()
a.createAssemblyGraphVertices()
a.writeAssemblyGraph("AssemblyGraph-Final.dot")

# Compute optimal repeat counts for each vertex of the marker graph.
a.assembleMarkerGraphVertices()

# Optionally compute coverage data for marker graph vertices.
storeCoverageData = ast.literal_eval(config['Assembly']['storeCoverageData'])
if storeCoverageData:
    a.computeMarkerGraphVerticesCoverageData()

# Compute consensus sequence for marker graph edges to be used for assembly.
a.assembleMarkerGraphEdges(
    markerGraphEdgeLengthThresholdForConsensus =
    int(config['Assembly']['markerGraphEdgeLengthThresholdForConsensus']),
    useMarginPhase = useMarginPhase,
    storeCoverageData = storeCoverageData)

# Use the assembly graph for global assembly.
a.assemble()

a.computeAssemblyStatistics()
a.writeGfa1('Assembly.gfa')
a.writeFasta('Assembly.fasta')
#endif
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
