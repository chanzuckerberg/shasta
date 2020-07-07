#include "AssemblerOptions.hpp"
#include "buildId.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/tokenizer.hpp>

// Standard library.
#include "fstream.hpp"
#include"stdexcept.hpp"



// Constructor.
AssemblerOptions::AssemblerOptions(int argumentCount, const char** arguments) :
    commandLineOnlyOptionsDescription("Options allowed only on the command line"),
    configurableOptionsDescription("Options allowed on the command line and in the config file")
{
    using boost::program_options::positional_options_description;
    using boost::program_options::value;
    using boost::program_options::variables_map;
    using boost::program_options::command_line_parser;

    addCommandLineOnlyOptions();
    addConfigurableOptions();

    // Define allOptionsDescription as the union of
    // commandLineOnlyOptionsDescription and configurableOptionsDescription.
    allOptionsDescription.add(commandLineOnlyOptionsDescription);
    allOptionsDescription.add(configurableOptionsDescription);

    // Machinery to capture invalid positional options.
    allOptionsIncludingInvalidDescription.add(allOptionsDescription);
    allOptionsIncludingInvalidDescription.add_options()
        ("invalidOption",
            value< vector<string> >(&invalidPositionalOptions));
    positional_options_description positionalOptionsDescription;
    positionalOptionsDescription.add("invalidOption", -1);


    // Get options from the command line.
    // These take precedence over values entered in the config file.
    variables_map variablesMap;
    store(command_line_parser(argumentCount, arguments).
          options(allOptionsIncludingInvalidDescription).
          positional(positionalOptionsDescription).
          run(),
          variablesMap);
    notify(variablesMap);

    // If any invalid positional options were specified,
    // stop here.
    if(!invalidPositionalOptions.empty()) {
        string message = "Positional options are not allowed. "
            "The following positional options were used:";
        for(const string& s: invalidPositionalOptions) {
            message += " ";
            message += s;
        }
        throw runtime_error(message);
    }

    // If help was requested, write the help message and exit.
    if (variablesMap.count("help")) {
        cout << allOptionsDescription << endl;
        ::exit(0);
    }

    // If version was requested, write the build id and exit.
    if (variablesMap.count("version")) {
        cout << buildId() << endl;
#ifdef __linux__
        cout << "Linux version" << endl;
#else
        cout << "MacOS version" << endl;
#endif
        ::exit(0);
    }

    // Get options from the config file, if one was specified.
    if(!commandLineOnlyOptions.configFileName.empty()) {
        ifstream configFile(commandLineOnlyOptions.configFileName);
        if (!configFile) {
            throw runtime_error("Unable to open open config file " +
                commandLineOnlyOptions.configFileName);
        }
        store(parse_config_file(configFile, configurableOptionsDescription), variablesMap);
        notify(variablesMap);
    }

    // Parse MarkerGraph.simplifyMaxLength.
    markerGraphOptions.parseSimplifyMaxLength();

}



// Add non-configurable options to the Boost option description object.
// These are options that can only be used on the command line
// (not in the configuration file).
void AssemblerOptions::addCommandLineOnlyOptions()
{
    using boost::program_options::value;
    using boost::program_options::bool_switch;

    commandLineOnlyOptionsDescription.add_options()

        ("help,h",
        "Write a help message.")

        ("version,v",
        "Identify the Shasta version.")

        ("config",
        value<string>(&commandLineOnlyOptions.configFileName),
        "Configuration file name.")

        ("input",
        value< vector<string> >(&commandLineOnlyOptions.inputFileNames)->multitoken(),
        "Names of input files containing reads. Specify at least one.")

        ("assemblyDirectory",
        value<string>(&commandLineOnlyOptions.assemblyDirectory)->
        default_value("ShastaRun"),
        "Name of the output directory. If command is assemble, this directory must not exist.")

        ("command",
        value<string>(&commandLineOnlyOptions.command)->
        default_value("assemble"),
        "Command to run. Must be one of: "
        "assemble, saveBinaryData, cleanupBinaryData, explore, createBashCompletionScript")

#ifdef __linux__
        ("memoryMode",
        value<string>(&commandLineOnlyOptions.memoryMode)->
        default_value("anonymous"),
        "Specify whether allocated memory is anonymous or backed by a filesystem. "
        "Allowed values: anonymous, filesystem.")

        ("memoryBacking",
        value<string>(&commandLineOnlyOptions.memoryBacking)->
        default_value("4K"),
        "Specify the type of pages used to back memory.\n"
        "Allowed values: disk, 4K , 2M (for best performance). "
        "All combinations (memoryMode, memoryBacking) are allowed "
        "except for (anonymous, disk).\n"
        "Some combinations require root privilege, which is obtained using sudo "
        "and may result in a password prompting depending on your sudo set up.")
#endif

        ("threads",
        value<uint32_t>(&commandLineOnlyOptions.threadCount)->
        default_value(0),
        "Number of threads, or 0 to use one thread per virtual processor.")
        
#ifdef SHASTA_HTTP_SERVER
        ("exploreAccess",
        value<string>(&commandLineOnlyOptions.exploreAccess)->
        default_value("user"),
        "Specify allowed access for --command explore. "
        "Allowed values: user, local, unrestricted. "
        "DO NOT CHANGE FROM DEFAULT VALUE WITHOUT UNDERSTANDING THE "
        "SECURITY IMPLICATIONS."
        )

        ("port",
        value<uint16_t>(&commandLineOnlyOptions.port)->
        default_value(17100),
        "Port to be used by the http server (command --explore).")
#endif
        ;


    // For reasons not completely understood and that there was no time to investigate,
    // the only combination that works on MacOS is
    // "--memoryMode filesystem --memoryBacking disk".
    // This incurs a performance price but this is not too much of a big deal
    // as macOS  is only to be used for small test runs.
#ifndef __linux__
    commandLineOnlyOptions.memoryMode = "filesystem";
    commandLineOnlyOptions.memoryBacking = "disk";
#endif
}



// Add configurable options to the Boost option description object.
// These are options that can be used on the command line
// and in the configuration file.
// If used in both places, the value given on the command line
// takes precedence.
void AssemblerOptions::addConfigurableOptions()
{
    using boost::program_options::value;
    using boost::program_options::bool_switch;

    configurableOptionsDescription.add_options()

        ("Reads.minReadLength",
        value<int>(&readsOptions.minReadLength)->
        default_value(10000),
        "Read length cutoff. Shorter reads are discarded.")

        ("Reads.noCache",
        bool_switch(&readsOptions.noCache)->
        default_value(false),
        "If set, skip the Linux cache when loading reads. "
        "This is done by specifying the O_DIRECT flag when opening "
        "input files containing reads.")

        ("Reads.palindromicReads.skipFlagging",
        bool_switch(&readsOptions.palindromicReads.skipFlagging)->
        default_value(false),
        "Skip flagging palindromic reads. Oxford Nanopore reads should be flagged for better results.")

        ("Reads.palindromicReads.maxSkip",
        value<int>(&readsOptions.palindromicReads.maxSkip)->
        default_value(100),
        "Used for palindromic read detection.")

        ("Reads.palindromicReads.maxDrift",
        value<int>(&readsOptions.palindromicReads.maxDrift)->
        default_value(100),
        "Used for palindromic read detection.")

        ("Reads.palindromicReads.maxMarkerFrequency",
        value<int>(&readsOptions.palindromicReads.maxMarkerFrequency)->
        default_value(10),
        "Used for palindromic read detection.")

        ("Reads.palindromicReads.alignedFractionThreshold",
        value<double>(&readsOptions.palindromicReads.alignedFractionThreshold)->
        default_value(0.1, "0.1"),
        "Used for palindromic read detection.")

        ("Reads.palindromicReads.nearDiagonalFractionThreshold",
        value<double>(&readsOptions.palindromicReads.nearDiagonalFractionThreshold)->
        default_value(0.1, "0.1"),
        "Used for palindromic read detection.")

        ("Reads.palindromicReads.deltaThreshold",
         value<int>(&readsOptions.palindromicReads.deltaThreshold)->
         default_value(100),
         "Used for palindromic read detection.")

         ("Kmers.generationMethod",
         value<int>(&kmersOptions.generationMethod)->
         default_value(0),
         "Method to generate marker k-mers: "
         "0 = random, "
         "1 = random, excluding globally overenriched,"
         "2 = random, excluding overenriched even in a single read,"
         "3 = read from file.")

         ("Kmers.k",
         value<int>(&kmersOptions.k)->
         default_value(10),
         "Length of marker k-mers (in run-length space).")

         ("Kmers.probability",
         value<double>(&kmersOptions.probability)->
         default_value(0.1, "0.1"),
         "Fraction k-mers used as a marker.")

        ("Kmers.enrichmentThreshold",
        value<double>(&kmersOptions.enrichmentThreshold)->
        default_value(100., "100."),
        "Enrichment threshold for Kmers.generationMethod 1 and 2.")

        ("Kmers.file",
        value<string>(&kmersOptions.file),
        "The absolute path of a file containing the k-mers "
        "to be used as markers, one per line. "
        "A relative path is not accepted. "
        "Only used if Kmers.generationMethod is 3.")

        ("MinHash.version",
        value<int>(&minHashOptions.version)->
        default_value(0),
        "Controls the version of the LowHash algorithm to use. Can be 0 (default) "
        "or 1.(experimental).")

        ("MinHash.m",
        value<int>(&minHashOptions.m)->
        default_value(4),
        "The number of consecutive markers that define a MinHash/LowHash feature.")

        ("MinHash.hashFraction",
        value<double>(&minHashOptions.hashFraction)->
        default_value(0.01, "0.01"),
        "Defines how low a hash has to be to be used with the LowHash algorithm.")

        ("MinHash.minHashIterationCount",
        value<int>(&minHashOptions.minHashIterationCount)->
        default_value(10),
        "The number of MinHash/LowHash iterations, or 0 to let "
        "--MinHash.alignmentCandidatesPerRead control the number of iterations.")

        ("MinHash.alignmentCandidatesPerRead",
        value<double>(&minHashOptions.alignmentCandidatesPerRead)->
        default_value(20.),
        "If --MinHash.minHashIterationCount is 0, MinHash iteration is stopped "
        "when the average number of alignment candidates that each read is involved in "
        "reaches this value. If --MinHash.minHashIterationCount is not 0, "
        "this is not used.")

        ("MinHash.minBucketSize",
        value<int>(&minHashOptions.minBucketSize)->
        default_value(0),
        "The minimum bucket size to be used by the LowHash algoritm.")

        ("MinHash.maxBucketSize",
        value<int>(&minHashOptions.maxBucketSize)->
        default_value(10),
        "The maximum bucket size to be used by the LowHash algoritm.")

        ("MinHash.minFrequency",
        value<int>(&minHashOptions.minFrequency)->
        default_value(2),
        "The minimum number of times a pair of reads must be found by the MinHash/LowHash algorithm "
        "in order to be considered a candidate alignment.")

        ("MinHash.allPairs",
        bool_switch(&minHashOptions.allPairs)->
        default_value(false),
        "Skip the MinHash algorithm and mark all pairs of reads as alignment"
        "candidates with both orientation. This should only be used for experimentation "
        "on very small runs because it is very time consuming.")

        ("Align.alignMethod",
        value<int>(&alignOptions.alignMethod)->
        default_value(3),
        "The alignment method to be used to create the read graph & the marker graph. "
        "0 = old Shasta method, 1 = SeqAn (slow), 3 = banded SeqAn.")

        ("Align.maxSkip",
        value<int>(&alignOptions.maxSkip)->
        default_value(30),
        "The maximum number of markers that an alignment is allowed to skip.")

        ("Align.maxDrift",
        value<int>(&alignOptions.maxDrift)->
        default_value(30),
        "The maximum amount of marker drift that an alignment is allowed to tolerate "
        "between successive markers.")

        ("Align.maxTrim",
        value<int>(&alignOptions.maxTrim)->
        default_value(30),
        "The maximum number of unaligned markers tolerated at the beginning and end of an alignment.")

        ("Align.maxMarkerFrequency",
        value<int>(&alignOptions.maxMarkerFrequency)->
        default_value(10),
        "Marker frequency threshold. Markers more frequent than this value in either of "
        "two oriented reads being aligned are discarded and not used to compute "
        "the alignment.")

        ("Align.minAlignedMarkerCount",
        value<int>(&alignOptions.minAlignedMarkerCount)->
        default_value(100),
        "The minimum number of aligned markers for an alignment to be used.")

        ("Align.minAlignedFraction",
        value<double>(&alignOptions.minAlignedFraction)->
        default_value(0.),
        "The minimum fraction of aligned markers for an alignment to be used.")

        ("Align.matchScore",
        value<int>(&alignOptions.matchScore)->
        default_value(6),
        "Match score for marker alignments (only used for alignment methods 1 and 3).")

        ("Align.mismatchScore",
        value<int>(&alignOptions.mismatchScore)->
        default_value(-1),
        "Mismatch score for marker alignments (only used for alignment methods 1 and 3).")

        ("Align.gapScore",
        value<int>(&alignOptions.gapScore)->
        default_value(-1),
        "Gap score for marker alignments (only used for alignment methods 1 and 3).")

        ("Align.downsamplingFactor",
        value<double>(&alignOptions.downsamplingFactor)->
        default_value(0.1),
        "Downsampling factor (only used for alignment method 3).")

        ("Align.bandExtend",
        value<int>(&alignOptions.bandExtend)->
        default_value(10),
        "Amount to extend the downsampled band "
        "(only used for alignment method 3).")

        ("Align.sameChannelReadAlignment.suppressDeltaThreshold",
        value<int>(&alignOptions.sameChannelReadAlignmentSuppressDeltaThreshold)->
        default_value(0),
        "If not zero, alignments between reads from the same nanopore channel "
        "and close in time are suppressed. The \"read\" meta data fields "
        "from the FASTA or FASTQ header are checked. If their difference, in "
        "absolute value, is less than the value of this option, the alignment "
        "is suppressed. This can help avoid assembly artifact. "
        "This check is only done if the two reads have identical meta data fields "
        "\"runid\", \"sampleid\", and \"ch\". "
        "If any of these meta data fields are missing, this check is suppressed and this "
        "option has no effect.")

        ("Align.suppressContainments",
        bool_switch(&alignOptions.suppressContainments)->
        default_value(false),
        "Suppress containment alignments, that is alignments in which "
        "one read is entirely contained in another read, "
        "except possibly for up to maxTrim markers at the beginning and end.")

        ("ReadGraph.creationMethod",
        value<int>(&readGraphOptions.creationMethod)->
        default_value(0),
        "The method used to create the read graph (0 = undirected, default, 1 = directed, experimental).")

        ("ReadGraph.maxAlignmentCount",
        value<int>(&readGraphOptions.maxAlignmentCount)->
        default_value(6),
        "The maximum number of alignments to be kept for each read.")

        ("ReadGraph.minComponentSize",
        value<int>(&readGraphOptions.minComponentSize)->
        default_value(100),
        "The minimum size (number of oriented reads) of a connected component "
        "of the read graph to be kept. This is currently ignored.")

        ("ReadGraph.maxChimericReadDistance",
        value<int>(&readGraphOptions.maxChimericReadDistance)->
        default_value(2),
        "Used for chimeric read detection.")

        ("ReadGraph.crossStrandMaxDistance",
        value<int>(&readGraphOptions.crossStrandMaxDistance)->
        default_value(6),
        "Maximum distance (edges) for flagCrossStrandReadGraphEdges. "
        "Set this to zero to entirely suppress flagCrossStrandReadGraphEdges.")

        ("ReadGraph.containedNeighborCount",
        value<int>(&readGraphOptions.containedNeighborCount)->
        default_value(6),
        "Maximum number of alignments to be kept for each contained read "
        "(only used when creationMethod is 1).")

        ("ReadGraph.uncontainedNeighborCountPerDirection",
        value<int>(&readGraphOptions.uncontainedNeighborCountPerDirection)->
        default_value(3),
        "Maximum number of alignments to be kept in each direction "
        "(forward, backward) for each uncontained read (only used when creationMethod is 1).")

        ("ReadGraph.removeConflicts",
        bool_switch(&readGraphOptions.removeConflicts)->
        default_value(false),
        "Remove conflicts from the read graph. Experimental - do not use.")

        ("MarkerGraph.minCoverage",
        value<int>(&markerGraphOptions.minCoverage)->
        default_value(10),
        "Minimum coverage (number of supporting oriented reads) "
        "for a marker graph vertex to be created."
        "Specifying 0 causes a suitable value of this parameter "
        "to be selected automatically.")

        ("MarkerGraph.maxCoverage",
        value<int>(&markerGraphOptions.maxCoverage)->
        default_value(100),
        "Maximum coverage (number of supporting oriented reads) "
        "for a marker graph vertex.")

        ("MarkerGraph.minCoveragePerStrand",
        value<int>(&markerGraphOptions.minCoveragePerStrand)->
        default_value(0),
        "Minimum coverage (number of supporting oriented reads) "
        "for each strand for a marker graph vertex.")

        ("MarkerGraph.lowCoverageThreshold",
        value<int>(&markerGraphOptions.lowCoverageThreshold)->
        default_value(0),
        "Used during approximate transitive reduction. Marker graph edges with coverage "
        "lower than this value are always marked as removed regardless of reachability.")

        ("MarkerGraph.highCoverageThreshold",
        value<int>(&markerGraphOptions.highCoverageThreshold)->
        default_value(256),
        "Used during approximate transitive reduction. Marker graph edges with coverage "
        "higher than this value are never marked as removed regardless of reachability.")

        ("MarkerGraph.maxDistance",
        value<int>(&markerGraphOptions.maxDistance)->
        default_value(30),
        "Used during approximate transitive reduction.")

        ("MarkerGraph.edgeMarkerSkipThreshold",
        value<int>(&markerGraphOptions.edgeMarkerSkipThreshold)->
        default_value(100),
        "Used during approximate transitive reduction.")

        ("MarkerGraph.pruneIterationCount",
        value<int>(&markerGraphOptions.pruneIterationCount)->
        default_value(6),
        "Number of prune iterations.")

        ("MarkerGraph.simplifyMaxLength",
        value<string>(&markerGraphOptions.simplifyMaxLength)->
        default_value("10,100,1000"),
        "Maximum lengths (in markers) used at each iteration of simplifyMarkerGraph.")

        ("MarkerGraph.crossEdgeCoverageThreshold",
        value<double>(&markerGraphOptions.crossEdgeCoverageThreshold)->
        default_value(0.),
        "Experimental. Cross edge coverage threshold. If this is not zero, assembly graph cross-edges "
        "with average edge coverage less than this value are removed, together with the "
        "corresponding marker graph edges. A cross edge is defined as an edge v0->v1 "
        "with out-degree(v0)>1, in-degree(v1)>1.")

        ("MarkerGraph.refineThreshold",
        value<uint64_t>(&markerGraphOptions.refineThreshold)->
        default_value(0),
        "Experimental. Length threshold, in markers, for the marker graph refinement step, "
        "or 0 to turn off the refinement step.")

        ("MarkerGraph.reverseTransitiveReduction",
        bool_switch(&markerGraphOptions.reverseTransitiveReduction)->
        default_value(false),
        "Perform approximate reverse transitive reduction of the marker graph.")

        ("Assembly.crossEdgeCoverageThreshold",
        value<int>(&assemblyOptions.crossEdgeCoverageThreshold)->
        default_value(3),
        "Maximum average edge coverage for a cross edge "
        "of the assembly graph to be removed.")

        ("Assembly.markerGraphEdgeLengthThresholdForConsensus",
        value<int>(&assemblyOptions.markerGraphEdgeLengthThresholdForConsensus)->
        default_value(1000),
        "Controls assembly of long marker graph edges.")

        ("Assembly.consensusCaller",
        value<string>(&assemblyOptions.consensusCaller)->
        default_value("Bayesian:guppy-2.3.5-a"),
        "Selects the consensus caller for repeat counts. "
        "See the documentation for available choices.")

        ("Assembly.storeCoverageData",
        bool_switch(&assemblyOptions.storeCoverageData)->
        default_value(false),
        "Used to request storing coverage data in binary format.")

        ("Assembly.storeCoverageDataCsvLengthThreshold",
        value<int>(&assemblyOptions.storeCoverageDataCsvLengthThreshold)->
        default_value(0),
        "Used to specify the minimum length of an assembled segment "
        "for which coverage data in csv format should be stored. "
        "If 0, no coverage data in csv format is stored.")

        ("Assembly.writeReadsByAssembledSegment",
        bool_switch(&assemblyOptions.writeReadsByAssembledSegment)->
        default_value(false),
        "Used to request writing the reads that contributed to assembling each segment.")

        ("Assembly.detangleMethod",
        value<int>(&assemblyOptions.detangleMethod)->
        default_value(0),
        "Experimental. Specify the method used to detangle the assembly graph. "
        "0 = no detangling, 1 = basic detangling.")
        ;
}



void AssemblerOptions::ReadsOptions::PalindromicReadOptions::write(ostream& s) const
{
    s << "palindromicReads.skipFlagging = " << convertBoolToPythonString(skipFlagging) << "\n";
    s << "palindromicReads.maxSkip = " << maxSkip << "\n";
    s << "palindromicReads.maxDrift = " << maxDrift << "\n";
    s << "palindromicReads.maxMarkerFrequency = " << maxMarkerFrequency << "\n";
    s << "palindromicReads.alignedFractionThreshold = " << alignedFractionThreshold << "\n";
    s << "palindromicReads.nearDiagonalFractionThreshold = " << nearDiagonalFractionThreshold << "\n";
    s << "palindromicReads.deltaThreshold = " << deltaThreshold << "\n";
}



void AssemblerOptions::ReadsOptions::write(ostream& s) const
{
    s << "[Reads]\n";
    s << "minReadLength = " << minReadLength << "\n";
    s << "noCache = " <<
        convertBoolToPythonString(noCache) << "\n";
    palindromicReads.write(s);
}



void AssemblerOptions::KmersOptions::write(ostream& s) const
{
    s << "[Kmers]\n";
    s << "generationMethod = " << generationMethod << "\n";
    s << "k = " << k << "\n";
    s << "probability = " << probability << "\n";
    s << "enrichmentThreshold = " << enrichmentThreshold << "\n";
    s << "file = " << file << "\n";
}



void AssemblerOptions::MinHashOptions::write(ostream& s) const
{
    s << "[MinHash]\n";
    s << "version = " << version << "\n";
    s << "m = " << m << "\n";
    s << "hashFraction = " << hashFraction << "\n";
    s << "minHashIterationCount = " << minHashIterationCount << "\n";
    s << "alignmentCandidatesPerRead = " << alignmentCandidatesPerRead << "\n";
    s << "minBucketSize = " << minBucketSize << "\n";
    s << "maxBucketSize = " << maxBucketSize << "\n";
    s << "minFrequency = " << minFrequency << "\n";
    s << "allPairs = " <<
        convertBoolToPythonString(allPairs) << "\n";
}



void AssemblerOptions::AlignOptions::write(ostream& s) const
{
    s << "[Align]\n";
    s << "alignMethod = " << alignMethod << "\n";
    s << "maxSkip = " << maxSkip << "\n";
    s << "maxDrift = " << maxDrift << "\n";
    s << "maxTrim = " << maxTrim << "\n";
    s << "maxMarkerFrequency = " << maxMarkerFrequency << "\n";
    s << "minAlignedMarkerCount = " << minAlignedMarkerCount << "\n";
    s << "minAlignedFraction = " << minAlignedFraction << "\n";
    s << "matchScore = " << matchScore << "\n";
    s << "mismatchScore = " << mismatchScore << "\n";
    s << "gapScore = " << gapScore << "\n";
    s << "downsamplingFactor = " << downsamplingFactor << "\n";
    s << "bandExtend = " << bandExtend << "\n";
    s << "sameChannelReadAlignment.suppressDeltaThreshold = " <<
        sameChannelReadAlignmentSuppressDeltaThreshold << "\n";
    s << "suppressContainments = " <<
        convertBoolToPythonString(suppressContainments) << "\n";
}



void AssemblerOptions::ReadGraphOptions::write(ostream& s) const
{
    s << "[ReadGraph]\n";
    s << "creationMethod = " << creationMethod << "\n";
    s << "maxAlignmentCount = " << maxAlignmentCount << "\n";
    s << "minComponentSize = " << minComponentSize << "\n";
    s << "maxChimericReadDistance = " << maxChimericReadDistance << "\n";
    s << "crossStrandMaxDistance = " << crossStrandMaxDistance << "\n";
    s << "containedNeighborCount = " << containedNeighborCount << "\n";
    s << "uncontainedNeighborCountPerDirection = " << uncontainedNeighborCountPerDirection << "\n";
    s << "removeConflicts = " <<
        convertBoolToPythonString(removeConflicts) << "\n";

}



void AssemblerOptions::MarkerGraphOptions::write(ostream& s) const
{
    s << "[MarkerGraph]\n";
    s << "minCoverage = " << minCoverage << "\n";
    s << "maxCoverage = " << maxCoverage << "\n";
    s << "minCoveragePerStrand = " << minCoveragePerStrand << "\n";
    s << "lowCoverageThreshold = " << lowCoverageThreshold << "\n";
    s << "highCoverageThreshold = " << highCoverageThreshold << "\n";
    s << "maxDistance = " << maxDistance << "\n";
    s << "edgeMarkerSkipThreshold = " << edgeMarkerSkipThreshold << "\n";
    s << "pruneIterationCount = " << pruneIterationCount << "\n";
    s << "simplifyMaxLength = " << simplifyMaxLength << "\n";
    s << "crossEdgeCoverageThreshold = " << crossEdgeCoverageThreshold << "\n";
    s << "refineThreshold = " << refineThreshold << "\n";
    s << "reverseTransitiveReduction = " <<
        convertBoolToPythonString(reverseTransitiveReduction) << "\n";
}



void AssemblerOptions::AssemblyOptions::write(ostream& s) const
{
    s << "[Assembly]\n";
    s << "crossEdgeCoverageThreshold = " << crossEdgeCoverageThreshold << "\n";
    s << "markerGraphEdgeLengthThresholdForConsensus = " <<
        markerGraphEdgeLengthThresholdForConsensus << "\n";
    s << "consensusCaller = " <<
        consensusCaller << "\n";
    s << "storeCoverageData = " <<
        convertBoolToPythonString(storeCoverageData) << "\n";
    s << "storeCoverageDataCsvLengthThreshold = " <<
        storeCoverageDataCsvLengthThreshold << "\n";
    s << "writeReadsByAssembledSegment = " <<
        convertBoolToPythonString(writeReadsByAssembledSegment) << "\n";
    s << "detangleMethod = " << detangleMethod << "\n";
}



void AssemblerOptions::write(ostream& s) const
{
    readsOptions.write(s);
    s << "\n";
    kmersOptions.write(s);
    s << "\n";
    minHashOptions.write(s);
    s << "\n";
    alignOptions.write(s);
    s << "\n";
    readGraphOptions.write(s);
    s << "\n";
    markerGraphOptions.write(s);
    s << "\n";
    assemblyOptions.write(s);
    s << endl;
}



void AssemblerOptions::MarkerGraphOptions::parseSimplifyMaxLength()
{
    simplifyMaxLengthVector.clear();

    boost::tokenizer< boost::char_separator<char> > tokenizer(
        simplifyMaxLength, boost::char_separator<char>(","));
    for(const string token: tokenizer) {
        try {
            size_t numberEndsHere;
            const size_t value = std::stoi(token, &numberEndsHere);
            if(numberEndsHere != token.size()) {
                throw runtime_error("Error parsing MarkerGraph.simplifyMaxLength " +
                    simplifyMaxLength);
            }
            simplifyMaxLengthVector.push_back(value);
        } catch(const std::invalid_argument& e) {
            throw runtime_error("Error parsing MarkerGraph,simplifyMaxLength " +
                simplifyMaxLength);
        }
    }

}



// Function to convert a bool to True or False for better
// compatibility with Python scripts.
string AssemblerOptions::convertBoolToPythonString(bool flag)
{
    return flag ? "True" : "False";
}

