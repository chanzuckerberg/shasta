#ifndef SHASTA_ASSEMBLER_OPTIONS_HPP
#define SHASTA_ASSEMBLER_OPTIONS_HPP


/*******************************************************************************

CLASSES DESCRIBING ASSEMBLER OPTIONS

There are two types of options:
- Options that can be used both on the command line and in a configuration file
  ("configurable options").
- Options that can only be used in a configuration file
  ("non-configurable options").

Configuration files are divided in sections formatted like this:

[SectionName]
optionName = optionValue

The command line syntax corresponding to the above is
--SectionName.optionName optionValue

Each SectionName corresponds to a class defined below. For example,
section [Align] corresonds to class AlignOptions.

If the option is a Boolean switch, use True or False as the optionValue.



ADDING A NEW CONFIGURABLE OPTION

1. Add the option to the class corresponding to the desired section.
2. Modify the write function to that class to also write the newly added option.
3. Modify AssemblerOptions::addCommandLineOnlyOptions to reflect the new option,
   making sure to include a default value and at least a minimal help message.
4. Document the option in shasta/docs/CommandLineOptions.html.
5. If the option requires validation add it at the appropriate place in
   shasta/srcMain/main.cpp.

For options not ready for end users, it is fine in steps 3 4 5 to use a
comment just saying "Experimental - leave at default value"
or something to that effect.



ADDING A NEW NON-CONFIGURABLE OPTION

1. Add the option to class CommandLineOnlyOptions.
2. Modify the write function to that class to also write the newly added option
3. Modify AssemblerOptions::addConfigurableOptions to reflect the new option,
   making sure to include a default value and at least a minimal help message.
4. Document the option in shasta/docs/CommandLineOptions.html.
5. If the option requires validation add it at tha appropriate place in
   shasta/srcMain/main.cpp.

For options not ready for end users, it is fine in steps 3 4 5 to use a
comment just saying "Experimental - leave at default value"
or something to that effect.

*******************************************************************************/

// Boost libraries.
#include <boost/program_options.hpp>

// Standard library.
#include "iostream.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta {
    class AlignOptions;
    class AssemblerOptions;
    class AssemblyOptions;
    class CommandLineOnlyOptions;
    class KmersOptions;
    class MarkerGraphOptions;
    class MinHashOptions;
    class Mode2AssemblyOptions;
    class PalindromicReadOptions;
    class ReadsOptions;
    class ReadGraphOptions;

    // Function to convert a bool to True or False for better
    // compatibility with Python scripts.
    string convertBoolToPythonString(bool);
}



// Options only allowed on the command line and not in the configuration file.
class shasta::CommandLineOnlyOptions {
public:
    string configName;
    vector <string> inputFileNames;
    string assemblyDirectory;
    string command;
    string memoryMode;
    string memoryBacking;
    uint32_t threadCount;
    bool suppressStdoutLog;
    string exploreAccess;
    uint16_t port;
    string alignmentsPafFile;
};


class shasta::PalindromicReadOptions {
public:
    bool skipFlagging;
    int maxSkip;
    int maxDrift;
    int maxMarkerFrequency;
    double alignedFractionThreshold;
    double nearDiagonalFractionThreshold;
    int deltaThreshold;
    void write(ostream&) const;
};


// Options in the [Reads] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Reads.".
class shasta::ReadsOptions {
public:
    uint64_t representation;    // 0 = Raw, 1=RLE
    int minReadLength;
    bool noCache;
    string desiredCoverageString;
    uint64_t desiredCoverage;
    PalindromicReadOptions palindromicReads;

    void write(ostream&) const;

    void parseDesiredCoverageString();
};



// Options in the [Kmers] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Kmers.".
class shasta::KmersOptions {
public:
    int generationMethod;
    int k;
    double probability;
    double enrichmentThreshold;
    uint64_t distanceThreshold;
    string file;
    void write(ostream&) const;
};



// Options in the [MinHash] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "MinHash.".
class shasta::MinHashOptions {
public:
    int version;
    int m;
    double hashFraction;
    int minHashIterationCount;
    double alignmentCandidatesPerRead;
    int minBucketSize;
    int maxBucketSize;
    int minFrequency;
    bool allPairs;
    void write(ostream&) const;
};



// Options in the [Align] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Align.".
class shasta::AlignOptions {
public:
    int alignMethod;
    int maxSkip;
    int maxDrift;
    int maxTrim;
    int maxMarkerFrequency;
    int minAlignedMarkerCount;
    double minAlignedFraction;
    int matchScore;
    int mismatchScore;
    int gapScore;
    double downsamplingFactor;
    int bandExtend;
    int maxBand;
    int sameChannelReadAlignmentSuppressDeltaThreshold;
    bool suppressContainments;
    uint64_t align4DeltaX;
    uint64_t align4DeltaY;
    uint64_t align4MinEntryCountPerCell;
    uint64_t align4MaxDistanceFromBoundary;
    void write(ostream&) const;
};



// Options in the [ReadGraph] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "ReadGraph.".
class shasta::ReadGraphOptions {
public:
    int creationMethod;
    int maxAlignmentCount;
    int maxChimericReadDistance;
    uint64_t strandSeparationMethod;
    int crossStrandMaxDistance;
    bool removeConflicts;
    double markerCountPercentile;
    double alignedFractionPercentile;
    double maxSkipPercentile;
    double maxDriftPercentile;
    double maxTrimPercentile;
    bool flagInconsistentAlignments;
    uint64_t flagInconsistentAlignmentsTriangleErrorThreshold;
    uint64_t flagInconsistentAlignmentsLeastSquareErrorThreshold;
    uint64_t flagInconsistentAlignmentsLeastSquareMaxDistance;
    void write(ostream& ) const;
};



// Options in the [MarkerGraph] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "MarkerGraph.".
class shasta::MarkerGraphOptions {
public:
    int minCoverage;
    int maxCoverage;
    int minCoveragePerStrand;
    uint64_t minEdgeCoverage;
    uint64_t minEdgeCoveragePerStrand;
    bool allowDuplicateMarkers;
    bool cleanupDuplicateMarkers;
    double duplicateMarkersPattern1Threshold;
    int lowCoverageThreshold;
    int highCoverageThreshold;
    int maxDistance;
    int edgeMarkerSkipThreshold;
    int pruneIterationCount;
    string simplifyMaxLength;
    double crossEdgeCoverageThreshold;
    vector<size_t> simplifyMaxLengthVector;
    bool reverseTransitiveReduction;
    double peakFinderMinAreaFraction;
    uint64_t peakFinderAreaStartIndex;
    void parseSimplifyMaxLength();
    void write(ostream&) const;
};



// Assembly options that are specific to Mode 2 assembly.
// See class AssemblyGraph2 for more information.
class shasta::Mode2AssemblyOptions {
public:

    // Threshold that defines a strong branch.
    // A branch is strong if it is supported by at least this number of
    // distinct oriented reads.
    // Weak branches are subject to removal by removeWeakBranches
    // (but at least one branch in each bubble will always be kept).
    uint64_t strongBranchThreshold;

    // Epsilon for the Bayesian model used for phasing and for bubble removal.
    // This is the probability that a read appears on the wrong branch.
    double epsilon;

    // Parameters for bubble removal.
    uint64_t minConcordantReadCountForBubbleRemoval;
    uint64_t maxDiscordantReadCountForBubbleRemoval;
    double minLogPForBubbleRemoval;
    uint64_t componentSizeThresholdForBubbleRemoval;

    // Parameters for phasing.
    uint64_t minConcordantReadCountForPhasing;
    uint64_t maxDiscordantReadCountForPhasing;
    double minLogPForPhasing;

    // Parameters for superbubble removal.
    uint64_t maxSuperbubbleSize;
    uint64_t maxSuperbubbleChunkSize;
    uint64_t maxSuperbubbleChunkPathCount;
    uint64_t superbubbleEdgeLengthThreshold;

    // Parameters to suppress output.
    bool suppressGfaOutput;
    bool suppressFastaOutput;
    bool suppressDetailedOutput;
    bool suppressPhasedOutput;
    bool suppressHaploidOutput;

    void write(ostream&) const;

};



// Options in the [Assembly] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Assembly.".
class shasta::AssemblyOptions {
public:
    uint64_t mode;
    int crossEdgeCoverageThreshold;
    int markerGraphEdgeLengthThresholdForConsensus;
    string consensusCallerString;
    string consensusCaller;
    bool storeCoverageData;
    int storeCoverageDataCsvLengthThreshold;
    bool writeReadsByAssembledSegment;
    uint64_t pruneLength;

    // Options that control detangling.
    int detangleMethod;
    uint64_t detangleDiagonalReadCountMin;
    uint64_t detangleOffDiagonalReadCountMax;
    double detangleOffDiagonalRatio;

    // Options that control iterative assembly.
    bool iterative;
    uint64_t iterativeIterationCount;
    int64_t iterativePseudoPathAlignMatchScore;
    int64_t iterativePseudoPathAlignMismatchScore;
    int64_t iterativePseudoPathAlignGapScore;
    double iterativeMismatchSquareFactor;
    double iterativeMinScore;
    uint64_t iterativeMaxAlignmentCount;
    uint64_t iterativeBridgeRemovalIterationCount;
    uint64_t iterativeBridgeRemovalMaxDistance;

    // Mode 2 assembly options.
    Mode2AssemblyOptions mode2Options;

    void write(ostream&) const;

    // If a relative path is provided for a Bayesian consensus caller
    // replace it with its absolute path.
    void parseConsensusCallerString();
};



class shasta::AssemblerOptions {
public:

    // Object containing the options.
    CommandLineOnlyOptions commandLineOnlyOptions;
    ReadsOptions readsOptions;
    KmersOptions kmersOptions;
    MinHashOptions minHashOptions;
    AlignOptions alignOptions;
    ReadGraphOptions readGraphOptions;
    MarkerGraphOptions markerGraphOptions;
    AssemblyOptions assemblyOptions;

    // Constructor.
    AssemblerOptions(int argumentCount, const char** arguments);

    // Add configurable options to the Boost option description object.
    void addCommandLineOnlyOptions();
    void addConfigurableOptions();

    // Write the options as a config file.
    void write(ostream&) const;

    // Boost program_options library objects.
    boost::program_options::options_description commandLineOnlyOptionsDescription;
    boost::program_options::options_description configurableOptionsDescription;
    boost::program_options::options_description allOptionsDescription;

    // This one is the same as allOptionsDescription, with
    // "--invalidOption" added to capture invalid positional options.
    vector<string> invalidPositionalOptions;
    boost::program_options::options_description allOptionsIncludingInvalidDescription;

};

#endif

