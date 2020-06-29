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
2. Modify the write function to that class to also write the newly added option.3.
3. Modify AssemblerOptions::addCommandLineOnlyOptions to reflect the new option,
   making sure to include a default value and at least a minimal help message.
4. Add the option to shasta/conf/shasta.conf with a short comment
   and its default value.
5. Document the option in shasta/docs/CommandLineOptions.html.
6. If the option requires validation add it at the appropriate place in
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
    class AssemblerOptions;
}



// Class AssemblyOptions contains one nested class
// corresponding to each group of options.
class shasta::AssemblerOptions {
public:

    // Options only allowed on the command line and not in the configuration file.
    class CommandLineOnlyOptions {
    public:
        string configFileName;
        vector <string> inputFileNames;
        string assemblyDirectory;
        string command;
        string memoryMode;
        string memoryBacking;
        uint32_t threadCount;
#ifdef SHASTA_HTTP_SERVER
        string exploreAccess;
        uint16_t port;
#endif
    };
    CommandLineOnlyOptions commandLineOnlyOptions;



    // Options in the [Reads] section of the configuration file.
    // Can also be entered on the command line with option names
    // beginning with "Reads.".
    class ReadsOptions {
    public:
        int minReadLength;
        bool noCache;
        class PalindromicReadOptions {
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
        PalindromicReadOptions palindromicReads;

        void write(ostream&) const;
    };
    ReadsOptions readsOptions;



    // Options in the [Kmers] section of the configuration file.
    // Can also be entered on the command line with option names
    // beginning with "Kmers.".
    class KmersOptions {
    public:
        int generationMethod;
        int k;
        double probability;
        double enrichmentThreshold;
        string file;
        void write(ostream&) const;
    };
    KmersOptions kmersOptions;



    // Options in the [MinHash] section of the configuration file.
    // Can also be entered on the command line with option names
    // beginning with "MinHash.".
    class MinHashOptions {
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
    MinHashOptions minHashOptions;



    // Options in the [Align] section of the configuration file.
    // Can also be entered on the command line with option names
    // beginning with "Align.".
    class AlignOptions {
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
        int sameChannelReadAlignmentSuppressDeltaThreshold;
        bool suppressContainments;
        void write(ostream&) const;
    };
    AlignOptions alignOptions;



    // Options in the [ReadGraph] section of the configuration file.
    // Can also be entered on the command line with option names
    // beginning with "ReadGraph.".
    class ReadGraphOptions {
    public:
        int creationMethod;
        int maxAlignmentCount;
        int minComponentSize;
        int maxChimericReadDistance;
        int crossStrandMaxDistance;
        int containedNeighborCount;
        int uncontainedNeighborCountPerDirection;
        bool removeConflicts;
        void write(ostream& ) const;
    };
    ReadGraphOptions readGraphOptions;



    // Options in the [MarkerGraph] section of the configuration file.
    // Can also be entered on the command line with option names
    // beginning with "MarkerGraph.".
    class MarkerGraphOptions {
    public:
        int minCoverage;
        int maxCoverage;
        int minCoveragePerStrand;
        int lowCoverageThreshold;
        int highCoverageThreshold;
        int maxDistance;
        int edgeMarkerSkipThreshold;
        int pruneIterationCount;
        string simplifyMaxLength;
        double crossEdgeCoverageThreshold;
        uint64_t refineThreshold;
        vector<size_t> simplifyMaxLengthVector;
        bool reverseTransitiveReduction;
        void parseSimplifyMaxLength();
        void write(ostream&) const;
    };
    MarkerGraphOptions markerGraphOptions;



    // Options in the [Assembly] section of the configuration file.
    // Can also be entered on the command line with option names
    // beginning with "Assembly.".
    class AssemblyOptions {
    public:
        int crossEdgeCoverageThreshold;
        int markerGraphEdgeLengthThresholdForConsensus;
        string consensusCaller;
        bool storeCoverageData;
        int storeCoverageDataCsvLengthThreshold;
        bool writeReadsByAssembledSegment;
        int detangleMethod;
        void write(ostream&) const;
    };
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

    // Function to convert a bool to True or False for better
    // compatibility with Python scripts.
    static string convertBoolToPythonString(bool);

};

#endif

