#ifndef SHASTA_ASSEMBLY_OPTIONS_HPP
#define SHASTA_ASSEMBLY_OPTIONS_HPP

// Boost libraries.
#include <boost/program_options.hpp>

// Standard library.
#include "iostream.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace shasta {
    class AssemblyOptions;
}



// Class AssemblyOptions contains one nested class
// corresponding to each group of options.
class shasta::AssemblyOptions {
public:

    class ReadsOptions {
    public:
        int minReadLength;
        class PalindromicReadOptions {
        public:
            int maxSkip;
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



    class KmersOptions {
    public:
        int k;
        double probability;
        void write(ostream&) const;
    };
    KmersOptions kmersOptions;



    class MinHashOptions {
    public:
        int m;
        double hashFraction;
        int minHashIterationCount;
        int maxBucketSize;
        int minFrequency;
        void write(ostream&) const;
    };
    MinHashOptions minHashOptions;



    class AlignOptions {
    public:
        int maxSkip;
        int maxTrim;
        int maxMarkerFrequency;
        int minAlignedMarkerCount;
        void write(ostream&) const;
    };
    AlignOptions alignOptions;



    class ReadGraphOptions {
    public:
        int maxAlignmentCount;
        int minComponentSize;
        int maxChimericReadDistance;
        void write(ostream& ) const;
    };
    ReadGraphOptions readGraphOptions;



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
        vector<size_t> simplifyMaxLengthVector;
        void parseSimplifyMaxLength();
        void write(ostream&) const;
    };
    MarkerGraphOptions markerGraphOptions;



    class AssemblyOptionsInner {
    public:
        int markerGraphEdgeLengthThresholdForConsensus;
        string consensusCaller;
        string useMarginPhase;      // False or True
        string storeCoverageData;   // False or True
        void write(ostream&) const;
    };
    AssemblyOptionsInner Assembly;

    // Add these options to a Bost option description object.
    // Note that this cannot be const as the Boost option descriptions
    // stores pointers to members that will be filled in with option values.
    void add(boost::program_options::options_description&);

    // Write the options as a config file.
    void write(ostream&) const;
};

#endif

