#include "AssemblyOptions.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/tokenizer.hpp>

// Standard library.
#include"stdexcept.hpp"



// Add the AssemblyOptions to a Boost option description object.
void AssemblyOptions::add(boost::program_options::options_description& options)
{
    using boost::program_options::value;

    options.add_options()

        ("Reads.minReadLength",
        value<int>(&readsOptions.minReadLength)->
        default_value(10000),
        "Read length cutoff.")

        ("Reads.palindromicReads.maxSkip",
        value<int>(&readsOptions.palindromicReads.maxSkip)->
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

        ("Kmers.k",
        value<int>(&kmersOptions.k)->
        default_value(10),
        "Length of marker k-mers (in run-length space).")

        ("Kmers.probability",
        value<double>(&kmersOptions.probability)->
        default_value(0.1, "0.1"),
        "Probability that a k-mer is used as a marker.")

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
        "The number of MinHash/LowHash iterations.")

        ("MinHash.maxBucketSize",
        value<int>(&minHashOptions.maxBucketSize)->
        default_value(10),
        "The maximum bucket size to be used by the MinHash/LowHash algoritm.")

        ("MinHash.minFrequency",
        value<int>(&minHashOptions.minFrequency)->
        default_value(2),
        "The minimum number of times a pair of reads must be found by the MinHash/LowHash algorithm "
        "in order to be considered a candidate alignment.")

        ("Align.maxSkip",
        value<int>(&alignOptions.maxSkip)->
        default_value(30),
        "The maximum number of markers that an alignment is allowed to skip.")

        ("Align.maxTrim",
        value<int>(&alignOptions.maxTrim)->
        default_value(30),
        "The maximum number of trim markers tolerated at the beginning and end of an alignment.")

        ("Align.maxMarkerFrequency",
        value<int>(&alignOptions.maxMarkerFrequency)->
        default_value(10),
        "Marker frequency threshold.")

        ("Align.minAlignedMarkerCount",
        value<int>(&alignOptions.minAlignedMarkerCount)->
        default_value(100),
        "The minimum number of aligned markers for an alignment to be used.")

        ("ReadGraph.maxAlignmentCount",
        value<int>(&readGraphOptions.maxAlignmentCount)->
        default_value(6),
        "The maximum alignments to be kept for each read.")

        ("ReadGraph.minComponentSize",
        value<int>(&readGraphOptions.minComponentSize)->
        default_value(100),
        "The minimum size (number of oriented reads) of a connected component "
        "of the read graph to be kept.")

        ("ReadGraph.maxChimericReadDistance",
        value<int>(&readGraphOptions.maxChimericReadDistance)->
        default_value(2),
        "Used for chimeric read detection.")

        ("MarkerGraph.minCoverage",
        value<int>(&markerGraphOptions.minCoverage)->
        default_value(10),
        "Minimum number of markers for a marker graph vertex.")

        ("MarkerGraph.maxCoverage",
        value<int>(&markerGraphOptions.maxCoverage)->
        default_value(100),
        "Maximum number of markers for a marker graph vertex.")

        ("MarkerGraph.lowCoverageThreshold",
        value<int>(&markerGraphOptions.lowCoverageThreshold)->
        default_value(0),
        "Used during approximate transitive reduction.")

        ("MarkerGraph.highCoverageThreshold",
        value<int>(&markerGraphOptions.highCoverageThreshold)->
        default_value(256),
        "Used during approximate transitive reduction.")

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

        ("Assembly.markerGraphEdgeLengthThresholdForConsensus",
        value<int>(&Assembly.markerGraphEdgeLengthThresholdForConsensus)->
        default_value(1000),
        "Controls assembly of long marker graph edges.")

        ("Assembly.consensusCaller",
        value<string>(&Assembly.consensusCaller)->
        default_value("SimpleConsensusCaller"),
        "Selects the consensus caller for repeat counts.\n"
        "SimpleConsensusCaller is the only choice currently\n"
        "supported by the Shasta executable.\n"
        "Other choices are available with the Shasta library.")

        ("Assembly.useMarginPhase",
        value<string>(&Assembly.useMarginPhase)->
        default_value("False"),
        "Used to turn on margin phase.")

        ("Assembly.storeCoverageData",
        value<string>(&Assembly.storeCoverageData)->
        default_value("False"),
        "Used to request storing coverage data.")
        ;
}



void AssemblyOptions::ReadsOptions::PalindromicReadOptions::write(ostream& s) const
{
    s << "palindromicReads.maxSkip = " << maxSkip << "\n";
    s << "palindromicReads.maxMarkerFrequency = " << maxMarkerFrequency << "\n";
    s << "palindromicReads.alignedFractionThreshold = " << alignedFractionThreshold << "\n";
    s << "palindromicReads.nearDiagonalFractionThreshold = " << nearDiagonalFractionThreshold << "\n";
    s << "palindromicReads.deltaThreshold = " << deltaThreshold << "\n";
}



void AssemblyOptions::ReadsOptions::write(ostream& s) const
{
    s << "[Reads]\n";
    s << "minReadLength = " << minReadLength << "\n";
    palindromicReads.write(s);
}



void AssemblyOptions::KmersOptions::write(ostream& s) const
{
    s << "[Kmers]\n";
    s << "k = " << k << "\n";
    s << "probability = " << probability << "\n";
}



void AssemblyOptions::MinHashOptions::write(ostream& s) const
{
    s << "[MinHash]\n";
    s << "m = " << m << "\n";
    s << "hashFraction = " << hashFraction << "\n";
    s << "minHashIterationCount = " << minHashIterationCount << "\n";
    s << "maxBucketSize = " << maxBucketSize << "\n";
    s << "minFrequency = " << minFrequency << "\n";
}



void AssemblyOptions::AlignOptions::write(ostream& s) const
{
    s << "[Align]\n";
    s << "maxSkip = " << maxSkip << "\n";
    s << "maxTrim = " << maxTrim << "\n";
    s << "maxMarkerFrequency = " << maxMarkerFrequency << "\n";
    s << "minAlignedMarkerCount = " << minAlignedMarkerCount << "\n";
}



void AssemblyOptions::ReadGraphOptions::write(ostream& s) const
{
    s << "[ReadGraph]\n";
    s << "maxAlignmentCount = " << maxAlignmentCount << "\n";
    s << "minComponentSize = " << minComponentSize << "\n";
    s << "maxChimericReadDistance = " << maxChimericReadDistance << "\n";
}


void AssemblyOptions::MarkerGraphOptions::write(ostream& s) const
{
    s << "[MarkerGraph]\n";
    s << "minCoverage = " << minCoverage << "\n";
    s << "maxCoverage = " << maxCoverage << "\n";
    s << "lowCoverageThreshold = " << lowCoverageThreshold << "\n";
    s << "highCoverageThreshold = " << highCoverageThreshold << "\n";
    s << "maxDistance = " << maxDistance << "\n";
    s << "edgeMarkerSkipThreshold = " << edgeMarkerSkipThreshold << "\n";
    s << "pruneIterationCount = " << pruneIterationCount << "\n";
    s << "simplifyMaxLength = " << simplifyMaxLength << "\n";
}



void AssemblyOptions::AssemblyOptionsInner::write(ostream& s) const
{
    s << "[Assembly]\n";
    s << "markerGraphEdgeLengthThresholdForConsensus = " <<
        markerGraphEdgeLengthThresholdForConsensus << "\n";
    s << "consensusCaller = " <<
        consensusCaller << "\n";
    s << "useMarginPhase = " <<
        useMarginPhase << "\n";
    s << "storeCoverageData = " <<
        storeCoverageData << "\n";
}



void AssemblyOptions::write(ostream& s) const
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
    Assembly.write(s);
    s << endl;
}



void AssemblyOptions::MarkerGraphOptions::parseSimplifyMaxLength()
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



