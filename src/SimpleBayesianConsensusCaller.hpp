#ifndef CZI_SHASTA_SIMPLE_BAYESIAN_CONSENSUS_CALLER_HPP
#define CZI_SHASTA_SIMPLE_BAYESIAN_CONSENSUS_CALLER_HPP

/*******************************************************************************

A SimpleBayesianConsensusCaller uses a simple Bayesian approach
to compute the "best" base and repeat count at a position of an alignment.

Based on initial work by Ryan Lorig-Roach at UCSC, the method works as follows.
Here, n is the true repeat count, m is the observed repeat count,
and mi the observed repeat counts in a set of reads.

- Once for a given sequencing technology, estimate conditional probabilities
P(m | n, base read) by mapping reads to portions of a reference
known not to contain variants.

- Use Bayes theorem to estimate

P(n | mi, base) proportional to P(n) times product over i P(m | n, base read)

Where P(n) is the prior probability of a homopolymer run of length n
and can be initially neglected.

Note that in the above, the base read at a given alignment position
must take into account which strand each read is on.

*******************************************************************************/

#include "ConsensusCaller.hpp"
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <array>
#include <string>
#include <limits>
#include <map>

using ChanZuckerberg::shasta::Consensus;
using std::ifstream;
using std::vector;
using std::array;
using std::string;
using std::pair;
using std::map;

namespace ChanZuckerberg {
    namespace shasta {
        class SimpleBayesianConsensusCaller;
    }
}


const double INF = std::numeric_limits<double>::infinity();;


// Given a set of observations (repeat, strand, base), predict the true repeat count
class ChanZuckerberg::shasta::SimpleBayesianConsensusCaller:
    public ChanZuckerberg::shasta::ConsensusCaller {
public:
    /// ----- Methods ----- ///

    // The constructor does not have any parameters.
    // All data is read from file SimpleBayesianConsensusCaller.csv
    // in the run directory. We will update the documentation accordingly.
    SimpleBayesianConsensusCaller();

    // Given a coverage object, return the most likely run length, and the normalized log likelihood vector for all run
    // lengths as a pair
    uint16_t predictRunlength(const Coverage &coverage, AlignedBase consensusBase, vector<double>& logLikelihoodY) const;

    AlignedBase predictConsensusBase(const Coverage& coverage) const;

    // This is the primary function of this class. Given a coverage object and consensus base, predict the true
    // run length of the aligned bases at a position
    virtual Consensus operator()(const Coverage&) const;

private:

    /// ---- Attributes ---- ///

    // The name specified under the field ">Name" in the configuration file
    string configurationName;

    // Defined at initialization, this is the size of the probability matrix generated from the data
    uint16_t maxRunlength;

    // Boolean flags for configuration
    bool ignoreNonConsensusBaseRepeats;
    bool predictGapRunlengths;
    bool countGapsAsZeros;

    // p(X|Y) normalized for each Y, where X = observed and Y = True run length
    array<vector<vector<double> >, 4> probabilityMatrices;

    // priors p(Y) normalized for each Y, where X = observed and Y = True run length
    array<vector<double>, 2> priors;

    /// ----- Methods ----- ///

    // For parsing any character separated file format
    void splitAsDouble(string s, char separatorChar, vector<double>& tokens);
    void splitAsString(string s, char separatorChar, vector<string>& tokens);

    // Read each probability matrix from its file and store them in a vector (assuming decibel units, aka base 10)
    // Each delimited table in text should be preceded by a fasta-like header e.g.: ">A" for the base it corresponds to.
    // This converts each line to a vector of doubles, appending probabilityMatrices according to the matrix header.
    void loadConfiguration(ifstream& matrixFile);

    // Parsing functions broken out for readability
    void parseName(ifstream& matrixFile, string& line);
    void parsePrior(ifstream& matrixFile, string& line, vector<string>& tokens);
    void parseLikelihood(ifstream& matrixFile, string& line, vector<string>& tokens);

    // For a given vector of likelihoods over each Y value, normalize by the maximum
    void normalizeLikelihoods(vector<double>& x, double xMax) const;

    // Count the number of times each unique repeat was observed, to reduce redundancy in calculating log likelihoods
    void factorRepeats(array<map<uint16_t,uint16_t>,2>& factoredRepeats, const Coverage& coverage) const;
    void factorRepeats(array<map<uint16_t,uint16_t>,2>& factoredRepeats, const Coverage& coverage, AlignedBase consensusBase) const;

    // For debugging or exporting
    void printPriors(char separator);
    void printProbabilityMatrices(char separator=',');
    void printLogLikelihoodVector(vector<double>& logLikelihoods);
};

void testSimpleBayesianConsensusCaller();

#endif
