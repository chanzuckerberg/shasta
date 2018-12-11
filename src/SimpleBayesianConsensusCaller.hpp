#ifndef CZI_SHASTA_SIMPLE_BAYESIAN_CONSENSUS_CALLER_HPP
#define CZI_SHASTA_SIMPLE_BAYESIAN_CONSENSUS_CALLER_HPP

/*******************************************************************************

A SimpleBayesianConsensusCaller uses a simple Bayesian approach
to compute the "best" base and repeat count at a position of an alignment.

Based on initial work by Ryan Lorigro at UCSC, the method works as follows.
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
const vector<char> BASES = {'A', 'C', 'G', 'T'};
const map<char,int> BASE_INDEXES = {{'A',0}, {'C',1}, {'G',2}, {'T',3}};

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
    pair<int, vector<double> > predict_runlength(const Coverage &coverage, AlignedBase base) const;

    AlignedBase predict_consensus_base(const Coverage& coverage) const;

    // This is the primary function of this class. Given a coverage object and consensus base, predict the true
    // run length of the aligned bases at a position
    virtual Consensus operator()(const Coverage&) const;

private:

    /// ---- Attributes ---- ///

    // Defined at initialization, this is the size of the probability matrix generated from the data
    int max_runlength;

    // Boolean flags for configuration
    bool ignore_non_consensus_base_repeats;
    bool predict_gap_runlengths;
    bool count_gaps_as_zeros;

    // p(X|Y) normalized for each Y, where X = observed and Y = True run length
    array<vector<vector<double> >, 4> probability_matrices;


    /// ----- Methods ----- ///

    // Read a single csv containing run length probability matrix
    vector<vector<double> > load_probability_matrix(string file_path);

    // Read each probability matrix from its file and store them in a vector (assuming decibel units, aka base 10)
    // Each delimited table in text should be preceded by a fasta-like header e.g.: ">A" for the base it corresponds to.
    // This converts each line to a vector of doubles, appending probability_matrices according to the matrix header.
    void load_probability_matrices(ifstream& matrix_file);

    // For a given vector of likelihoods over each Y value, normalize by the maximum
    vector<double> normalize_likelihoods(vector<double> x, double x_max) const;

        // Count the number of times each unique repeat was observed, to reduce redundancy in calculating log likelihoods
    map<int,int> factor_repeats(const Coverage& coverage) const;
    map<int,int> factor_repeats(const Coverage& coverage, AlignedBase consensus_base) const;

    // For debugging or exporting
    void print_probability_matrices(char separator=',');
    void print_log_likelihood_vector(vector<double> log_likelihoods);

};

#endif
