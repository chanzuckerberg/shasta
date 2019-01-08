#include <boost/tokenizer.hpp>
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <map>
#include "SimpleBayesianConsensusCaller.hpp"
#include "Coverage.hpp"
#include "ConsensusCaller.hpp"

using ChanZuckerberg::shasta::Consensus;
using Separator = boost::char_separator<char>;
using Tokenizer = boost::tokenizer<Separator>;
using std::runtime_error;
using std::ifstream;
using std::vector;
using std::string;
using std::printf;
using std::pair;
using std::make_pair;
using std::cout;
using std::pow;
using std::map;

using namespace ChanZuckerberg;
using namespace shasta;


// Helper function
void SimpleBayesianConsensusCaller::split(string s, char separator_char, vector<double>& tokens){
    Separator separator(&separator_char);
    Tokenizer tok{s, separator};

    for (string token: tok) {
        tokens.push_back(stod(token));
    }
}


SimpleBayesianConsensusCaller::SimpleBayesianConsensusCaller(){
    max_runlength = 50;
    ignore_non_consensus_base_repeats = false;
    predict_gap_runlengths = false;
    count_gaps_as_zeros = false;

    const string fileName = "SimpleBayesianConsensusCaller-1.csv";
    ifstream matrix_file(fileName);
    if (not matrix_file.good()) {
        const string error_message = "Error opening file: " + fileName;
        throw runtime_error(error_message);
    }

    load_probability_matrices(matrix_file);
}


void SimpleBayesianConsensusCaller::print_probability_matrices(char separator){
    const unsigned int length = uint(probability_matrices[0].size());
    unsigned int n_bases = 4;

    for (unsigned long b=0; b<n_bases; b++){
        cout << '>' << Base::fromInteger(b).character() << '\n';

        for (unsigned long i=0; i<length; i++){
            for (unsigned long j=0; j<length; j++){
                // Print with exactly 9 decimal values
                printf("%.9f",probability_matrices[b][i][j]);
                if (j != length-1)
                    cout << separator;
            }
            cout << '\n';
        }
        if (b != n_bases-1)
            cout << '\n';
    }
}


void SimpleBayesianConsensusCaller::load_probability_matrices(ifstream& matrix_file){
    string line;
    char base;
    unsigned long base_index;

    // Matrices are labeled via fasta-like headers
    while (getline(matrix_file, line)){

        // Header line
        if (line[0] == '>'){
            base = line[1];
            base_index = ulong(Base::fromCharacter(base).value);
        }

        // Data line
        else if (line[0] != '>' && not line.empty()) {
            // Initialize empty vector to fill with tokens from csv
            vector<double> row;

            // Assume csv format
            split(line, ',', row);
            probability_matrices[base_index].push_back(row);
        }
    }
}


void SimpleBayesianConsensusCaller::print_log_likelihood_vector(vector<double>& log_likelihoods){
    int i = 0;
    for (auto& item: log_likelihoods){
        cout << i << " " << pow(10, item) << '\n';
        i++;
    }
}


void SimpleBayesianConsensusCaller::normalize_likelihoods(vector<double>& x, double x_max) const{
    for (unsigned long i=0; i<x.size(); i++){
        x[i] = x[i]-x_max;
    }
}


map<int,int> SimpleBayesianConsensusCaller::factor_repeats(const Coverage& coverage) const{
    map<int,int> factored_repeats;

    // Store counts for each unique observation
    for (auto& observation: coverage.getReadCoverageData() ){
        if (not observation.base.isGap()) {
            factored_repeats[int(observation.repeatCount)]++;
        }
        else if (count_gaps_as_zeros){
            factored_repeats[0]++;
        }
    }

    return factored_repeats;
}


map<int,int> SimpleBayesianConsensusCaller::factor_repeats(const Coverage& coverage, AlignedBase consensus_base) const{
    map<int,int> factored_repeats;

    // Store counts for each unique observation
    for (auto& observation: coverage.getReadCoverageData() ){
        // Ignore non consensus repeat values
        if (observation.base.value == consensus_base.value){
            if (not observation.base.isGap()) {
                factored_repeats[int(observation.repeatCount)]++;
            }
            else if (count_gaps_as_zeros){
                factored_repeats[0]++;
            }
        }
    }

    return factored_repeats;
}


int SimpleBayesianConsensusCaller::predict_runlength(const Coverage &coverage, AlignedBase consensusBase, vector<double>& log_likelihood_y) const{
    // Count the number of times each unique repeat was observed, to reduce redundancy in calculating log likelihoods/
    // Depending on class boolean "ignore_non_consensus_base_repeats" filter out observations
    map<int, int> factored_repeats;

    if (ignore_non_consensus_base_repeats) {
        factored_repeats = factor_repeats(coverage, consensusBase);
    }
    else {
        factored_repeats = factor_repeats(coverage);
    }

    unsigned long x_i;    // Element of X = {x_0, x_1, ..., x_i} observed repeats
    unsigned long c_i;    // Number of times x_i was observed
    unsigned long y_j;    // Element of Y = {y_0, y_1, ..., y_j} true repeat between 0 and j=max_runlength (50)

    double log_sum;                                                   // Product (in logspace) of P(x_i|y_j) for each i

    double y_max_likelihood = -INF;
    int y_max = -1;

    // Iterate all possible Y from 0 to j to calculate p(Y_j|X) where X is all observations 0 to i,
    // assuming i and j are less than max_runlength
    for (y_j = 0; y_j <= ulong(max_runlength); y_j++){
        // Initialize log_sum for this Y value
        log_sum = 0;

        for (auto& item: factored_repeats){
            x_i = ulong(item.first);
            c_i = ulong(item.second);

            // In the case that observed runlength is too large for the matrix, cap it at max_runlength
            if (x_i > ulong(max_runlength)){
                x_i = ulong(max_runlength);
            }

            // Increment log likelihood for this y_j
            log_sum += double(c_i)*probability_matrices[consensusBase.value][y_j][x_i];
        }

        log_likelihood_y[y_j] = log_sum;

        // Update max Y value if log likelihood is greater than previous maximum... Should do this outside this loop?
        if (log_sum > y_max_likelihood){
            y_max_likelihood = log_sum;
            y_max = int(y_j);
        }
    }

    normalize_likelihoods(log_likelihood_y, y_max_likelihood);

    return y_max;
}


AlignedBase SimpleBayesianConsensusCaller::predict_consensus_base(const Coverage& coverage) const{
    const vector<CoverageData>& coverage_data_vector = coverage.getReadCoverageData();
    vector<int> base_counts(5,0);
    int max_base_count = 0;
    uint8_t max_base = 4;   // Default to gap in case coverage is empty (is this possible?)
    unsigned long key;

    // Count bases. If it's a gap increment placeholder 4 in base_count vector
    for(const CoverageData& observation: coverage_data_vector) {
        if (not observation.base.isGap()) {
            key = observation.base.value;
            base_counts[key]++;
        }
        else{
            key = 4;
            base_counts[key]++;
        }
    }

    // Determine most represented base (consensus)
    for (unsigned long i=0; i<5; i++){
        if (base_counts[i] > max_base_count){
            max_base_count = base_counts[i];
            max_base = uint8_t(i);
        }
    }

    return AlignedBase::fromInteger(max_base);
}


Consensus SimpleBayesianConsensusCaller::operator()(const Coverage& coverage) const{
    // TODO: test that coverage is not empty?
    AlignedBase consensus_base;
    int consensus_repeat;
    vector<double> log_likelihoods(u_long(max_runlength), -INF);    // initialize as zeros in log space

    consensus_base = predict_consensus_base(coverage);

    if (predict_gap_runlengths) {
        // Predict all run lengths regardless of whether consensus base is a gap
        consensus_repeat = predict_runlength(coverage, consensus_base, log_likelihoods);
    }
    else {
        if (not consensus_base.isGap()) {
            // Consensus is NOT a gap character, and the configuration forbids predicting gaps
            consensus_repeat = predict_runlength(coverage, consensus_base, log_likelihoods);
        } else {
            // Consensus IS a gap character, and the configuration forbids predicting gaps
            consensus_repeat = 0;
        }
    }

    return Consensus(AlignedBase::fromInteger(consensus_base.value), consensus_repeat);
}
