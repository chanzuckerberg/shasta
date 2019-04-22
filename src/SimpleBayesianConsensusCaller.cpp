#include <boost/tokenizer.hpp>
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <string>
#include <cstdio>
#include <array>
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
using std::array;
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
    ignore_non_consensus_base_repeats = true;
    predict_gap_runlengths = false;
    count_gaps_as_zeros = false;

    const string fileName = "SimpleBayesianConsensusCaller.csv";
    ifstream matrix_file(fileName);
    if (not matrix_file.good()) {
        const string error_message = "Error opening file: " + fileName;
        throw runtime_error(error_message);
    }

    load_probability_matrices(matrix_file);
}


void SimpleBayesianConsensusCaller::print_probability_matrices(char separator){
    const uint32_t length = uint(probability_matrices[0].size());
    uint32_t n_bases = 4;

    for (uint32_t b=0; b<n_bases; b++){
        cout << '>' << Base::fromInteger(b).character() << '\n';

        for (uint32_t i=0; i<length; i++){
            for (uint32_t j=0; j<length; j++){
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
    uint32_t base_index = 0;

    // Matrices are labeled via fasta-like headers
    while (getline(matrix_file, line)){

        // Header line
        if (line[0] == '>'){
            base = line[1];
            base_index = uint32_t(Base::fromCharacter(base).value);
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
    for (uint32_t i=0; i<x.size(); i++){
        x[i] = x[i]-x_max;
    }
}


void SimpleBayesianConsensusCaller::factor_repeats(array<map<uint16_t,uint16_t>,2>& factored_repeats, const Coverage& coverage) const{
    // Store counts for each unique observation
    for (auto& observation: coverage.getReadCoverageData() ){
        // If NOT a gap, always increment
        if (not observation.base.isGap()) {
            factored_repeats[uint16_t(observation.strand)][uint16_t(observation.repeatCount)]++;
        // If IS a gap only increment if "count_gaps_as_zeros" is true
        }else if (count_gaps_as_zeros){
            factored_repeats[uint16_t(observation.strand)][0]++;
        }
    }
}


void SimpleBayesianConsensusCaller::factor_repeats(array<map<uint16_t,uint16_t>,2>& factored_repeats, const Coverage& coverage, AlignedBase consensus_base) const{
    // Store counts for each unique observation
    for (auto& observation: coverage.getReadCoverageData() ){
        // Ignore non consensus repeat values
        if (observation.base.value == consensus_base.value){
            // If NOT a gap, always increment
            if (not observation.base.isGap()) {
                factored_repeats[uint16_t(observation.strand)][uint16_t(observation.repeatCount)]++;
            // If IS a gap only increment if "count_gaps_as_zeros" is true
            }else if (count_gaps_as_zeros){
                factored_repeats[uint16_t(observation.strand)][0]++;
            }
        }
    }
}


uint16_t SimpleBayesianConsensusCaller::predict_runlength(const Coverage &coverage, AlignedBase consensusBase, vector<double>& log_likelihood_y) const{
    // Count the number of times each unique repeat was observed, to reduce redundancy in calculating log likelihoods/
    // Depending on class boolean "ignore_non_consensus_base_repeats" filter out observations
    array <map <uint16_t,uint16_t>, 2> factored_repeats;

    if (ignore_non_consensus_base_repeats) {
        factor_repeats(factored_repeats, coverage, consensusBase);
    }
    else {
        factor_repeats(factored_repeats, coverage);
    }

    uint16_t x_i;    // Element of X = {x_0, x_1, ..., x_i} observed repeats
    uint16_t c_i;    // Number of times x_i was observed
    uint16_t y_j;    // Element of Y = {y_0, y_1, ..., y_j} true repeat between 0 and j=max_runlength (50)

    double log_sum;  // Product (in logspace) of P(x_i|y_j) for each i

    double y_max_likelihood = -INF;
    uint16_t y_max = 0;

    // Iterate all possible Y from 0 to j to calculate p(Y_j|X) where X is all observations 0 to i,
    // assuming i and j are less than max_runlength
    for (y_j = 0; y_j <= max_runlength; y_j++){
        // Initialize log_sum for this Y value.
        // Use a prior (penalty) of a factor 1/4 for each new base.
        log_sum = (y_j==0) ? -20. : (-y_j * log10(4.));

        for (uint16_t strand = 0; strand <= factored_repeats.size() - 1; strand++){
            for (auto& item: factored_repeats[strand]){
                x_i = item.first;
                c_i = item.second;

                // In the case that observed runlength is too large for the matrix, cap it at max_runlength
                if (x_i > max_runlength){
                    x_i = max_runlength;
                }

                // Increment log likelihood for this y_j
                log_sum += double(c_i)*probability_matrices[consensusBase.value][y_j][x_i];
            }
        }

        log_likelihood_y[y_j] = log_sum;

        // Update max Y value if log likelihood is greater than previous maximum... Should do this outside this loop?
        if (log_sum > y_max_likelihood){
            y_max_likelihood = log_sum;
            y_max = y_j;
        }
    }

    normalize_likelihoods(log_likelihood_y, y_max_likelihood);

    return y_max;
}


AlignedBase SimpleBayesianConsensusCaller::predict_consensus_base(const Coverage& coverage) const{
    const vector<CoverageData>& coverage_data_vector = coverage.getReadCoverageData();
    vector<uint32_t> base_counts(5,0);
    uint32_t max_base_count = 0;
    uint8_t max_base = 4;   // Default to gap in case coverage is empty (is this possible?)
    uint32_t key;

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
    for (uint32_t i=0; i<5; i++){
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
    uint16_t consensus_repeat;
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


void testSimpleBayesianConsensusCaller(){
    SimpleBayesianConsensusCaller classifier;
    Coverage coverage;

    coverage.addRead(AlignedBase::fromInteger((uint8_t)1), 1, 1);    // Arguments are base, strand, repeat count.
    coverage.addRead(AlignedBase::fromInteger((uint8_t)1), 0, 2);
    coverage.addRead(AlignedBase::fromInteger((uint8_t)2), 1, 3);
    coverage.addRead(AlignedBase::fromInteger((uint8_t)1), 0, 2);
    coverage.addRead(AlignedBase::fromInteger((uint8_t)1), 1, 2);
    coverage.addRead(AlignedBase::fromInteger((uint8_t)2), 0, 2);
    coverage.addRead(AlignedBase::fromInteger((uint8_t)4), 0, 0);

    AlignedBase consensus_base;

    consensus_base = coverage.mostFrequentBase();

    cout << "CONSENSUS BASE = " << consensus_base << "\n";

    const Consensus consensus = classifier(coverage);

    cout << consensus.base << " " << consensus.repeatCount << '\n';
}
