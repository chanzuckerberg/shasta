#include <boost/tokenizer.hpp>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include "SimpleBayesianConsensusCaller.hpp"
#include "Coverage.hpp"
#include "ConsensusCaller.hpp"

using ChanZuckerberg::shasta::Consensus;
using Separator = boost::char_separator<char>;
using Tokenizer = boost::tokenizer<Separator>;
using std::ifstream;
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::pow;
using std::map;

using namespace ChanZuckerberg;
using namespace shasta;


// The constructor does not have any parameters.
// All data should be read from a file with fixed name
// in the run directory. We will update the documentation accordingly.
SimpleBayesianConsensusCaller::SimpleBayesianConsensusCaller(){
    this->load_probability_matrices(file_paths);
    this->max_runlength = 50;
    this->ignore_non_consensus_base_repeats = false;
}


void SimpleBayesianConsensusCaller::print_probability_matrices(char separator){
    int length = this->probability_matrices[0].size();
    int n_bases = 4;

    for (int b=0; b<n_bases; b++){
        cout << '>' << bases[b] << '\n';
        for (int i=0; i<length; i++){
            for (int j=0; j<length; j++){
                cout << this->probability_matrices[b][i][j];
                if (j != length-1)
                    cout << separator;
            }
            cout << '\n';
        }
        if (b != n_bases-1)
            cout << '\n';
    }
}


vector<double> split(string s, char separator_char){
    vector<double> tokens;
    Separator separator(&separator_char);
    Tokenizer tok{s, separator};

    for (string token: tok) {
        tokens.push_back(stod(token));
    }

    return tokens;
}


void SimpleBayesianConsensusCaller::load_probability_matrices(vector<string> matrix_file_paths){
    // Assume that matrix file paths are in order of base convention A=0, C=1, G=2, T=3
    for (string& file_path: matrix_file_paths){
        this->probability_matrices.push_back(this->load_probability_matrix(file_path));
    }
}


vector<vector<double> > SimpleBayesianConsensusCaller::load_probability_matrix(string file_path){
    ifstream matrix_file(file_path);
    vector<vector<double> > matrix;
    vector<double> row;
    string line;

    // Ensure that file pointer is not null
    if (not matrix_file.good()){
        cout << "ERROR: file read error: " << file_path << '\n';
        exit(1);
    }

    // Read all lines in file and convert each to a vector of doubles, appending each to matrix
    while (getline(matrix_file, line)){
        // Assume csv format
        row = split(line, ',');
        matrix.push_back(row);
    }

    return matrix;
}


void SimpleBayesianConsensusCaller::print_log_likelihood_vector(vector<double> log_likelihoods){
    int i = 0;
    for (auto& item: log_likelihoods){
        cout << i << " " << pow(10, item) << '\n';
        i++;
    }
}


vector<double> SimpleBayesianConsensusCaller::normalize_likelihoods(vector<double> x, double x_max) const{
    vector<double> x_normalized(x.size());

    for (int i=0; i<x.size(); i++){
        x_normalized[i] = x[i]-x_max;
    }

    return x_normalized;
}


map<int,int> SimpleBayesianConsensusCaller::factor_repeats(const Coverage& coverage) const{
    // Store counts for each unique observation
    map<int,int> factored_repeats;

    for (auto& observation: coverage.getReadCoverageData() ){
        if (not observation.base.isGap()) {
            factored_repeats[observation.repeatCount]++;
        }
        else{
            // TODO: address the issue of whether to count gaps (zeros) in runlength prediction...
            factored_repeats[0]++;
        }
    }

    return factored_repeats;
}


map<int,int> SimpleBayesianConsensusCaller::factor_repeats(const Coverage& coverage, AlignedBase consensus_base) const{
    // Store counts for each unique observation
    map<int,int> factored_repeats;

    for (auto& observation: coverage.getReadCoverageData() ){
        if (observation.base.value == consensus_base.value){
            if (not observation.base.isGap()) {
                factored_repeats[observation.repeatCount]++;
            }
            else {
                // TODO: address the issue of whether to count gaps (zeros) in runlength prediction...
                factored_repeats[0]++;
            }
        }
    }

    return factored_repeats;
}


//TODO: convert "base" to proper Base class?
pair<int, vector<double> > SimpleBayesianConsensusCaller::predict_runlength(const Coverage &coverage, AlignedBase consensusBase) const{
    // Count the number of times each unique repeat was observed, to reduce redundancy in calculating log likelihoods/
    // Depending on class boolean "ignore_non_consensus_base_repeats" filter out observations
    map<int, int> factored_repeats;

    if (this->ignore_non_consensus_base_repeats)
        factored_repeats = this->factor_repeats(coverage, consensusBase);
    else
        factored_repeats = this->factor_repeats(coverage);

    int x_i;    // Member of X such that X = {x_0, x_1, ..., x_i} observed repeats
    int c_i;    // Number of times x_i was observed
    int y_j;    // Member of Y such that Y = {y_0, y_1, ..., y_j} true repeat between 0 and j=max_runlength (50)

    double log_sum;                                                 // Product (in logspace) of P(x_i|y_j) for each i
    vector<double> log_likelihood_y(this->max_runlength, -inf);     // Loglikelihoods for all possible y_j

    double y_max_likelihood = -inf;
    int y_max = -1;

    // Iterate all possible Y from 0 to j to calculate p(Y_j|X) where X is all observations 0 to i,
    // assuming i and j are less than max_runlength
    for (y_j = 0; y_j <= this->max_runlength; y_j++){
        // Initialize log_sum for this Y value
        log_sum = 0;

        for (auto& item: factored_repeats){
            x_i = item.first;
            c_i = item.second;

            // In the case that observed runlength is too large for the matrix, cap it at max_runlength
            if (x_i > this->max_runlength){
                x_i = this->max_runlength;
            }

            // Increment log likelihood for this y_j
            log_sum += c_i*this->probability_matrices[consensusBase.value][y_j][x_i];
        }

        log_likelihood_y[y_j] = log_sum;

        // Update max Y value if log likelihood is greater than previous maximum... Should do this outside this loop?
        if (log_sum > y_max_likelihood){
            y_max_likelihood = log_sum;
            y_max = y_j;
        }
    }

    log_likelihood_y = this->normalize_likelihoods(log_likelihood_y, y_max_likelihood);

    return make_pair(y_max, log_likelihood_y);
}


AlignedBase SimpleBayesianConsensusCaller::predict_consensus_base(const Coverage& coverage) const{
    const vector<CoverageData>& coverage_data_vector = coverage.getReadCoverageData();
    vector<int> base_keys = {0, 1, 2, 3, 4};
    vector<int> base_counts(base_keys.size(),0);
    int max_base_count = 0;
    int max_base;
    int key;

    // Count bases. If it's a gap increment placeholder 4 in base_count vector
    for(const CoverageData& observation: coverage_data_vector) {
        if (not observation.base.isGap()) {
            key = (int)observation.base.value;
            base_counts[key]++;
        }
        else{
            key = 4;
            base_counts[key]++;
        }
    }

    // Determine most represented base (consensus)
    for (auto& key: base_keys){
        if (base_counts[key] > max_base_count){
            max_base_count = base_counts[key];
            max_base = key;
        }
    }

    return AlignedBase::fromInteger((uint8_t)max_base);
}


Consensus SimpleBayesianConsensusCaller::operator()(const Coverage& coverage) const{
    // TODO: test that coverage is not empty?
    // TODO: test that consensus base is not gap
    AlignedBase consensus_base;
    int consensus_repeat;
    vector<double> log_likelihoods;

    consensus_base = this->predict_consensus_base(coverage);
    tie(consensus_repeat, log_likelihoods) = this->predict_runlength(coverage, consensus_base);

    // This returns a A with repeat count 3.
    return Consensus(AlignedBase::fromCharacter(consensus_base.character()), consensus_repeat);
}
