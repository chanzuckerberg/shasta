#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
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
using namespace shasta;

using Separator = boost::char_separator<char>;
using Tokenizer = boost::tokenizer<Separator>;





// Helper function
void SimpleBayesianConsensusCaller::splitAsDouble(string s, string& separators, vector<double>& tokens){
    Separator separator(separators.c_str());
    Tokenizer tok{s, separator};

    for (string token: tok) {
        tokens.push_back(stod(token));
    }
}


// Helper function
void SimpleBayesianConsensusCaller::splitAsString(string s, string& separators, vector<string>& tokens){
    Separator separator(separators.c_str());
    Tokenizer tok{s, separator};

    for (string token: tok) {
        tokens.push_back(token);
    }
}


void validate_text_file(string configPath){
    ///
    /// Check if all characters in a string are printable, used for checking if input files are valid. If any
    /// non-printables exist, throw an error
    ///
    ifstream file(configPath);
    char character;
    uint64_t c = 0;

    while(file.get(character)) {
        if (not std::isgraph(character) and not std::isspace(character)) {
            throw runtime_error("ERROR: unprintable character detected in bayesian config file"
                                "\n in config file: " + configPath +
                                "\n at position " + to_string(c));
        }
        c++;
    }
}


void SimpleBayesianConsensusCaller::validateMatrixDimensions(string configPath){
    size_t ySize = 0;
    size_t xSize = 0;

    for (size_t baseIndex=0; baseIndex < probabilityMatrices.size(); baseIndex++){
        if (baseIndex > 0){
            if (ySize != probabilityMatrices[baseIndex].size()){
                string base = string(1, AlignedBase::fromInteger(uint8_t(baseIndex)).character());
                throw runtime_error("ERROR: matrix size conflict detected. Matrix " + base +
                " does not match previous base.");
            }
        }

        ySize = probabilityMatrices[baseIndex].size();

        for (size_t yIndex=0; yIndex < ySize; yIndex++) {
            if (yIndex > 0) {
                if (xSize != probabilityMatrices[baseIndex][yIndex].size()) {
                    string base = string(1, AlignedBase::fromInteger(uint8_t(baseIndex)).character());
                    throw runtime_error("ERROR: matrix row size conflict in matrix " + base +
                    " at row " + to_string(yIndex));
                }
            }

            xSize = probabilityMatrices[baseIndex][yIndex].size();
        }
    }

    if (priors[0].empty() || priors[1].empty()){
        throw runtime_error("ERROR: no priors in bayesian config file, or possible error parsing priors"
        "\n in config file: " + configPath);
    }

    for (auto& baseMatrix: probabilityMatrices){
        if (baseMatrix.empty()){
            throw runtime_error("ERROR: no likelihoods in bayesian config file, or possible error parsing likelihoods"
            "\n in config file: " + configPath);
        }
    }

    if (priors[0].size() != priors[1].size()){
        throw runtime_error("ERROR: prior probability vector sizes do not match."
        "\n in config file: " + configPath);
    }

    if (probabilityMatrices[0].size() != priors[0].size()){
        throw runtime_error("ERROR: prior probability vector size (" + to_string(priors[0].size()) +
        ") does not match y (true) dimension size (" + to_string(probabilityMatrices[0].size()) +
        ") of likelihood matrix." +
        "\n in config file: " + configPath);
    }
}


// The constructor string can be either:
// - A name identifying one of the built-in configurations.
// - A path to a configuration file.

SimpleBayesianConsensusCaller::SimpleBayesianConsensusCaller(
    const string& constructorString){
    ignoreNonConsensusBaseRepeats = true;
    predictGapRunlengths = false;
    countGapsAsZeros = false;

    // Try to interpret the constructor string
    // as a name of a built-in configuration.
    const bool isBuiltin = constructBuiltin(constructorString);

    // If it was not a built-in name,
    // interpret the constructor string as a path to
    // a configuration file.
    if(isBuiltin) {
        cout << "Using predefined Bayesian consensus caller " << constructorString << endl;
    } else {

        if(constructorString.size()==0 || constructorString[0]!='/') {
            const string errorMessage = constructorString + " is not the name of a built-in Bayesian model "
                "and therefore it must be an absolute path to a configuration file. "
                "A relative path is not accepted. "
                "Valid built-in choices are: guppy-2.3.1-a, guppy-2.3.5-a, guppy-3.0.5-a, "
                "guppy-3.4.4-a, guppy-3.6.0-a, r10-guppy-3.4.8-a.";
            throw runtime_error(errorMessage);
        }
        ifstream matrixFile(constructorString);
        if (not matrixFile.good()) {
            const string errorMessage = constructorString + " is not a built-in Bayesian model "
                "and could not be open as a configuration file. "
                "Valid built-in choices are: guppy-2.3.1-a, guppy-2.3.5-a, guppy-3.0.5-a";
            throw runtime_error(errorMessage);
        }
        validate_text_file(constructorString);
        loadConfiguration(matrixFile);
        cout << "Loaded Bayesian consensus caller from configuration file " <<
            constructorString << endl;
    }

    validateMatrixDimensions(constructorString);

    maxInputRunlength = uint16_t(probabilityMatrices[0][0].size() - 1);
    maxOutputRunlength = uint16_t(probabilityMatrices[0].size() - 1);

    cout << "Bayesian consensus caller configuration name is " <<
        configurationName << endl;
}



void SimpleBayesianConsensusCaller::printProbabilityMatrices(char separator){
    const uint32_t length = uint(probabilityMatrices[0].size());
    uint32_t nBases = 4;

    for (uint32_t b=0; b<nBases; b++){
        cout << '>' << Base::fromInteger(b).character() << " " << probabilityMatrices[b].size() << '\n';

        for (uint32_t i=0; i<length; i++){
            for (uint32_t j=0; j<length; j++){
                // Print with exactly 9 decimal values
                printf("%.9f",probabilityMatrices[b][i][j]);
                if (j != length-1){
                    cout << separator;
                }
            }
            cout << '\n';
        }
        if (b != nBases-1){
            cout << '\n';
        }
    }
}


void SimpleBayesianConsensusCaller::printPriors(char separator){
    const uint32_t length = uint(priors[0].size());
    uint32_t nBases = 2;

    for (uint32_t b=0; b<nBases; b++){
        cout << '>' << Base::fromInteger(b).character() << " " << priors[b].size() << '\n';

        for (uint32_t i=0; i<length; i++){
            printf("%d %.9f",int(i), priors[b][i]);
            if (i != length-1){
                cout << separator;
            }
        }
        if (b != nBases-1){
            cout << '\n';
        }
    }
}


void SimpleBayesianConsensusCaller::parseName(ifstream& matrixFile, string& line){
    // Expect only one line to follow
    getline(matrixFile, line);
    configurationName = line;
}


void SimpleBayesianConsensusCaller::parsePrior(ifstream& matrixFile, string& line, vector<string>& tokens){
    // Expect only one line to follow
    getline(matrixFile, line);

    // Initialize empty vector to fill with tokens from csv
    vector<double> row;

    // Assume csv format
    string separators = ",";
    splitAsDouble(line, separators, row);

    // Two prior distributions exist. One for AT and one for GC, since observed reads are bidirectional
    if (tokens[0] == "AT"){
        priors[0] = row;
    }
    else if (tokens[0] == "GC"){
        priors[1] = row;
    }
}


void SimpleBayesianConsensusCaller::parseLikelihood(ifstream& matrixFile, string& line, vector<string>& tokens){
    char base;
    uint32_t baseIndex = 0;
    string separators = ",";

    // Expect many lines (usually 51)
    while (getline(matrixFile, line)){

        // Stop iterating lines when blank line is encountered
        if (line.empty()){
            break;
        }

        // Initialize empty vector to fill with tokens from csv
        vector<double> row;

        // Assume csv format
        splitAsDouble(line, separators, row);

        base = tokens[0][0];
        baseIndex = uint32_t(Base::fromCharacter(base).value);

        probabilityMatrices[baseIndex].push_back(row);
    }
}


void SimpleBayesianConsensusCaller::loadConfiguration(ifstream& matrixFile){
    string line;
    string separators = " ";

    while (getline(matrixFile, line)){
        // Ignore empty lines
        if (line.empty()){
            continue;
        }

        // Header line (labeled via fasta-like headers)
        if (line[0] == '>'){
            vector<string> tokens;

            // Store the header
            line = line.substr(1, line.size()-1);
            splitAsString(line, separators, tokens);

            if (tokens[0] == "Name"){
                parseName(matrixFile, line);

            }else if (tokens.size() > 1 and not tokens[0].empty() and tokens[1] == "prior"){
                parsePrior(matrixFile, line, tokens);

            }else if (tokens.size() > 1 and not tokens[0].empty() and tokens[1] == "likelihood"){
                parseLikelihood(matrixFile, line, tokens);
            }
        }
    }
}


void SimpleBayesianConsensusCaller::printLogLikelihoodVector(vector<double>& logLikelihoods){
    int i = 0;
    for (auto& item: logLikelihoods){
        cout << i << " " << pow(10, item) << '\n';
        i++;
    }
}


void SimpleBayesianConsensusCaller::normalizeLikelihoods(vector<double>& x, double xMax) const{
    for (uint32_t i=0; i<x.size(); i++){
        x[i] = x[i]-xMax;
    }
}


void SimpleBayesianConsensusCaller::factorRepeats(
    array<std::map<uint16_t,uint16_t>,2>& factoredRepeats,
    const Coverage& coverage) const{

    // Store counts for each unique observation
    for (auto& observation: coverage.getReadCoverageData() ){
        // If NOT a gap, always increment
        if (not observation.base.isGap()) {
            factoredRepeats[uint16_t(observation.strand)][uint16_t(observation.repeatCount)]++;
        // If IS a gap only increment if "countGapsAsZeros" is true
        }else if (countGapsAsZeros){
            factoredRepeats[uint16_t(observation.strand)][0]++;
        }
    }
}


void SimpleBayesianConsensusCaller::factorRepeats(
    array<std::map<uint16_t,uint16_t>,2>& factoredRepeats,
    const Coverage& coverage,
    AlignedBase consensusBase) const{

    // Store counts for each unique observation
    for (auto& observation: coverage.getReadCoverageData() ){
        // Ignore non consensus repeat values
        if (observation.base.value == consensusBase.value){
            // If NOT a gap, always increment
            if (not observation.base.isGap()) {
                factoredRepeats[uint16_t(observation.strand)][uint16_t(observation.repeatCount)]++;
            // If IS a gap only increment if "countGapsAsZeros" is true
            }else if (countGapsAsZeros){
                factoredRepeats[uint16_t(observation.strand)][0]++;
            }
        }
    }
}


uint16_t SimpleBayesianConsensusCaller::predictRunlength(const Coverage &coverage, AlignedBase consensusBase, vector<double>& logLikelihoodY) const{
    array <std::map <uint16_t,uint16_t>, 2> factoredRepeats;    // Repeats grouped by strand and length

    size_t priorIndex = -1;   // Used to determine which prior probability vector to access (AT=0 or GC=1)
    uint16_t x;               // Element of X = {x_0, x_1, ..., x_i} observed repeats
    uint16_t c;               // Number of times x_i was observed
    uint16_t y;               // Element of Y = {y_0, y_1, ..., y_j} true repeat between 0 and j=max_runlength
    double logSum;            // Product (in logspace) of P(x_i|y_j) for each i

    double yMaxLikelihood = -INF;     // Probability of most probable true repeat length
    uint16_t yMax = 0;                 // Most probable repeat length

    // Determine which index to use for this->priors
    if (consensusBase.character() == 'A' || consensusBase.character() == 'T'){
        priorIndex = 0;
    }
    else if (consensusBase.character() == 'G' || consensusBase.character() == 'C'){
        priorIndex = 1;
    }

    // Count the number of times each unique repeat was observed, to reduce redundancy in calculating log likelihoods/
    // Depending on class boolean "ignoreNonConsensusBaseRepeats" filter out observations
    if (ignoreNonConsensusBaseRepeats) {
        factorRepeats(factoredRepeats, coverage, consensusBase);
    }
    else {
        factorRepeats(factoredRepeats, coverage);
    }

    // Iterate all possible Y from 0 to j to calculate p(Y_j|X) where X is all observations 0 to i,
    // assuming i and j are less than maxRunlength
    for (y = 0; y <= maxOutputRunlength; y++){
        // Initialize logSum for this Y value using empirically determined priors
        logSum = priors[priorIndex][y];

        for (uint16_t strand = 0; strand <= factoredRepeats.size() - 1; strand++){
            for (auto& item: factoredRepeats[strand]){
                x = item.first;
                c = item.second;

                // In the case that observed runlength is too large for the matrix, cap it at maxRunlength
                if (x > maxInputRunlength){
                    x = maxInputRunlength;
                }

                // Increment log likelihood for this y_j
                logSum += double(c)*probabilityMatrices[consensusBase.value][y][x];
            }
        }

        logLikelihoodY[y] = logSum;

        // Update max Y value if log likelihood is greater than previous maximum... Should do this outside this loop?
        if (logSum > yMaxLikelihood){
            yMaxLikelihood = logSum;
            yMax = y;
        }
    }

    normalizeLikelihoods(logLikelihoodY, yMaxLikelihood);

    return max(uint16_t(1), yMax);   // Don't allow zeroes...
}


AlignedBase SimpleBayesianConsensusCaller::predictConsensusBase(const Coverage& coverage) const{
    const vector<CoverageData>& coverageDataVector = coverage.getReadCoverageData();
    vector<uint32_t> baseCounts(5,0);
    uint32_t maxBaseCount = 0;
    uint8_t maxBase = 4;   // Default to gap in case coverage is empty (is this possible?)
    uint32_t key;

    // Count bases. If it's a gap increment placeholder 4 in baseCount vector
    for(const CoverageData& observation: coverageDataVector) {
        if (not observation.base.isGap()) {
            key = observation.base.value;
            baseCounts[key]++;
        }
        else{
            key = 4;
            baseCounts[key]++;
        }
    }

    // Determine most represented base (consensus)
    for (uint32_t i=0; i<5; i++){
        if (baseCounts[i] > maxBaseCount){
            maxBaseCount = baseCounts[i];
            maxBase = uint8_t(i);
        }
    }

    return AlignedBase::fromInteger(maxBase);
}


Consensus SimpleBayesianConsensusCaller::operator()(const Coverage& coverage) const{
    AlignedBase consensusBase;
    uint16_t consensusRepeat;

    vector<double> logLikelihoods(u_long(maxOutputRunlength+1), -INF);    // initialize as zeros in log space

    consensusBase = predictConsensusBase(coverage);

    if (predictGapRunlengths) {
        // Predict all run lengths regardless of whether consensus base is a gap
        consensusRepeat = predictRunlength(coverage, consensusBase, logLikelihoods);
    }
    else{
        if (not consensusBase.isGap()) {
            // Consensus is NOT a gap character, and the configuration forbids predicting gaps
            consensusRepeat = predictRunlength(coverage, consensusBase, logLikelihoods);
        }
        else{
            // Consensus IS a gap character, and the configuration forbids predicting gaps
            consensusRepeat = 0;
        }
    }

    return Consensus(AlignedBase::fromInteger(consensusBase.value), consensusRepeat);
}


// The test program expects as input two csv files:
// - SimpleBayesianConsensusCaller.csv:
//       The configuration file to create the caller.
// - Observations.csv:
//       The observed bases, strands, repeat counts, in this format:
//       A,0,3
//       A,1,4
void shasta::testSimpleBayesianConsensusCaller(
    const string& configurationFileName)
{
    // Read the observations.
    Coverage coverage;
    ifstream input("Observations.csv");
    while(true) {

        // Get a line.
        string line;
        getline(input, line);
        if(!input) {
            break;
        }
        cout << line << endl;

        // Parse it.
        vector<string> tokens;
        boost::algorithm::split(tokens, line, boost::algorithm::is_any_of(","));
        SHASTA_ASSERT(tokens.size() == 3);

        // Add it to our Coverage object.
        const char baseCharacter = tokens[0][0];
        coverage.addRead(AlignedBase::fromCharacter(baseCharacter),
            boost::lexical_cast<int>(tokens[1]),
            boost::lexical_cast<int>(tokens[2]));
    }



    // Construct the caller.
    // This expects the configuration file to be named
    // SimpleBayesianConsensusCaller.csv.
    SimpleBayesianConsensusCaller caller(configurationFileName);

    // Compute consensus.
    const Consensus consensus = caller(coverage);

    // Write it out.
    cout << "Consensus: " << consensus.base << " " << consensus.repeatCount << '\n';
}
