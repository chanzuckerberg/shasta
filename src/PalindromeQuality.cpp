#include "PalindromeQuality.hpp"
#include "span.hpp"
#include <cmath>

using std::runtime_error;
using std::cout;
using std::cerr;
using std::pow;

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

using namespace boost::accumulators;
using namespace shasta;

typedef accumulator_set<float, features<tag::count, tag::mean, tag::variance>> stats_accumulator;


double qualityCharToErrorProbability(char q) {
    return pow(10, double(q - 33) / -10.0);
}


bool shasta::isPalindromic(span<char> qualities){
    double relativeMeanDifference = 0.09;
    double minimumMean = 0.15;
    double minimumVariance = 0.025;

    bool isPalindromic = false;

    stats_accumulator leftStats;
    stats_accumulator rightStats;

    auto length = qualities.size();

    // Don't bother classifying any reads that would have less than 3 scores per side
    if (length < 6){
        return false;
    }

    size_t midpoint = length/2;

    for (size_t i=0; i<midpoint; i++){
        auto q = qualities[i];
        auto p = qualityCharToErrorProbability(q);

        leftStats(p);
    }

    for (size_t i=midpoint; i<length; i++){
        auto q = qualities[i];
        auto p = qualityCharToErrorProbability(q);

        rightStats(p);
    }

    float leftMean = mean(leftStats);
    float leftVariance = variance(leftStats);

    float rightMean = mean(rightStats);
    float rightVariance = variance(rightStats);

    // Compare the mean and variance using thresholds derived empirically from some palindromic reads
    if (rightMean - leftMean > relativeMeanDifference and rightMean >= minimumMean){
        if (rightVariance > leftVariance and rightVariance > minimumVariance){
            isPalindromic = true;
        }
    }

    return isPalindromic;
}
