#include "PalindromeQuality.hpp"

#include <iostream>
#include <stdexcept>
#include "span.hpp"

using std::runtime_error;
using std::cout;
using std::cerr;

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>


bool shasta::classify_palindromic_q_scores(span<char> qualities){
    bool is_palindromic = false;

    stats_accumulator left_stats;
    stats_accumulator right_stats;

    auto length = qualities.size();

    // Don't bother classifying any reads that would have less than 3 scores per side
    if (length < 6){
        return false;
    }

    size_t midpoint = length/2;

    for (size_t i=0; i<midpoint; i++){
        auto q = qualities[i];
        auto p = quality_char_to_error_probability(q);

        left_stats(p);
    }

    for (size_t i=midpoint; i<length; i++){
        auto q = qualities[i];
        auto p = quality_char_to_error_probability(q);

        right_stats(p);
    }

    float left_mean = mean(left_stats);
    float left_variance = variance(left_stats);

    float right_mean = mean(right_stats);
    float right_variance = variance(right_stats);

    // Compare the mean and variance using thresholds derived empirically from some palindormic reads
    if (right_mean - left_mean > 0.09 and right_mean >= 0.15){
        if (right_variance > left_variance and right_variance > 0.025){
            is_palindromic = true;
        }
    }

    return is_palindromic;
}
