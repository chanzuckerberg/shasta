#ifndef SHASTA_PALINDROMEQUALITY_HPP
#define SHASTA_PALINDROMEQUALITY_HPP

#include "span.hpp"

namespace shasta {

bool isPalindromic(span<char> qualities,
                   double relativeMeanDifference,
                   double minimumMean,
                   double minimumVariance);

}

#endif //SHASTA_PALINDROMEQUALITY_HPP
