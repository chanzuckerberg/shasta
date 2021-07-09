#ifndef SHASTA_HSV_HPP
#define SHASTA_HSV_HPP

#include "algorithm.hpp"
#include "utility.hpp"

namespace shasta {
    pair<double, double> hsvToHsl(double SV, double V);
}


// HSV to HSL conversion.
// Given S and V for an HSV color, returns the S and L
// values for the corresponding HSL color.
// See https://en.wikipedia.org/wiki/HSL_and_HSV#HSV_to_HSL
inline std::pair<double, double> shasta::hsvToHsl(double SV, double V)
{
    const double L = V * (1. - 0.5 * SV);
    const double denominator = min(L, 1. - L);
    const double SL = (denominator == 0.) ? 0. : (V - L)/denominator;
    return make_pair(SL, L);
}

#endif
