#ifndef SHASTA_FISHER_HPP
#define SHASTA_FISHER_HPP

#include "cstdint.hpp"

// Fisher test for a 2 by 2 matrix of frequencies.
// See https://en.wikipedia.org/wiki/Fisher%27s_exact_test.
// Returns log(P) in decibels (dB). High is good.

namespace shasta {

    double logFisher(
        uint32_t a,
        uint32_t b,
        uint32_t c,
        uint32_t d);
}

#endif

