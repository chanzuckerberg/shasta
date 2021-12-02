#ifndef SHASTA_DIPLOID_BAYESIAN_PHASE_HPP
#define SHASTA_DIPLOID_BAYESIAN_PHASE_HPP

#include "array.hpp"
#include "cstdint.hpp"
#include "utility.hpp"

/*******************************************************************************

Function diploidBayesianPhase uses a Bayesian model to evaluate the
phasing of two bubbles relative to each other.

Call the two bubbles bubble0 and bubble1.
Each bubble has two sides (branches), 0 and 1.

On input:

matrix[side0][side1] is the number of common reads between the two bubbles
that fall on side side0 of bubble0 and on side side1 of bubble1.

epsilon is the fraction of reads in error assumed by the Bayesian model.

The Bayesian model consider three hypotheses:

- The random hypothesis in which the two bubbles are uncorrelated.
  This usually means that one of the two bubbles is the result of errors.

- The in-phase hypothesis in which the two bubbles are in phase
  relative to each other. Under this hypothesis, all common reads that
  are not errors are either:
  * on side 0 of bubble0 and on side 0 of bubble1, or
  * on side 1 of bubble0 and on side 1 of bubble1.

- The out-of-phase hypothesis in which the two bubbles are out of phase
  relative to each other. Under this hypothesis, all common reads that
  are not errors are either:
  * on side 0 of bubble0 and on side 1 of bubble1, or
  * on side 1 of bubble0 and on side 0 of bubble1.

The Bayesian model computes the probability ratios
Pin/Prandom and Pout/Prandom conditional to the observed
read distribution (and assuming neutral priors for the
three hypotheses).

On exit, diploidBayesianPhase returns a pair containing
log(Pin/Prandom) and log(Pout/Prandom) expressed in decibels (dB).

*******************************************************************************/

namespace shasta {


    pair<double, double> diploidBayesianPhase(
        const array<array<uint64_t, 2>, 2>& matrix,
        double epsilon
    );

    void testDiploidBayesianPhase(
        double epsilon,
        uint64_t m00,
        uint64_t m01,
        uint64_t m10,
        uint64_t m11);

}


#endif
