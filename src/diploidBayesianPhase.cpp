#include "diploidBayesianPhase.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;

#include "iostream.hpp"
#include <cmath>
#include "tuple.hpp"

// See comments in diploidBayesianPhase.hpp for the meaning
// of the arguments and the return value.


static void writeMatrix(const array<array<double, 2>, 2>& m)
{
    for(uint64_t side0=0; side0<2; side0++) {
        for(uint64_t side1=0; side1<2; side1++) {
            cout << m[side0][side1] << " ";
        }
        cout << endl;
    }

}


pair<double, double> shasta::diploidBayesianPhase(
    const array<array<uint64_t, 2>, 2>& matrix,
    double epsilon)
{
    const bool debug = false;

    // Construct a version of the matrix containing doubles,
    // so we don't need to keep converting.
    array<array<double, 2>, 2> m;
    for(uint64_t side0=0; side0<2; side0++) {
        for(uint64_t side1=0; side1<2; side1++) {
            m[side0][side1] = double(matrix[side0][side1]);
        }
    }

    // Contract the matrix over side1 to compute n0.
    // n0[side0] is the number of common reads seen on
    // side side0 of bubble0 regardless of where they are on bubble1.
    array<double, 2> n0;
    for(uint64_t side0=0; side0<2; side0++) {
        n0[side0] = 0.;
        for(uint64_t side1=0; side1<2; side1++) {
            n0[side0] += m[side0][side1];
        }
    }

    // Contract the matrix over side0 to compute n1.
    // n1[side1] is the number of common reads seen on
    // side side1 of bubble1 regardless of where they are on bubble0.
    array<double, 2> n1;
    for(uint64_t side1=0; side1<2; side1++) {
        n1[side1] = 0.;
        for(uint64_t side0=0; side0<2; side0++) {
            n1[side1] += m[side0][side1];
        }
    }

    // Compute the total number of common reads, n.
    const double n = n0[0] + n0[1];
    SHASTA_ASSERT(n == n1[0] + n1[1]);
    const double nm2 = 1. / (n * n);

    // Compute Prandom[side0][side1], the probability that a
    // read is on side side0 of bubble0 and on side side1 of bubble1
    // under the random hypothesis.
    array<array<double, 2>, 2> Prandom;
    for(uint64_t side0=0; side0<2; side0++) {
        for(uint64_t side1=0; side1<2; side1++) {
            Prandom[side0][side1] = nm2 * n0[side0] * n1[side1];
        }
    }

    // Compute Pin[side0][side1], the probability that a
    // read is on side side0 of bubble0 and on side side1 of bubble1
    // under the in phase hypothesis.
    array<array<double, 2>, 2> Pin;
    for(uint64_t side0=0; side0<2; side0++) {
        for(uint64_t side1=0; side1<2; side1++) {
            Pin[side0][side1] = epsilon * Prandom[side0][side1];
        }
    }
    const double factorIn = (1. - epsilon) / (n0[0] * n1[0] + n0[1] * n1[1]);
    Pin[0][0] += factorIn * (n0[0] * n1[0]);
    Pin[1][1] += factorIn * (n0[1] * n1[1]);

    // Compute Pout[side0][side1], the probability that a
    // read is on side side0 of bubble0 and on side side1 of bubble1
    // under the out of phase hypothesis.
    array<array<double, 2>, 2> Pout;
    for(uint64_t side0=0; side0<2; side0++) {
        for(uint64_t side1=0; side1<2; side1++) {
            Pout[side0][side1] = epsilon * Prandom[side0][side1];
        }
    }
    const double factorOut = (1. - epsilon) / (n0[0] * n1[1] + n0[1] * n1[0]);
    Pout[0][1] += factorOut * (n0[0] * n1[1]);
    Pout[1][0] += factorOut * (n0[1] * n1[0]);

    if(debug) {
        cout << "Prandom" << endl;
        writeMatrix(Prandom);
        cout << "Pin" << endl;
        writeMatrix(Pin);
        cout << "Pout" << endl;
        writeMatrix(Pout);
    }



    // Now compute probability ratios, conditional to the observed
    // distribution of reads.
    double logPin = 0.;
    double logPout = 0.;
    for(uint64_t side0=0; side0<2; side0++) {
        for(uint64_t side1=0; side1<2; side1++) {
            logPin += m[side0][side1] * 10. * std::log10(Pin[side0][side1] / Prandom[side0][side1]);
            logPout += m[side0][side1] * 10. * std::log10(Pout[side0][side1] / Prandom[side0][side1]);
        }
    }

    return make_pair(logPin, logPout);
}



void shasta::testDiploidBayesianPhase(
    double epsilon,
    uint64_t m00,
    uint64_t m01,
    uint64_t m10,
    uint64_t m11
    )
{
    const array<array<uint64_t, 2>, 2> matrix = {m00, m01, m10, m11};

    double logPin, logPout;
    tie(logPin, logPout) = diploidBayesianPhase(matrix, epsilon);
    cout << logPin << " " << logPout << endl;
}

