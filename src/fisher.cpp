#include "fisher.hpp"
using namespace shasta;

#include <boost/math/distributions/hypergeometric.hpp>
#include "algorithm.hpp"



// Fisher test for a 2 by 2 matrix of frequencies.
// See https://en.wikipedia.org/wiki/Fisher%27s_exact_test.
// Returns log(P) in decibels (dB). High is good.
double shasta::logFisher(
    uint32_t a,
    uint32_t b,
    uint32_t c,
    uint32_t d)
{

    // Construct the probability distribution for c.
    // Names are as in the definition of boost::math::hypergeometric_distribution.
    const uint32_t N = a + b + c + d;
    const uint32_t r = a + c;
    const uint32_t n = c + d;
    boost::math::hypergeometric_distribution<> h(r, n, N);

    // Probability of the observed values.
    double P0 = boost::math::pdf(h, c);

    // Sum the probabilities of all possibilities
    // with probability not greater than this.
    double sum = 0.;
    uint32_t kFirst = uint32_t(max(0, int32_t(r + n - N)));
    uint32_t kLast = uint32_t(min(r, uint32_t(n)));
    for(uint32_t k=kFirst; k<=kLast; k++) {
        const double P = boost::math::pdf(h, k);
        if(P <= P0) {
            sum += P;
        }
    }
    return -10. * log10(sum);
}
