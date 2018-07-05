#ifndef CZI_SHASTA_CHRONO_HPP
#define CZI_SHASTA_CHRONO_HPP

#include <chrono>

/*******************************************************************************

Usage pattern (from code in Nanopore2 namespace):

const auto t0 = steady_clock::now();
const auto t1 = steady_clock::now();
const double t01 = seconds(t1-t0);   // Can use auto instead of double.

*******************************************************************************/

namespace ChanZuckerberg {
    namespace Nanopore2 {
        using std::chrono::steady_clock;

        template<class Duration> double seconds(Duration duration)
        {
            return std::chrono::duration<double>(duration).count();
        }
    }
}

#endif
