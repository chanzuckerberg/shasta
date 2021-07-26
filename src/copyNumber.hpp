#ifndef SHASTA_COPY_NUMBER_HPP
#define SHASTA_COPY_NUMBER_HPP

#include "prefixLength.hpp"
#include "span.hpp"
#include "cstdint.hpp"

namespace shasta {

    // Figure out if two sequences differ only by copy numbers in
    // a repeat with given period, 2 <= period <= maxPeriod.
    // If this is the case, returns the shortest period for which this is true.
    // Otherwise, returns 0.
    // Container must be a random access container (e. g. vector, span, or deque).
    template<class Container> uint64_t isCopyNumberDifference(
        const Container& x,
        const Container& y,
        uint64_t maxPeriod);

}



template<class Container> uint64_t shasta::isCopyNumberDifference(
    const Container& x,
    const Container& y,
    uint64_t maxPeriod)
{

    // Get the lengths and their differences.
    const uint64_t nx = x.size();
    const uint64_t ny = y.size();

    // When they have the same length return 0.
    if(nx == ny) {
        return 0;
    }

    // Recursive call so x is shorter than y.
    if(ny < nx) {
        return isCopyNumberDifference(y, x, maxPeriod);
    }
    SHASTA_ASSERT(nx < ny);

    // If the length difference is not a multiple of one of the allowed periods,
    // return 0.
    const uint64_t dn = ny - nx;
    bool found = false;
    for(uint64_t period=2; period<=maxPeriod; period++) {
        if((dn % period) == 0) {
            found = true;
            break;
        }
    }
    if(not found) {
        return 0;
    }


    const uint64_t prefixLength = commonPrefixLength(x, y);
    const uint64_t suffixLength = commonSuffixLength(x, y);

    // Find the portion of y that is not in x.
    uint64_t ix = prefixLength;
    uint64_t iy = prefixLength;
    uint64_t jx = nx - suffixLength;
    uint64_t jy = ny - suffixLength;

    // Reduce the suffix side to avoid overlap between the prefixes and suffixes.
    while((jx<ix) or (jy<iy)) {
        ++jx;
        ++jy;
    }

    if(ix != jx) {
        // There is more than just an insertion.
        return 0;
    }


    // If getting here, x and y differ by an insertion in iy of range [iy, jy).
    SHASTA_ASSERT(ix == jx);
    SHASTA_ASSERT(jy - iy == dn);



    // Check for k base repeat.
    // We kept the entire common prefix, so we can check just to the left of the insertion.
    for(uint64_t period=2; period<=maxPeriod; period++) {
        if((dn % period) != 0) {
            continue;
        }

        // Check that the inserted bases are a repeat with this period.
        const uint64_t m = dn / period;
        bool repeatViolationFound = false;
        for(uint64_t i=0; i<m; i++) {
            for(uint64_t j=0; j<period; j++) {
                if(y[iy + i*period + j] != y[iy + j]) {
                    repeatViolationFound = true;
                    break;
                }
            }
        }
        if(repeatViolationFound) {
            // The inserted portion is not a repeat with this period.
            continue;
        }

        // Check the previous period bases in both x and y.
        if(ix < period) {
            continue;
        }
        if(iy < period) {
            continue;
        }
        for(uint64_t j=0; j<period; j++) {
            if(y[iy - period + j] != y[ix + j]) {
                repeatViolationFound = true;
                break;
            }
            if(x[ix - period + j] != y[ix + j]) {
                repeatViolationFound = true;
                break;
            }
        }
        if(repeatViolationFound) {
            continue;
        }

        // It is an insertion in y with this period.
        return period;
    }



    // None of the periods we tried worked.
    return 0;

}

#endif

