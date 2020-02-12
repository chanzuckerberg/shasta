#ifndef SHASTA_COMPRESS_ALIGNMENT_HPP
#define SHASTA_COMPRESS_ALIGNMENT_HPP

/*******************************************************************************

Compress/decompress a marker alignment to bytes.

An alignment is a vector of pairs of marker ordinals, not necessarily sorted
in any partricular way. However, in most cases each pair differs from
the previous pair by a smal increment of most ordinals.

The increment is usually positive and often equal to 1.

The functions defined here implement an ad hoc compression scheme for
marker alignments that takes advantage of these characteristics.

A marker alignment is decomposed as a sequence of streaks.
In each streak, each pair of ordinals can be obtained from the
previous pair by simply incrementing both ordinals by 1.

A streak can be completely described by the number of marker pairs
in the streak, n, plus the number of markers skipped relative to the
previous streak, skip0 and skip1 (or relative to the origin, for the first streak).
Note that n is always unsigned, but skipx and skipy can in general
be negative, and so they generally need to be represented using
a signed integer.

For example, consider the following alignment, which consists of three streaks:

300 200 The first streak begins here
301 201
302 202
305 206 The second streak begins here
306 207
320 250 The third streak begins here
321 251
322 252
323 253

This can be described by the following tuples (skip0, skip1, n):

300 200 3
  3   4 2
 14  43 4

The compressed format stores a sequence of such tuples,
each describing a streak. For space economy, each tuple can be stored
in a variable number of bytes according to a small number of
possible formats described below. Each streak is described using
the small format that can be used for the streak.

The least significant few bits of the first byte of each streak identify
the format used.

The table below summarizes the formats used:

Format                                     0       1       2       3       4
Streak size (bytes)                        1       2       4       8      16
Streak size (bits)                         8      16      32      64     128
Least significant bits of first byte,
used to identify the format                0     001     011     101     111
Number of bits not used to identify
the format                                 7      13      29      61     125
Number of bits used to represent n-1       3       5       9      21      32
Number of bits used to represent
each of skipx and skipy                    2       4      10      20      32
skipx and skipy are signed                 No      Yes    Yes     Yes     Yes
Minimum value of skipx and skipy
that can be represented.                   0      -8     -512    2^19-1  2^31-1
Maximum value of skipx and skipy
that can be represented                    3       7      511   -2^19   -2^31

*******************************************************************************/

// Shasta.
#include "Alignment.hpp"

// Standard library.
#include "string.hpp"
#include "span.hpp"


namespace shasta {
    void compress(const Alignment&, string&);
    void decompress(span<const char>, Alignment);

    namespace compressAlignment {
        class Format0;
        class Format1;
        class Format2;
        class Format3;
        class Format4;
    }
}



class shasta::compressAlignment::Format0 {
public:
    uint8_t formatIdentifier: 1;
    uint8_t skipx: 2;
    uint8_t skipy: 2;
    uint8_t nMinus1: 3;
    uint32_t n() const
    {
        return uint32_t(nMinus1) + 1;
    }
    Format0(
        int32_t skipxArgument,
        int32_t skipyArgument,
        uint32_t nArgument)
    {
        SHASTA_ASSERT(skipxArgument >= 0);
        SHASTA_ASSERT(skipxArgument <= 3);
        SHASTA_ASSERT(skipyArgument >= 0);
        SHASTA_ASSERT(skipyArgument <= 3);
        SHASTA_ASSERT(nArgument <= 8);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
        skipx = skipxArgument;
        skipy = skipyArgument;
        nMinus1 = nArgument - 1;
#pragma GCC diagnostic pop
    }
};
static_assert(sizeof(shasta::compressAlignment::Format0) == 1,
    "Unexpected size for shasta::compressAlignment::Format0");



class shasta::compressAlignment::Format1 {
public:
    uint16_t formatIdentifier: 3;
    int16_t skipx: 4;
    int16_t skipy: 4;
    uint16_t nMinus1: 5;
    uint32_t n() const
    {
        return uint32_t(nMinus1) + 1;
    }
    Format1(
         int32_t skipxArgument,
         int32_t skipyArgument,
         uint32_t nArgument)
     {
         SHASTA_ASSERT(skipxArgument >= -8);
         SHASTA_ASSERT(skipxArgument <= 7);
         SHASTA_ASSERT(skipyArgument >= -8);
         SHASTA_ASSERT(skipyArgument <= 7);
         SHASTA_ASSERT(nArgument <= 32);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
         skipx = skipxArgument;
         skipy = skipyArgument;
         nMinus1 = nArgument - 1;
#pragma GCC diagnostic pop
     }
};
static_assert(sizeof(shasta::compressAlignment::Format1) == 2,
    "Unexpected size for shasta::compressAlignment::Format1");



class shasta::compressAlignment::Format2 {
public:
    uint32_t formatIdentifier: 3;
    int32_t skipx: 10;
    int32_t skipy: 10;
    uint32_t nMinus1: 9;
    uint32_t n() const
    {
        return nMinus1 + 1;
    }
    Format2(
         int32_t skipxArgument,
         int32_t skipyArgument,
         uint32_t nArgument)
     {
         SHASTA_ASSERT(skipxArgument >= -512);
         SHASTA_ASSERT(skipxArgument <= 511);
         SHASTA_ASSERT(skipyArgument >= -512);
         SHASTA_ASSERT(skipyArgument <= 511);
         SHASTA_ASSERT(nArgument <= 512);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
         skipx = skipxArgument;
         skipy = skipyArgument;
         nMinus1 = nArgument - 1;
#pragma GCC diagnostic pop
     }
};
static_assert(sizeof(shasta::compressAlignment::Format2) == 4,
    "Unexpected size for shasta::compressAlignment::Format2");



class shasta::compressAlignment::Format3 {
public:
    uint64_t formatIdentifier: 3;
    int64_t skipx: 20;
    int64_t skipy: 20;
    uint64_t nMinus1: 21;
    uint32_t n() const
    {
        return uint32_t(nMinus1) + 1;
    }
    Format3(
         int32_t skipxArgument,
         int32_t skipyArgument,
         uint32_t nArgument)
     {
         SHASTA_ASSERT(skipxArgument >= -524288);
         SHASTA_ASSERT(skipxArgument <= 524287);
         SHASTA_ASSERT(skipyArgument >= -524288);
         SHASTA_ASSERT(skipyArgument <= 524287);
         SHASTA_ASSERT(nArgument <= 2097152);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
         skipx = skipxArgument;
         skipy = skipyArgument;
         nMinus1 = nArgument - 1;
#pragma GCC diagnostic pop
     }
};
static_assert(sizeof(shasta::compressAlignment::Format3) == 8,
    "Unexpected size for shasta::compressAlignment::Format3");



class shasta::compressAlignment::Format4 {
public:
    uint32_t formatIdentifier: 3;
    int32_t skipx;
    int32_t skipy;
    uint32_t nMinus1;
    uint32_t n() const
    {
        return nMinus1 + 1;
    }
    Format4(
         int32_t skipx,
         int32_t skipy,
         uint32_t n):
         skipx(skipx),
         skipy(skipy),
         nMinus1(n-1)
     {
     }
};
static_assert(sizeof(shasta::compressAlignment::Format4) == 16,
    "Unexpected size for shasta::compressAlignment::Format4");

#endif
