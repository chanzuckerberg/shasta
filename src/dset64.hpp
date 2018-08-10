#if !defined(__DSET64_HPP)
#define __DSET64_HPP

#include <atomic>
#include <stdexcept>

/**
 * Lock-free parallel disjoint set data structure (aka UNION-FIND)
 * with path compression and union by rank
 *
 * Supports concurrent find(), same() and unite() calls as described
 * in the paper
 *
 * "Wait-free Parallel Algorithms for the Union-Find Problem"
 * by Richard J. Anderson and Heather Woll
 *
 * This file is a modified version of dset.h from GitHub repository
 * wjakob/dset by Wenzel Jacob.
 *
 * The original implementation by Wenzel Jacob uses 64-bit
 * atomic primitives, and implements the union-find
 * algorithm for 32-bit item ids, which allows up to
 * 2^32^ items. This modified version
 * uses 128-bit primitives for 64-bit item ids,
 * which brings the maximum number of items to 2^64^.
 *
 * See the LICENSE file for licensing information
 * specific to this file.
 *
 * \author Wenzel Jakob
 *
 *
 * IMPORTANT COMPILATION AND PORTABILITY INFORMATION
 *
 * 128-bit integers and 128-bit atomic primitives
 * are not available on all platforms.
 * This code uses __uint128_t as the 128 integer type,
 * and synchronization primitives available in C++11
 * via std::atomic<uint128_t>
 *
 * Both __uint128_t and std::atomic<__uint128_t>
 * are available with g++ on 64-bit x86 Linux,
 * and std::atomic<uint128_t>::is_lock_free returns true,
 * if the following compilation option is used to enable the
 * 16-byte (128 bit) compare-and-swap instruction, CMPXCHG16B:
 * -mcx16
 * This instruction is available on most modern 64-bit x86 processors.
 * Some older processors that don't implement this instruction
 * will crash with an "Illegal instruction" error
 * upon attempting to run this code.
 *
 * Unfortunately, however, beginning with gcc 7, a bug was introduced
 * that causes std::atomic<uint128_t>::is_lock_free to return false,
 * even when compile option -mcx16 is used (gcc bug 80878):
 * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=80878
 *
 * Ubuntu 18.04 uses gcc 7.3, so here are the alternatives
 * to port to Ubuntu 18.04:
 * - The gcc bug gets fixed.
 * - We change the code to use gcc __sync primitives instead of std::atomic.
 * - We use an older version of gcc on ubuntu 18.04 (messy).
 *
 */


class DisjointSets {
public:

    // Integer type used for the item ids.
    // This determines the maximum number of items that can be handled.
    using Uint = uint64_t;
    static_assert(sizeof(Uint) == 8, "Unexpected size of DisjointSets::Uint.");

    // Integer type used for synchronization primitives.
    // This must have 128 bits and its atomic type must be lock-free
    // (this is checked in the constructor).
    // See the Compilation/Portability comment above.
    using Aint = __uint128_t;
    static_assert(sizeof(Aint) == 16, "Unexpected size of DisjointSets::Aint.");

    // We use the 128 bits of Aint to hold the parent in the
    // 64 least significant bits and the rank in the most significant bits.
    // Define two masks for these two sets of bits.
    static const Aint parentMask = Aint(~0ULL);
    static const Aint rankMask = parentMask << 64;

    // For memory allocation flexibility, the memory is allocated
    // and ownded by the caller.
    DisjointSets(std::atomic<Aint>* mData, Uint size) : mData(mData), n(size) {
        if(!mData->is_lock_free()) {
            // If this happens with g++ on 64-bit x86 Linux, use
            // compile option -mcx16.
            // This throw is commented out to permit running on Ubuntu 18.04
            // (see comments above), but at a performance penalty,
            // possibly significant.
            // throw std::runtime_error("DisjointSets::Aint is not lock-free.");
        }
        for (Uint i=0; i<size; ++i)
            mData[i] = Aint(i);
    }

    Uint find(Uint id) const {
        while (id != parent(id)) {
            Aint value = mData[id];
            Uint new_parent = parent((Uint) value);
            Aint new_value =
                (value & rankMask) | new_parent;
            /* Try to update parent (may fail, that's ok) */
            if (value != new_value)
                mData[id].compare_exchange_weak(value, new_value);
            id = new_parent;
        }
        return id;
    }

    bool same(Uint id1, Uint id2) const {
        for (;;) {
            id1 = find(id1);
            id2 = find(id2);
            if (id1 == id2)
                return true;
            if (parent(id1) == id1)
                return false;
        }
    }

    Uint unite(Uint id1, Uint id2) {
        for (;;) {
            id1 = find(id1);
            id2 = find(id2);

            if (id1 == id2)
                return id1;

            Uint r1 = rank(id1), r2 = rank(id2);

            if (r1 > r2 || (r1 == r2 && id1 < id2)) {
                std::swap(r1, r2);
                std::swap(id1, id2);
            }

            Aint oldEntry = ((Aint) r1 << 64) | id1;
            Aint newEntry = ((Aint) r1 << 64) | id2;

            if (!mData[id1].compare_exchange_strong(oldEntry, newEntry))
                continue;

            if (r1 == r2) {
                oldEntry = ((Aint) r2 << 64) | id2;
                newEntry = ((Aint) (r2+1) << 64) | id2;
                /* Try to update the rank (may fail, that's ok) */
                mData[id2].compare_exchange_weak(oldEntry, newEntry);
            }

            break;
        }
        return id2;
    }

    Uint size() const { return n; }

    Uint rank(Uint id) const {
        return ((Uint) (mData[id] >> 64)) & parentMask;
    }

    Uint parent(Uint id) const {
        return (Uint) mData[id];
    }

    // Use memory supplied by the caller, rather than an owned vector.
    // This provides more flexibility in allocating the memory.
    std::atomic<Aint>* mData;
    Uint n;
};

#endif /* __DSET64_HPP */
