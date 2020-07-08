#if !defined(__DSET64_GCC_ATOMIC_HPP)
#define __DSET64_GCC_ATOMIC_HPP

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
 * The original implementation by Wenzel Jakob uses 64-bit
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
 * USAGE OF GCC ATOMIC PRIMITIVES
 *
 * The implementation in shasta/src/dset64.hpp uses std::atomic<__uint128_t>
 * for lock-free synchronization.
 * On older GCC versions, std::atomic<__uint128_t> is lock-free
 * if compilation is done with -mcx16, which enables the use of the
 * 16-byte (128 bit) compare-and-swap instruction, CMPXCHG16B.
 *
 * Unfortunately, on newer GCC versions, this is no longer true
 * because of gcc bug 80878:
 * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=80878
 *
 * As a result, there was a significant performance loss in
 * versions of Shasta built with gcc 7,
 * which is used by default on Ubuntu 18.04, when using
 * machines with large number of virtual processors.
 *
 * It is unlikely that this gcc bug will ever be fixed,
 * and to avoid this performance loss this implementation
 * uses gcc primitive __sync_bool_compare_and_swap instead
 * for lock-free synchronization. When compilation
 * is done with -mcx16 and optimization turned on,
 * this primitive uses the CMPXCHG16B instruction
 * and results in optimal speed.
 *
 * The CMPXCHG16B instruction is available on most modern 64-bit x86 processors.
 * If it is not available, and compilation is done without -mcx16, then gcc
 * generates slower code. Use of the gcc primitive __sync_bool_compare_and_swap
 * ensures that the code will execute on both x86_64 and aarch64.
 *
 */

// Sanity check that we are compiling on x86_64 or aarch64
#if !__x86_64__ && !__aarch64__
#error "Shasta can only be built on an x86_64 machine (64-bit Intel/AMD) or an ARM64 machine. "
#endif


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
    // and owned by the caller.
    DisjointSets(Aint* mData, Uint size) : mData(mData), n(size), parentUpdated(0) {
        for (Uint i=0; i<size; ++i)
            mData[i] = Aint(i);
    }


    // `find` is guaranteed to return the correct set representative. It is not guaranteed
    // that the parent information (lower 64 bits) will have the set representative stored in it.
    // This is because the `__sync_bool_compare_and_swap` could fail because of race conditions.
    // When `trackParentUpdated` is set to True, we bump a counter that tracks the number
    // of entries that don't have the set representative populated in their lower 64 bits.
    //
    // In a highly parallel environment, where N threads are calling `find` on different
    // entries in DisjointSets, it is very likely that several entries won't have the
    // set representative populated as the parent, at the end of the pass.
    //
    // If we do multiple such passes over all entries, then within a few iterations, the parent
    // information in all entries (lower 64 bits) will have the correct set representative.
    // `this->parentUpdated` can be used to track this convergence. 
    Uint find(Uint id, bool trackParentUpdated = false) {
        uint64_t parentUpdatedCount = 0;
        while (id != parent(id)) {
            Aint value = mData[id];
            Uint new_parent = parent((Uint) value);
            Aint new_value =
                (value & rankMask) | new_parent;
            /* Try to update parent (may fail, that's ok) */
            if (value != new_value) {
                bool swapped = __sync_bool_compare_and_swap(&mData[id], value, new_value);
                if (trackParentUpdated && swapped)
                    parentUpdatedCount++;
            }
            id = new_parent;
        }
        if (trackParentUpdated && parentUpdatedCount > 0) {
            __sync_fetch_and_add(&parentUpdated, parentUpdatedCount);
        }

        return id;
    }

    bool same(Uint id1, Uint id2) {
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

            if (!__sync_bool_compare_and_swap(&mData[id1], oldEntry, newEntry))
                continue;

            if (r1 == r2) {
                oldEntry = ((Aint) r2 << 64) | id2;
                newEntry = ((Aint) (r2+1) << 64) | id2;
                /* Try to update the rank (may fail, that's ok) */
                __sync_bool_compare_and_swap(&mData[id2], oldEntry, newEntry);
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
    Aint* mData;
    Uint n;

    uint64_t parentUpdated; // See comment on `find` method for an explanation.
};

#endif /* __DSET64_HPP */
