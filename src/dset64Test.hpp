#ifndef SHASTA_DSET_64_TEST_HPP
#define SHASTA_DSET_64_TEST_HPP

// Unit test for dset64.hpp/dset64-gccAtomic.hpp.
#include "dset64-gccAtomic.hpp"
#include "MultithreadedObject.hpp"
#include "utility.hpp"
#include <map>

namespace shasta {
    void dset64Test(
        uint64_t n,             // The number of items (vertices).
        uint64_t m,             // The number of union operations (edges).
        uint64_t threadCount,   // The number of threads to use.
        uint64_t batchSize,     // The number of union operations per batch.
        int seed                // The random seed.
        );
    class Dset64Test;

    extern template class MultithreadedObject<Dset64Test>;
}


// Class describing the overlap between a pair of oriented reads.
class shasta::Dset64Test : public MultithreadedObject<Dset64Test> {
public:

    Dset64Test(
        uint64_t n,             // The number of items (vertices).
        uint64_t m,             // The number of union operations (edges).
        uint64_t threadCount,   // The number of threads to use.
        uint64_t batchSize,     // The number of union operations per batch.
        int seed                // The random seed.
        );

private:

    vector< pair<uint64_t, uint64_t> > edges;

    static void getSortedComponents(
        const std::map<uint64_t, vector<uint64_t> >& componentTable,
        vector< vector<uint64_t> >& sortedComponents
    );

    DisjointSets* disjointSetsPointer;
    void threadFunction(size_t threadId);
};

#endif
