#include "dset64.hpp"
#include "dset64Test.hpp"
#include "vector.hpp"



void ChanZuckerberg::Nanopore2::dset64Test()
{
    // For the test, allocate the memory in a vactor.
    using Aint = DisjointSets::Aint;
    const uint64_t n = 10;
    vector< std::atomic<Aint> > data(n);

    // Construct the DisjointSets object.
    DisjointSets disjointSets(&data.front(), n);
}
