// Nanopore2.
#include "dset64Test.hpp"
#include "vector.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard libraries.
#include "chrono.hpp"
#include "iostream.hpp"
#include <random>



void ChanZuckerberg::Nanopore2::dset64Test(
    uint64_t n,             // The number of items (vertices).
    uint64_t m,             // The number of union operations (edges).
    uint64_t threadCount,   // The number of threads to use.
    uint64_t batchSize,     // The number of union operations per batch.
    int seed                // The random seed.
    )
{
    Dset64Test(n, m, threadCount, batchSize, seed);
}



Dset64Test::Dset64Test(
    uint64_t n,             // The number of items (vertices).
    uint64_t m,             // The number of union operations (edges).
    uint64_t threadCount,   // The number of threads to use.
    uint64_t batchSize,     // The number of union operations per batch.
    int seed                // The random seed.
    ) : MultithreadedObject(*this)
{

    // Randomly generate the union operations.
    // Prepare to generate uniformly distributed numbers between 0 and 1.
    std::mt19937 randomSource(seed);
    std::uniform_int_distribution<uint64_t> uniformDistribution(0, n-1);
    edges.resize(m);
    for(auto& p: edges) {
        p.first = uniformDistribution(randomSource);
        p.second = uniformDistribution(randomSource);
    }



    // First, do it sequentially using the disjoint sets implementation
    // in the Boost library.
    vector< vector<uint64_t> > sortedComponentsBoost;
    {
        vector<uint64_t> rank(n);
        vector<uint64_t> parent(n);
        boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
        for(uint64_t i=0; i<n; i++) {
            disjointSets.make_set(i);
        }
        const auto t0 = std::chrono::steady_clock::now();
        for(const auto& p: edges) {
            disjointSets.union_set(p.first, p.second);
        }
        const auto t1 = std::chrono::steady_clock::now();
        cout << "Boost implementation ran in " << seconds(t1-t0) << "s." << endl;

        // Gather the components.
        std::map<uint64_t, vector<uint64_t> > componentTable;
        for(uint64_t i=0; i<n; i++) {
            componentTable[disjointSets.find_set(i)].push_back(i);
        }
        getSortedComponents(componentTable, sortedComponentsBoost);
    }



    // Now, do it using dset64.hpp, sequentially.
    vector< vector<uint64_t> > sortedComponentsSequential;
    {
        using Aint = DisjointSets::Aint;
        vector< std::atomic<Aint> > data(n);
        DisjointSets disjointSets(&data.front(), n);
        const auto t0 = std::chrono::steady_clock::now();
        for(const auto& p: edges) {
            disjointSets.unite(p.first, p.second);
        }
        const auto t1 = std::chrono::steady_clock::now();
        cout << "Sequential dset64 ran in " << seconds(t1-t0) << "s." << endl;

        // Gather the components.
        std::map<uint64_t, vector<uint64_t> > componentTable;
        for(uint64_t i=0; i<n; i++) {
            componentTable[disjointSets.find(i)].push_back(i);
        }
        getSortedComponents(componentTable, sortedComponentsSequential);
    }
    CZI_ASSERT(sortedComponentsSequential == sortedComponentsBoost);



    // Now, do it using dset64.hpp, using the specified number of threads.
    vector< vector<uint64_t> > sortedComponentsParallel;
    {
        using Aint = DisjointSets::Aint;
        vector< std::atomic<Aint> > data(n);
        DisjointSets disjointSets(&data.front(), n);
        disjointSetsPointer = &disjointSets;
        const auto t0 = std::chrono::steady_clock::now();
        setupLoadBalancing(n, batchSize);
        runThreads(&Dset64Test::threadFunction, threadCount);
        const auto t1 = std::chrono::steady_clock::now();
        cout << "Parallel dset64 ran in " << seconds(t1-t0) << "s." << endl;

        // Gather the components.
        std::map<uint64_t, vector<uint64_t> > componentTable;
        for(uint64_t i=0; i<n; i++) {
            componentTable[disjointSets.find(i)].push_back(i);
        }
        getSortedComponents(componentTable, sortedComponentsParallel);
    }
    CZI_ASSERT(sortedComponentsParallel == sortedComponentsBoost);


    cout << "No error found. All algorithms found " << sortedComponentsBoost.size();
    cout << " identical connected components." << endl;
}



void Dset64Test::threadFunction(size_t threadId)
{
    uint64_t begin;
    uint64_t end;
    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; ++i) {
            const auto& p = edges[i];
            disjointSetsPointer->unite(p.first, p.second);
        }
    }

}



void Dset64Test::getSortedComponents(
    const std::map<uint64_t, vector<uint64_t> >& componentTable,
    vector< vector<uint64_t> >& sortedComponents
)
{
    sortedComponents.clear();
    for(const auto& p: componentTable) {
        sortedComponents.push_back(p.second);
    }
    sort(sortedComponents.begin(), sortedComponents.end());
}

