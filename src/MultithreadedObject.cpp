// Shasta.
#include "MultithreadedObject.hpp"
#include "MultithreadedObject.tpp"
#include "timestamp.hpp"
using namespace shasta;

// Linux.
#include <pthread.h>

// Standard library.
#include "chrono.hpp"



// Class used only by function testMultithreadedObject.
class shasta::MultithreadedObjectTestClass :
    public MultithreadedObject<MultithreadedObjectTestClass> {
public:

    MultithreadedObjectTestClass(size_t n) : MultithreadedObject(*this), n(n)
    {
        x.resize(n);
        y.resize(n);
        z.resize(n);
        for(uint64_t i=0; i<n; i++) {
            x[i] = i;
            y[i] = 2 * i;
            z[i] = 18;
        }
    }
    void compute(size_t threadId)
    {

        uint64_t begin, end;
        while(getNextBatch(begin, end)) {
            for(uint64_t i=begin; i!=end; i++) {
                uint64_t s = 0;
                for(uint64_t j=0; j<n; j++) {
                    s += x[i] * y[j];
                }
                z[i] = s;
            }
        }
    }
    void check() const
    {
        for(uint64_t i=0; i<x.size(); i++) {
            uint64_t s = 0;
            for(uint64_t j=0; j<n; j++) {
                s += x[i] * y[j];
            }
            SHASTA_ASSERT(z[i] == s);
        }
    }
    void run(
        uint64_t n,
        uint64_t batchSize,
        size_t threadCount
        )
    {
        setupLoadBalancing(n, batchSize);
        runThreads(&MultithreadedObjectTestClass::compute, threadCount);
    }

    uint64_t n;
    vector<uint64_t> x;
    vector<uint64_t> y;
    vector<uint64_t> z;
};



void shasta::MultithreadedObjectBaseClass::waitForThreads()
{
    for(std::shared_ptr<std::thread> thread: threads) {
        thread->join();
    }
    threads.clear();
    if(exceptionsOccurred) {
        throw runtime_error("Exceptions occurred in at least one thread.");
    }
    exceptionsOccurred = false;
    // __sync_synchronize (); A full memory barrier is probably not needed here.
}



void shasta::MultithreadedObjectBaseClass::setupLoadBalancing(
    uint64_t nArgument,
    uint64_t batchSizeArgument)
{
    n = nArgument;
    batchSize = batchSizeArgument;
    nextBatch = 0;
}



bool shasta::MultithreadedObjectBaseClass::getNextBatch(
    uint64_t& begin,
    uint64_t& end)
{
    begin = __sync_fetch_and_add(&nextBatch, batchSize);
    if(begin < n) {
        end = min(n, begin + batchSize);
        return true;
    } else {
        return false;
    }
}



void shasta::MultithreadedObjectBaseClass::killAllThreadsExceptMe(size_t me)
{
    for(size_t threadId=0; threadId<threads.size(); threadId++) {
        if(threadId == me) {
            continue;
        }
        const std::shared_ptr<std::thread> thread = threads[threadId];
        const auto h = thread->native_handle();
        if(h) {
            ::pthread_cancel(h);
        }
    }
}


void shasta::testMultithreadedObject()
{

    const uint64_t n = 32 * 1024;
    const uint64_t batchSize = 64;
    const uint64_t threadCount = 8;
    MultithreadedObjectTestClass x(n);
    for(int i=0; i<10; i++) {
        const auto t0 = std::chrono::steady_clock::now();
        x.run(n, batchSize, threadCount);
        const auto t1 = std::chrono::steady_clock::now();
        const double t01 = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count());
        cout << t01 << endl;
    }
    x.check();
}

