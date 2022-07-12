#ifndef SHASTA_MULTITHREADED_OBJECT_HPP
#define SHASTA_MULTITHREADED_OBJECT_HPP



// Template class MultithreadedObject can be used as a base class
// to provide basic multithreading functionality.

// Usage pattern:
// class A : public MultithreadedObject<A> {
// public:
//     A() : MultithreadedObject(*this) {}
// };

// Standard libraries.
#include <memory>
#include <mutex>
#include <thread>
#include "vector.hpp"

namespace shasta {
    class MultithreadedObjectBaseClass;
    template<class T> class MultithreadedObject;

    // Testing.
    void testMultithreadedObject();
    class MultithreadedObjectTestClass;
}



// The base class contains code that is not templated.
class shasta::MultithreadedObjectBaseClass {
protected:

    // The running threads.
    vector< std::shared_ptr<std::thread> > threads;

    // General purpose mutex, used when exclusive access is needed.
    // It is better to use atomic memory primitives for synchronization instead,
    // whenever possible.
    // Make it mutable so it can also be used by const functions.
    mutable std::mutex mutex;

    // Keep track of exceptions occurring in threads.
    bool exceptionsOccurred= false;

    // Dynamic load balancing.
    void setupLoadBalancing(
        uint64_t n,
        uint64_t batchSize);
    bool getNextBatch(
        uint64_t& begin,
        uint64_t& end);

    // Wait for the running threads to complete.
    void waitForThreads();

    // Kill all threads except the one passed as argument.
    void killAllThreadsExceptMe(size_t me);

private:
    uint64_t n = 0;
    uint64_t batchSize = 0;
    uint64_t nextBatch = 0;
};



template<class T> class shasta::MultithreadedObject : public MultithreadedObjectBaseClass {
public:

    // A function passed as argument to runThreads or startThreads
    // takes as input an integer thread id.
    // The function is called in parallel for all integer
    // thread ids in [0, threadCount-1].
    using ThreadFunction = void (T::*)(size_t threadId);

    // Start the threads and wait for them to complete.
    void runThreads(
        ThreadFunction,
        size_t threadCount);

    // Start the threads without waiting for them to complete.
    void startThreads(
        ThreadFunction,
        size_t threadCount);

protected:

    // The constructor stores a reference to *this.
    MultithreadedObject(T&);

private:
    T& t;

    // The function run by each thread.
    static void runThreadFunction(T& t, ThreadFunction f, size_t threadId);
};



#endif
