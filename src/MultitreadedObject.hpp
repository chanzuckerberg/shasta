#ifndef CZI_SHASTA_MULTITHREADED_OBJECT_HPP
#define CZI_SHASTA_MULTITHREADED_OBJECT_HPP



// Template class MultithreadedObject can be used as a base class
// to provide basic multithreading functionality.

// Usage pattern:
// class A : public MultithreadedObject<A> {
// public:
//     A() : MultithreadedObject(*this) {}
//     void compute(size_t threadId);
// };

// CZI.
#include "CZI_ASSERT.hpp"

// Standard libraries.
#include "algorithm.hpp"
#include "cstddef.hpp"
#include "fstream.hpp"
#include "iostream.hpp"
#include <memory>
#include <mutex>
#include "stdexcept.hpp"
#include "string.hpp"
#include <thread>
#include "utility.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        template<class T> class MultithreadedObject;
        void testMultithreadedObject();
        class MultithreadedObjectTestClass;
    }
}



template<class T> class ChanZuckerberg::shasta::MultithreadedObject {
public:

    // A function passed as argument to runThreads or startThreads
    // takes as input an integer thread id.
    // The function is called asynchronously for all integer
    // thread ids in [0, threadCount-1].
    using ThreadFunction = void (T::*)(size_t threadId);

    // Start the threads and wait for them to complete.
    void runThreads(
        ThreadFunction,
        size_t threadCount,
        const string& logFileNamePrefix = "");

    // Start the threads without waiting for them to complete.
    void startThreads(
        ThreadFunction,
        size_t threadCount,
        const string& logFileNamePrefix = "");

    // Wait for the running threads to complete.
    void waitForThreads();

    // Dynamic load balancing.
    void setupLoadBalancing(
        size_t n,
        size_t batchSize);

protected:

    // The constructor stores a reference to *this.
    MultithreadedObject(T&);

    bool getNextBatch(
        size_t& begin,
        size_t& end);

    ostream& getLog(size_t threadId)
    {
        CZI_ASSERT(threadId < threadLogs.size());
        ofstream& s = threadLogs[threadId];
        if(!s.is_open()) {
            throw runtime_error("Attempt to write to unopened log output for thread " + to_string(threadId));
        }
        return s;
    }

    // General purpose mutex, used when exclusive access is needed.
    // it is better to use atomic memort primitive for synchronizationinstead,
    // whenever possible.
    std::mutex mutex;

private:
    T& t;

    static void runThreadFunction(T& t, ThreadFunction f, size_t threadId)
    {
        try {
            (t.*f)(threadId);
        } catch(exception& e) {
            t.exceptionsOccurred = true;
            std::lock_guard<std::mutex> lock(t.mutex);
            cout << "A standard exception occurred in thread " << threadId << ": ";
            cout << e.what() << endl;
        } catch(...) {
            t.exceptionsOccurred = true;
            std::lock_guard<std::mutex> lock(t.mutex);
            cout << "A non-standard exception occurred in thread " << threadId << "." << endl;
        }
    }

    vector< std::shared_ptr<std::thread> > threads;

    vector<ofstream> threadLogs;

    bool exceptionsOccurred= false;

    // Load balancing.
    size_t n = 0;
    size_t batchSize = 0;
    size_t nextBatch = 0;
};



template<class T> inline ChanZuckerberg::shasta::MultithreadedObject<T>::MultithreadedObject(T& t) :
    t(t)
{
}



template<class T> inline void ChanZuckerberg::shasta::MultithreadedObject<T>::runThreads(
    ThreadFunction f,
    size_t threadCount,
    const string& logFileNamePrefix)
{
    startThreads(f, threadCount, logFileNamePrefix);
    waitForThreads();
}



template<class T> inline void ChanZuckerberg::shasta::MultithreadedObject<T>::startThreads(
    ThreadFunction f,
    size_t threadCount,
    const string& logFileNamePrefix)
{
    if(!threads.empty()) {
        throw runtime_error("Unsupported attempt to start new threads while other threads have not been joined.");
    }
    CZI_ASSERT(threadLogs.empty());

    // __sync_synchronize (); A full memory barrier is probably not needed here.
    exceptionsOccurred = false;
    threadLogs.resize(threadCount);
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        if(!logFileNamePrefix.empty()) {
            auto& log = threadLogs[threadId];
            const string fileName = logFileNamePrefix + "-" + to_string(threadId);
            log.open(fileName);
            if(!log) {
                throw runtime_error("Error opening thread log file " + fileName);
            }
            log.exceptions(ofstream::failbit | ofstream::badbit );
        }
        threads.push_back(std::make_shared<std::thread>(
            std::thread(
            &MultithreadedObject::runThreadFunction,
            std::ref(t),
            f,
            threadId)));
    }
}



template<class T> inline void ChanZuckerberg::shasta::MultithreadedObject<T>::waitForThreads()
{
    for(std::shared_ptr<std::thread> thread: threads) {
        thread->join();
    }
    threads.clear();
    threadLogs.clear();
    if(exceptionsOccurred) {
        throw runtime_error("Exceptions occurred in at least one thread.");
    }
    exceptionsOccurred = false;
    // __sync_synchronize (); A full memory barrier is probably not needed here.
}



template<class T> inline void ChanZuckerberg::shasta::MultithreadedObject<T>::setupLoadBalancing(
    size_t nArgument,
    size_t batchSizeArgument)
{
    n = nArgument;
    batchSize = batchSizeArgument;
    nextBatch = 0;
}
template<class T> inline bool ChanZuckerberg::shasta::MultithreadedObject<T>:: getNextBatch(
    size_t& begin,
    size_t& end)
{
    begin = __sync_fetch_and_add(&nextBatch, batchSize);
    if(begin < n) {
        end = min(n, begin + batchSize);
        return true;
    } else {
        return false;
    }
}

#endif
