#ifndef SHASTA_MULTITHREADED_OBJECT_HPP
#define SHASTA_MULTITHREADED_OBJECT_HPP



// Template class MultithreadedObject can be used as a base class
// to provide basic multithreading functionality.

// Usage pattern:
// class A : public MultithreadedObject<A> {
// public:
//     A() : MultithreadedObject(*this) {}
//     void compute(size_t threadId);
// };

// Shasta.
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"

// Linux.
#include <pthread.h>

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

namespace shasta {
    template<class T> class MultithreadedObject;
    void testMultithreadedObject();
    class MultithreadedObjectTestClass;
}



template<class T> class shasta::MultithreadedObject {
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
        uint64_t n,
        uint64_t batchSize);

protected:

    // The constructor stores a reference to *this.
    MultithreadedObject(T&);

    bool getNextBatch(
        uint64_t& begin,
        uint64_t& end);

    ostream& getLog(size_t threadId)
    {
        SHASTA_ASSERT(threadId < threadLogs.size());
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



    // The function run by each thread.
    // Note if an exception happens in a thread,
    // we don't want to wait for all remaining threads to finish.
    // Therefore, the catch block calls exit.
    // A slightly cleaner termination may be possible
    // via pthread_cancel, but even that would result in
    // lack of destruction of objects created by the other threads.
    // There is also the possibility of using std::exception_ptr
    // to propagate the exception, but even that way a completely
    // clean termination is not possible without waiting
    // for all threads to finish.
    static void runThreadFunction(T& t, ThreadFunction f, size_t threadId)
    {
        try {
            (t.*f)(threadId);
        } catch(const runtime_error& e) {
            t.exceptionsOccurred = true;
            std::lock_guard<std::mutex> lock(t.mutex);
            cout << timestamp << "A runtime error occurred in thread " << threadId << ":\n";
            cout << e.what() << endl;
            t.killAllThreadsExceptMe(threadId);
            ::exit(1);
        } catch(const std::bad_alloc& e) {
            t.exceptionsOccurred = true;
            std::lock_guard<std::mutex> lock(t.mutex);
            cout << timestamp << e.what() << endl;
            cout << "A memory allocation failure occurred in thread " << threadId << ":\n";
            cout << "This assembly requires more memory than available." << endl;
            cout << "Rerun on a larger machine." << endl;
            t.killAllThreadsExceptMe(threadId);
            ::exit(1);
        } catch(const exception& e) {
            t.exceptionsOccurred = true;
            std::lock_guard<std::mutex> lock(t.mutex);
            cout << "A standard exception occurred in thread " << threadId << ":\n";
            cout << e.what() << endl;
            t.killAllThreadsExceptMe(threadId);
            ::exit(1);
        } catch(...) {
            t.exceptionsOccurred = true;
            std::lock_guard<std::mutex> lock(t.mutex);
            cout << "A non-standard exception occurred in thread " << threadId << "." << endl;
            t.killAllThreadsExceptMe(threadId);
            ::exit(1);
        }
    }



    vector< std::shared_ptr<std::thread> > threads;

    void killAllThreadsExceptMe(size_t me)
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

    vector<ofstream> threadLogs;

    bool exceptionsOccurred= false;

    // Load balancing.
    uint64_t n = 0;
    uint64_t batchSize = 0;
    uint64_t nextBatch = 0;
};



template<class T> inline shasta::MultithreadedObject<T>::MultithreadedObject(T& t) :
    t(t)
{
}



template<class T> inline void shasta::MultithreadedObject<T>::runThreads(
    ThreadFunction f,
    size_t threadCount,
    const string& logFileNamePrefix)
{
    startThreads(f, threadCount, logFileNamePrefix);
    waitForThreads();
}



template<class T> inline void shasta::MultithreadedObject<T>::startThreads(
    ThreadFunction f,
    size_t threadCount,
    const string& logFileNamePrefix)
{
    if(!threads.empty()) {
        throw runtime_error("Unsupported attempt to start new threads while other threads have not been joined.");
    }
    SHASTA_ASSERT(threadLogs.empty());

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

        try {
            threads.push_back(std::make_shared<std::thread>(
                std::thread(
                &MultithreadedObject::runThreadFunction,
                std::ref(t),
                f,
                threadId)));
        } catch(const std::exception& e) {
            throw runtime_error(
                "The following error occurred while attempting to start thread " +
                to_string(threadId) + ":\n" + e.what() + "\n" +
                "You may have hit a limit imposed by your system on the maximum number of threads "
                "allowed. Rerunning with \"--threads " + to_string(threadId) + "\" may fix this problem "
                "at a cost in performance.");
        }
    }
}



template<class T> inline void shasta::MultithreadedObject<T>::waitForThreads()
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



template<class T> inline void shasta::MultithreadedObject<T>::setupLoadBalancing(
    uint64_t nArgument,
    uint64_t batchSizeArgument)
{
    n = nArgument;
    batchSize = batchSizeArgument;
    nextBatch = 0;
}
template<class T> inline bool shasta::MultithreadedObject<T>:: getNextBatch(
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

#endif
