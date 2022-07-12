// Implementation of template class MultithreadedObject.

// Shasta.
#include "MultithreadedObject.hpp"
#include "SHASTA_ASSERT.hpp"
#include "timestamp.hpp"

// Standard library.
#include "algorithm.hpp"
#include "iostream.hpp"
#include "stdexcept.hpp"
#include "string.hpp"



template<class T> shasta::MultithreadedObject<T>::MultithreadedObject(T& t) :
    t(t)
{
}



template<class T> void shasta::MultithreadedObject<T>::runThreads(
    ThreadFunction f,
    size_t threadCount)
{
    startThreads(f, threadCount);
    waitForThreads();
}



template<class T> void shasta::MultithreadedObject<T>::startThreads(
    ThreadFunction f,
    size_t threadCount)
{
    SHASTA_ASSERT(threadCount > 0);

    if(!threads.empty()) {
        throw runtime_error("Unsupported attempt to start new threads while other threads have not been joined.");
    }

    // __sync_synchronize (); A full memory barrier is probably not needed here.
    exceptionsOccurred = false;
    for(size_t threadId=0; threadId<threadCount; threadId++) {
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
template<class T> void shasta::MultithreadedObject<T>::runThreadFunction(
    T& t,
    ThreadFunction f,
    size_t threadId)
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

