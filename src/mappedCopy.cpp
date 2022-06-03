
// shasta.
#include "mappedCopy.hpp"
#include "timestamp.hpp"

// Standard library.
#include "algorithm.hpp"
#include <chrono>
#include "iostream.hpp"
#include "stdexcept.hpp"

// Linux.
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// This can be used to copy a file to the huge page filesystem.
// The regular cp command does not work (but it works to copy
// the other way around, from the huge page filesystem).
void shasta::mappedCopy(
    const string& inputPath,
    const string& outputPath)
{
    cout << timestamp << "Copying " << inputPath << " to " << outputPath << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    // Open the input file.
    const int inputFileDescriptor = ::open(inputPath.c_str(), O_RDONLY);
    if(inputFileDescriptor == -1) {
        throw runtime_error("Error opening " + inputPath);
    }

    // Let the system know that we will be accessing this file sequentially.
    // This improves performance in some cases.
    posix_fadvise(inputFileDescriptor, 0, 0, POSIX_FADV_SEQUENTIAL);

    // Open the output file.
    const int outputFileDescriptor = ::open(outputPath.c_str(),
        O_CREAT | O_TRUNC | O_RDWR,
        S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if(outputFileDescriptor == -1) {
        ::close(inputFileDescriptor);
        throw runtime_error("Error opening " + outputPath);
    }

    // Get the size of the input file.
    const off_t n = ::lseek(inputFileDescriptor, 0, SEEK_END);
    if(n == -1) {
        ::close(inputFileDescriptor);
        ::close(outputFileDescriptor);
        throw runtime_error("Error during lseek for " + inputPath);
    }

    // Set the size of the output file.
    if(::truncate(outputPath.c_str(), n) == -1) {
        ::close(inputFileDescriptor);
        ::close(outputFileDescriptor);
        throw runtime_error("Error setting file size for " + outputPath +
            ". Must be a multiple of page size on the target filesystem.");

    }

    // Memory map the input file.
    void* inputPointer = ::mmap(0, n, PROT_READ, MAP_SHARED, inputFileDescriptor, 0);
    ::close(inputFileDescriptor);
    if(inputPointer == reinterpret_cast<void*>(-1LL)) {
        ::close(outputFileDescriptor);
        throw runtime_error("Error mapping " + inputPath + " to memory: " +
            strerror(errno));
    }


    // Memory map the output file.
    void* outputPointer = ::mmap(0, n, PROT_WRITE, MAP_SHARED, outputFileDescriptor, 0);
    ::close(outputFileDescriptor);
    if(outputPointer == reinterpret_cast<void*>(-1LL)) {
        throw runtime_error("Error mapping " + outputPath + " to memory: " +
            strerror(errno));
    }

    // Make the copy.
    const char* inputBegin = static_cast<const char*>(inputPointer);
    const char* inputEnd = inputBegin + n;
    char* outputBegin = static_cast<char*>(outputPointer);
    copy(inputBegin, inputEnd, outputBegin);

    // Unmap.
    if(::munmap(inputPointer, n) == -1) {
        throw runtime_error("Error unmapping " + inputPath + " from memory: " +
            strerror(errno));
    }
    if(::munmap(outputPointer, n) == -1) {
        throw runtime_error("Error unmapping " + outputPath + " from memory: " +
            strerror(errno));
    }

    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    cout << timestamp << "Copied " << n << " bytes in " << tTotal << " s, " << double(n)/tTotal << " bytes/s." << endl;
}
