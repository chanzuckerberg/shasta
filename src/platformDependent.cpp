#include "platformDependent.hpp"
#ifdef __linux__
#include <stdlib.h>
#include "fstream.hpp"
#endif

// Return the path to a usable temporary directory, including the final "/".
std::string shasta::tmpDirectory()
{
#ifdef __linux__

    // For Linux, use /dev/shm so the temporary files never go to disk.
    return "/dev/shm/";
    
#else

    // For macOS, use environment variable TMPDIR.
    return string(::getenv("TMPDIR")) + "/";
    
#endif
}



// Return the name of a timeout command or equivalent.
std::string shasta::timeoutCommand()
{
#ifdef __linux__
    return "timeout";
#else

    // For macOS, return gtimeout.
    // Requires brew install coreutils.
    return "gtimeout";
    
#endif
    
}

uint64_t shasta::getPeakMemoryUsage() {
    uint64_t peakMemoryUsage = 0ULL;
#ifdef __linux__
    ifstream procStats("/proc/self/status");
    if (procStats) {
        string line;
        while (std::getline(procStats, line)) {
            if (string::npos == line.find("VmPeak")) {
                continue;
            }
            size_t pos = line.find(":");
            while (pos < line.size() && !isdigit(line[pos])) {
                pos++;
            }
            char* end;
            peakMemoryUsage = std::strtoull(line.c_str() + pos, &end, 10);
            // Convert from kB to bytes.
            peakMemoryUsage *= 1024;
            break;
        }
    }
#endif
    return(peakMemoryUsage);
}
