#include "platformDependent.hpp"
#include <stdlib.h>
#include "fstream.hpp"

// Return the path to a usable temporary directory, including the final "/".
std::string shasta::tmpDirectory()
{
    // Use /dev/shm, so the temporary files never go to disk.
    return "/dev/shm/";
}



// Return the name of a timeout command or equivalent.
std::string shasta::timeoutCommand()
{
    return "timeout";
}



uint64_t shasta::getPeakMemoryUsage() {
    uint64_t peakMemoryUsage = 0ULL;

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

    return peakMemoryUsage;
}



// Get total physical memory available, in bytes.
uint64_t shasta::getTotalPhysicalMemory()
{
    ifstream meminfo("/proc/meminfo");
    string s;
    uint64_t memoryKb;
    meminfo >> s >> memoryKb;
    return 1024 * memoryKb;
}

