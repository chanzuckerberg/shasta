#ifndef SHASTA_PLATFORM_DEPENDENT_HPP
#define SHASTA_PLATFORM_DEPENDENT_HPP

#include "string.hpp"

namespace shasta {
    
    // Return the path to a usable temporary directory, including the final "/".
    string tmpDirectory();
    
    // Return the name of a timeout command or equivalent.
    string timeoutCommand();

    // Get peak virtual memory utilization of the current process, in bytes.
    uint64_t getPeakMemoryUsage();

    // Get total physical memory available, in bytes.
    uint64_t getTotalPhysicalMemory();
}

#endif
