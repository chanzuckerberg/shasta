#ifndef SHASTA_PERFORMANCE_LOG_HPP
#define SHASTA_PERFORMANCE_LOG_HPP

#include "fstream.hpp"
#include "string.hpp"

// The performance log is used to write messages that are useful
// for performance analysis but mostly uninteresting to users.

namespace shasta {
    extern ofstream performanceLog;
    void openPerformanceLog(const string& fileName);
}


#endif
