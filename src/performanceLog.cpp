// The performance log is used to write messages that are useful
// for performance analysis but mostly uninteresting to users.

#include "performanceLog.hpp"

namespace shasta {
    ofstream performanceLog;
}



void shasta::openPerformanceLog(const string& fileName)
{
    performanceLog.open(fileName);
}
