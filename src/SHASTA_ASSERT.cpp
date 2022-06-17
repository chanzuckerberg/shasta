#include "SHASTA_ASSERT.hpp"

#include <stdexcept>
#include <string>

#if 0
#include "execinfo.h"
static void shasta::writeStackTrace()
{
    const int bufferSize = 64;  // To avoid extremely long, useless backtraces.
    void* buffer[bufferSize];
    ::backtrace(buffer, bufferSize);
    ::backtrace_symbols_fd(buffer, bufferSize, ::fileno(::stdout));
}
#endif


void shasta::handleFailedAssertion(
    const char* expression,
    const char* function,
    const char* file,
    int line)
{
    throw std::runtime_error(
        std::string("Assertion failed: ") + expression +
        " at " + function +
        " in " +  file +
        " line " + std::to_string(line));
}
