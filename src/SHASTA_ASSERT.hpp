// Definition of macro SHASTA_ASSERT.
// It is always compiled in, regardless of compilation settings.
// It throws a standard exception if the assertion fails.

#ifndef SHASTA_SHASTA_ASSERT_HPP
#define SHASTA_SHASTA_ASSERT_HPP

#include <stdexcept>
#include <string>

#if 0
// Gcc (for backtraces).
#include "execinfo.h"

namespace shasta {
    inline void writeBackTrace();
}
#endif


#define SHASTA_ASSERT(expression) ((expression) ? (static_cast<void>(0)) : \
    (/*writeBackTrace(),*/ throw std::runtime_error(std::string("Assertion failed: ") + #expression + " at " + __PRETTY_FUNCTION__ + " in " +  __FILE__ + " line " + std::to_string(__LINE__))))


#if 0
inline void shasta::writeBackTrace()
{
    const int bufferSize = 64;  // To avoid extremely long, useless backtraces.
    void* buffer[bufferSize];
    ::backtrace(buffer, bufferSize);
    ::backtrace_symbols_fd(buffer, bufferSize, ::fileno(::stdout));
}
#endif


#endif

