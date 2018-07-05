// Definition of macro CZI_ASSERT.
// It is always compiled in, regardless of compilation settings.
// It throws a standard exception if the assertion fails.

#ifndef CZI_SHASTA_CZI_ASSERT_HPP
#define CZI_SHASTA_CZI_ASSERT_HPP

#include <stdexcept>
#include <string>

// Gcc (for backtraces).
#include "execinfo.h"

namespace ChanZuckerberg {
    namespace Nanopore1 {
        inline void writeBackTrace();
    }
}

#define CZI_ASSERT(expression) ((expression) ? (static_cast<void>(0)) : \
    (/*writeBackTrace(),*/ throw std::runtime_error(std::string("Assertion failed: ") + #expression + " at " + __PRETTY_FUNCTION__ + " in " +  __FILE__ + " line " + std::to_string(__LINE__))))


#if 0
inline void ChanZuckerberg::Nanopore2::writeBackTrace()
{
    const int bufferSize = 64;  // To avoid extremely long, useless backtraces.
    void* buffer[bufferSize];
    ::backtrace(buffer, bufferSize);
    ::backtrace_symbols_fd(buffer, bufferSize, ::fileno(::stdout));
}
#endif


#endif

