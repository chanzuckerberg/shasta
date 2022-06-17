// Definition of macro SHASTA_ASSERT.
// It is always compiled in, regardless of compilation settings.
// It throws a standard exception if the assertion fails.

#ifndef SHASTA_SHASTA_ASSERT_HPP
#define SHASTA_SHASTA_ASSERT_HPP

namespace shasta {
    void handleFailedAssertion(
        const char* expression,
        const char* function,
        const char* file,
        int line
    ) __attribute__ ((__noreturn__));
}


#define SHASTA_ASSERT(expression) ((expression) ? (static_cast<void>(0)) : \
    (shasta::handleFailedAssertion(#expression, __PRETTY_FUNCTION__,  __FILE__ , __LINE__)))


#endif

