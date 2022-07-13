#ifndef SHASTA_INVALID_HPP
#define SHASTA_INVALID_HPP

// In many contexts, we use invalid<uint64_t> (or similar for other integer types)
// to indicate a value that is invalid, uninitialized, or unknown.

#include <concepts>
#include <numeric>

namespace shasta {
    template<std::integral Int> static const Int invalid = std::numeric_limits<Int>::max();
}

#endif
