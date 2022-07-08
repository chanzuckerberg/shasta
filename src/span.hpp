#ifndef SHASTA_SPAN_HPP
#define SHASTA_SPAN_HPP

#include "algorithm.hpp"
#include "iostream.hpp"
#include "iterator.hpp"
#include <span>
#include "stdexcept.hpp"
#include "string.hpp"

namespace shasta {
    using std::span;

    // Output a span<const char> as a string.
    inline ostream& operator<<(ostream&, const span<const char>&);

    // Convert a span<const char> to an integer.
    // This cannot be done using std::atol because the
    // span<const char> is not null terminated.
    uint64_t atoul(const span<const char>&);

    // Convert a span<const char> to a string.
    string convertToString(const span<const char>&);
    string convertToString(const span<char>&);

    // Comparison operators, which are missing from std::span.
    template<class T> inline bool operator==(const span<T>& x, const span<T>& y)
    {
        if(x.size() == y.size()) {
            return std::equal(x.begin(), x.end(), y.begin());
        } else {
            return false;
        }
    }
    template<class T> inline bool operator!=(const span<T>& x, const span<T>& y)
    {
        return not(x == y);
    }
    template<class T> inline bool operator<(const span<T>& x, const span<T>& y)
    {
        return std::lexicographical_compare(
            x.begin(), x.end(),
            y.begin(), y.end());
    }
}



// Output a span<const char> as a string.
inline std::ostream& shasta::operator<<(
    std::ostream& s,
    const std::span<const char>&  m)
{
    copy(m.begin(), m.end(), std::ostream_iterator<char>(s));
    return s;
}



// Convert a span<const char> to an integer.
// This cannot be done using std::atol because the
// span<const char> is not null terminated.
// The string can only contain numeric characters.
inline uint64_t shasta::atoul(const std::span<const char>& s)
{
    uint64_t n = 0;
    for(uint64_t i=0; ; i++) {
        const char c = s[i];
        if(not std::isdigit(c)) {
            throw runtime_error("Non-digit found in " + convertToString(s));
        }

        n += (c- '0');

        if(i ==  s.size()-1) {
            return n;
        }

        n *= 10;
    }
}



// Convert a span<const char> to an std::string.
inline std::string shasta::convertToString(const std::span<const char>& m)
{
    auto view = std::string_view(m.data(), m.size());
    return string(view);
}
inline std::string shasta::convertToString(const std::span<char>& m)
{
    auto view = std::string_view(m.data(), m.size());
    return string(view);
}



#endif
