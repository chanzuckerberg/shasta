#ifndef SHASTA_HTML_HPP
#define SHASTA_HTML_HPP

#include "iosfwd.hpp"
#include "string.hpp"

// Miscellaneous html related functions.

namespace shasta {

    void writeHtmlBegin(ostream&, const string& title);
    void writeHtmlEnd(ostream&);
    void writeStyle(ostream&);

}

#endif
