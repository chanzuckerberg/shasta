#ifndef SHASTA_FILESYSTEM_HPP
#define SHASTA_FILESYSTEM_HPP

#include "string.hpp"
#include "vector.hpp"

namespace shasta {
    namespace filesystem {

        // Return the contents of a directory. In case of failure, throw an exception.
        vector<string> directoryContents(const string&);

        // Return the extension of a path - that is, everything following
        // the last dot after the last slash.
        // If there is no dot after the last slash, throw an exception.
        string extension(const string&);

        // Find the absolute path.
        string getAbsolutePath(const string& path);

        // Find the absolute path of the executable. Works only for Linux.
        string executablePath();

    }
}



#endif
