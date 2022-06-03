#ifndef SHASTA_FILESYSTEM_HPP
#define SHASTA_FILESYSTEM_HPP

#include "string.hpp"
#include "vector.hpp"

namespace shasta {
    namespace filesystem {

        // Move (rename). In case of failure, throw an exception.
        void move(const string& oldPath, const string& newPath);

        // Copy a file.
        void copy(const string&, const string&);

        // Return the contents of a directory. In case of failure, throw an exception.
        vector<string> directoryContents(const string&);

        // Return the extension of a path - that is, everything following
        // the last dot after the last slash.
        // If there is no dot after the last slash, throw an exception.
        string extension(const string&);

        // Return everything up to the last dot following the last dash of a path.
        // If there is no dot following the last dash, throw an exception.
        string fileName(const string&);

        // Find the absolute path.
        string getAbsolutePath(const string& path);

        // Find the absolute path of the executable. Works only for Linux.
        string executablePath();

    }
}



#endif
