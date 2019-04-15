#ifndef CZI_SHASTA_FILESYSTEM_HPP
#define CZI_SHASTA_FILESYSTEM_HPP

/*******************************************************************************

Replacement for some basic functionality available in boost::filesystem
and std::filesystem.

We don't want to use boost::filesystem to avoid introducing a runtime
dependency on boost libraries.

We don't want to use std::filesystem because it is only available in C++17,
and gcc support for C++17 is still limited (particularly with gcc 4.8
which is the version used in CentOs 7).

*******************************************************************************/

#include "string.hpp"
#include "vector.hpp"

namespace ChanZuckerberg {
    namespace shasta {
        namespace filesystem {

            // Return true if the path exists.
            bool exists(const string&);

            // Return true if the path exists and is a regular file.
            bool isRegularFile(const string&);

            // Return true if the path exists and is a directory.
            bool isDirectory(const string&);

            // Create a directory. In case of failure, throw an exception.
            void createDirectory(const string&);

            // Return the current directory.
            string getCurrentDirectory();

            // Change the current directory.
            void changeDirectory(const string&);

            // Remove the specified path. In case of failure, throw an exception.
            void remove(const string&);

            // Return the contents of a directory. In case of failure, throw an exception.
            vector<string> directoryContents(const string&);

            // Return the extension of a path - that is, everything following
            // the last dot after the last slash.
            // If there is no dot after the last slash, throw an exception.
            string extension(const string&);

            // Return everything up to the last dot following the last dash of a path.
            // If there is no dot following the last dash, throw an exception.
            string fileName(const string&);

            // Find the size of a file.
            size_t fileSize(const string&);

            // Find the absolute path.
            string getAbsolutePath(const string& path);
        }
    }
}



#endif
