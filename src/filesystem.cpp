/*******************************************************************************

Replacement for some basic functionality available in boost::filesystem
and std::filesystem.

We don't want to use boost::filesystem to avoid introducing a runtime
dependency on boost libraries.

We don't want to use std::filesystem because it is only available in C++17,
and gcc support for C++17 is still limited (particularly with gcc 4.8
which is the version used in CentOs 7).

*******************************************************************************/

// Nanopore2
#include "filesystem.hpp"
#include "CZI_ASSERT.hpp"
#include "stdexcept.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;
using namespace ChanZuckerberg::Nanopore2::filesystem;

// Linux.
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>


// Return true if the path exists.
bool ChanZuckerberg::Nanopore2::filesystem::exists(const string& path)
{
    struct ::stat info;
    return ::stat(path.c_str(), &info) == 0;
}



// Return true if the path exists and is a regular file.
bool ChanZuckerberg::Nanopore2::filesystem::isRegularFile(const string& path)
{
    struct ::stat info;
    if(::stat(path.c_str(), &info) == -1) {
        return false;
    } else {
        return S_ISREG(info.st_mode);
    }
}



// Return true if the path exists and is a directory.
bool ChanZuckerberg::Nanopore2::filesystem::isDirectory(const string& path)
{
    struct ::stat info;
    if(::stat(path.c_str(), &info) == -1) {
        return false;
    } else {
        return S_ISDIR(info.st_mode);
    }
}



// Create a directory. In case of failure, throw an exception.
void ChanZuckerberg::Nanopore2::filesystem::createDirectory(const string& path)
{
    if(::mkdir(path.c_str(), -1) == -1) {
        throw runtime_error("Unable to create directory " + path);
    }
}


// Remove the specified path. In case of failure, throw an exception.
void ChanZuckerberg::Nanopore2::filesystem::remove(const string& path)
{
    if(::unlink(path.c_str()) == -1) {
        throw runtime_error("Unable to remove directory " + path);
    }
}



// Return the contents of a directory. In case of failure, throw an exception.
vector<string> ChanZuckerberg::Nanopore2::filesystem::directoryContents(const string& path)
{
    DIR* dir = opendir(path.c_str());
    if(!dir) {
        throw runtime_error("Error listing contents of directory " + path);
    }

    vector<string> directoryContents;
    ::dirent* entry = 0;
    while(true) {
        entry = ::readdir(dir);
        if(!entry) {
            break;
        }
        const string name(entry->d_name);
        if(name!="." && name!="..") {
            directoryContents.push_back(path + "/" + name);
        }
    }

    closedir(dir);
    return directoryContents;
}



// Return the extension of a path - that is, everything following
// the last dot after the last slash.
// If there is no dot after the last slash, throw an exception.
string ChanZuckerberg::Nanopore2::filesystem::extension(const string& path)
{
    // If the path is empty, throw an exception.
    if(path.empty()) {
        throw runtime_error("Cannot extract extension of empty path.");
    }

    // Loop backward beginning at the end.
    size_t i = path.size()-1;
    while(true) {
        const char c = path[i];

        // If we find a slash before a dot (looping from the end), there is no extension.
        if(c == '/') {
            throw runtime_error("Cannot extract extension from  path " + path);
        }

        // If we find a dot, return everything that follows it.
        if(c == '.') {
            return path.substr(i+1);
        }

        // If we reached the beginning of the string, there is no extension.
        if(i==0) {
            throw runtime_error("Cannot extract extension from  path " + path);
        }

        // Check the previous character.
        --i;
    }
}



// Return everything up to the last dot following the last dash of a path.
// If there is no dot following the last dash, throw an exception.
string ChanZuckerberg::Nanopore2::filesystem::fileName(const string& path)
{
    // If the path is empty, throw an exception.
    if(path.empty()) {
        throw runtime_error("Cannot extract file name of empty path.");
    }

    // Loop backward beginning at the end.
    size_t i = path.size()-1;
    while(true) {
        const char c = path[i];

        // If we find a slash before a dot (looping from the end), there is no extension.
        if(c == '/') {
            throw runtime_error("Cannot extract file name from  path " + path);
        }

        // If we find a dot, return everything that precedes it.
        if(c == '.') {
            return path.substr(0, i);
        }

        // If we reached the beginning of the string, there is no extension.
        if(i==0) {
            throw runtime_error("Cannot extract file name from  path " + path);
        }

        // Check the previous character.
        --i;
    }
}



// Find the size of a file.
size_t ChanZuckerberg::Nanopore2::filesystem::fileSize(const string& path)
{
    struct ::stat fileInformation;
    if(::stat(path.c_str(), &fileInformation) != 0) {
        throw runtime_error("Could not determine the size of file " + path);
    }
    return fileInformation.st_size;
}
