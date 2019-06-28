/*******************************************************************************

Replacement for some basic functionality available in boost::filesystem
and std::filesystem.

We don't want to use boost::filesystem to avoid introducing a runtime
dependency on boost libraries.

We don't want to use std::filesystem because it is only available in C++17,
and gcc support for C++17 is still limited (particularly with gcc 4.8
which is the version used in CentOs 7).

*******************************************************************************/

// shasta
#include "filesystem.hpp"
#include "SHASTA_ASSERT.hpp"
#include "stdexcept.hpp"
using namespace ChanZuckerberg;
using namespace ChanZuckerberg::shasta;
using namespace ChanZuckerberg::shasta::filesystem;

// Linux.
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

// Standard library.
#include "array.hpp"


// Return true if the path exists.
bool ChanZuckerberg::shasta::filesystem::exists(const string& path)
{
    struct ::stat info;
    return ::stat(path.c_str(), &info) == 0;
}



// Return true if the path exists and is a regular file.
bool ChanZuckerberg::shasta::filesystem::isRegularFile(const string& path)
{
    struct ::stat info;
    if(::stat(path.c_str(), &info) == -1) {
        return false;
    } else {
        return S_ISREG(info.st_mode);
    }
}



// Return true if the path exists and is a directory.
bool ChanZuckerberg::shasta::filesystem::isDirectory(const string& path)
{
    struct ::stat info;
    if(::stat(path.c_str(), &info) == -1) {
        return false;
    } else {
        return S_ISDIR(info.st_mode);
    }
}



// Create a directory. In case of failure, throw an exception.
void ChanZuckerberg::shasta::filesystem::createDirectory(const string& path)
{

    if(::mkdir(path.c_str(), 0777) == -1) {
        throw runtime_error("Unable to create directory " + path);
    }
}



// Return the current directory.
string ChanZuckerberg::shasta::filesystem::getCurrentDirectory()
{
    const size_t bufferSize = PATH_MAX;
    array<char, bufferSize> buffer;
    ::getcwd(buffer.data(), bufferSize);
    return string(buffer.data());
}



// Change the current directory.
void ChanZuckerberg::shasta::filesystem::changeDirectory(const string& path)
{
    if(::chdir(path.c_str()) == -1) {
        throw runtime_error("Unable to change directory to " + path);
    }

}



// Remove the specified path. In case of failure, throw an exception.
void ChanZuckerberg::shasta::filesystem::remove(const string& path)
{
    if(::unlink(path.c_str()) == -1) {
        throw runtime_error("Unable to remove directory " + path);
    }
}



// Return the contents of a directory. In case of failure, throw an exception.
vector<string> ChanZuckerberg::shasta::filesystem::directoryContents(const string& path)
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
string ChanZuckerberg::shasta::filesystem::extension(const string& path)
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
string ChanZuckerberg::shasta::filesystem::fileName(const string& path)
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
size_t ChanZuckerberg::shasta::filesystem::fileSize(const string& path)
{
    struct ::stat fileInformation;
    if(::stat(path.c_str(), &fileInformation) != 0) {
        throw runtime_error("Could not determine the size of file " + path);
    }
    return fileInformation.st_size;
}



// Find the absolute path.
string ChanZuckerberg::shasta::filesystem::getAbsolutePath(const string& path)
{
    const size_t bufferSize = PATH_MAX;
    array<char, bufferSize> buffer;
    ::realpath(path.c_str(), buffer.data());
    return string(buffer.data());
}
