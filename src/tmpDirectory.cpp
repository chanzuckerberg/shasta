#include "tmpDirectory.hpp"
#ifdef __linux__
#include <stdlib.h>
#endif

// Return the path to a usable temporary directory, including the final "/".
std::string shasta::tmpDirectory()
{
#ifdef __linux__

    // For Linux, use /dev/shm so the temporary files never go to disk.
    return "/dev/shm/";
    
#else

    // For macOS, use environment variable TMPDIR.
    return string(::getenv("TMPDIR")) + "/";
    
#endif
}
