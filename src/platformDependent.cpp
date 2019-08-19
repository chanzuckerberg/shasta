#include "platformDependent.hpp"
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



// Return the name of a timeout command or equivalent.
std::string shasta::timeoutCommand()
{
#ifdef __linux__
    return "timeout";
#else

    // For macOS, return gtimeout.
    // Requires brew install coreutils.
    return "gtimeout";
    
#endif
    
}
