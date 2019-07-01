#include "buildId.hpp"

// Some contortions forced by the preprocessor.
#define STRINGIZE2(s) #s
#define STRINGIZE(s) STRINGIZE2(s)

std::string shasta::buildId()
{
    return STRINGIZE(BUILD_ID);
}
