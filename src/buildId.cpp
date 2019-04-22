#include "buildId.hpp"

// Some contortions forced by the preprocessor.
#define STRINGIZE2(s) #s
#define STRINGIZE(s) STRINGIZE2(s)

const std::string ChanZuckerberg::shasta::buildId = STRINGIZE(BUILD_ID);
