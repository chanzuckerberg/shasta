#ifndef SHASTA_TMP_DIRECTRORY_HPP
#define SHASTA_TMP_DIRECTRORY_HPP

#include "string.hpp"

namespace shasta {
    
    // Return the path to a usable temporary directory, including the final "/".
    string tmpDirectory();
}

#endif
