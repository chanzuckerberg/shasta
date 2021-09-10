#ifndef SHASTA_CONFIGURATION_TABLE_HPP
#define SHASTA_CONFIGURATION_TABLE_HPP

#include <map>
#include "string.hpp"

namespace shasta {
    extern const std::map<string, string> configurationTable;

    // Return a pointer to the configuration string with the
    // given name, or 0 if not found.
    inline const string* getConfiguration(const string& name)
    {
        const auto it = configurationTable.find(name);
        if(it == configurationTable.end()) {
            return 0;
        } else {
            return &(it->second);
        }
    }
}

#endif
