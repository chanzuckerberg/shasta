#ifndef SHASTA_CONFIGURATION_TABLE_HPP
#define SHASTA_CONFIGURATION_TABLE_HPP

#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {

    // Make it a vector and not a map, so --command listConfigurations
    // lists configurations in the same order as they are entered in the
    // configuration table (as establisted by scripts/CreateConfigurationTable.py).
    extern const std::vector< pair<string, string> > configurationTable;

    // Return a pointer to the configuration string with the
    // given name, or 0 if not found.
    inline const string* getConfiguration(const string& name)
    {
        for(const auto& p:configurationTable) {
            if(p.first == name) {
                return &p.second;
            }
        }

        // We did not find it.
        return 0;
    }
}

#endif
