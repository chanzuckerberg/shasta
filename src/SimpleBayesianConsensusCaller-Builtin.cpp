#include "SimpleBayesianConsensusCaller.hpp"
using namespace shasta;

// Note: This needs to be kept in sync with constructBuiltin()
bool SimpleBayesianConsensusCaller::isBuiltIn(const string& constructorString) {
    return constructorString == "guppy-2.3.1-a" or
        constructorString == "guppy-3.0.5-a" or
        constructorString == "guppy-3.4.4-a" or
        constructorString == "guppy-3.6.0-a" or
        constructorString == "r10-guppy-3.4.8-a" or
        constructorString == "bonito-0.3.1-a" or
        constructorString == "guppy-5.0.7-a";
}

// Attempt to construct interpreting the constructor string as
// a built-in configuration name.
// Note: This needs to be kept in sync with isBuiltIn()
bool SimpleBayesianConsensusCaller::constructBuiltin(const string& constructorString)
{
    if(constructorString == "guppy-2.3.1-a"){
        // From SimpleBayesianConsensusCaller-3.csv
        #include "SimpleBayesianConsensusCaller-Builtin-guppy-2.3.1-a.hpp"
        return true;
    }

    if(constructorString == "guppy-3.0.5-a"){
        // From SimpleBayesianConsensusCaller-5.csv
        #include "SimpleBayesianConsensusCaller-Builtin-guppy-3.0.5-a.hpp"
        return true;
    }

    if(constructorString == "guppy-3.4.4-a"){
        // From SimpleBayesianConsensusCaller-6.csv
        #include "SimpleBayesianConsensusCaller-Builtin-guppy-3.4.4-a.hpp"
        return true;
    }

    if(constructorString == "guppy-3.6.0-a"){
        // From SimpleBayesianConsensusCaller-7.csv
        #include "SimpleBayesianConsensusCaller-Builtin-guppy-3.6.0-a.hpp"
        return true;
    }

    if(constructorString == "r10-guppy-3.4.8-a"){
        // From SimpleBayesianConsensusCaller-8.csv
        #include "SimpleBayesianConsensusCaller-Builtin-r10-guppy-3.4.8-a.hpp"
        return true;
    }

    if(constructorString == "bonito-0.3.1-a"){
        // From SimpleBayesianConsensusCaller-9.csv
        #include "SimpleBayesianConsensusCaller-Builtin-bonito-0.3.1-a.hpp"
        return true;
    }

    if(constructorString == "guppy-5.0.7-a"){
        // From SimpleBayesianConsensusCaller-10.csv
        #include "SimpleBayesianConsensusCaller-Builtin-guppy-5.0.7-a.hpp"
        return true;
    }

    return false;
}



