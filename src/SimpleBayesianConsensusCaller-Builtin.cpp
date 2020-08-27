#include "SimpleBayesianConsensusCaller.hpp"
using namespace shasta;

// Note: This needs to be kept in sync with constructBuiltin()
bool SimpleBayesianConsensusCaller::isBuiltIn(const string& constructorString) {
    return constructorString == "guppy-2.3.1-a" or
        constructorString == "guppy-2.3.5-a" or
        constructorString == "guppy-3.0.5-a" or
        constructorString == "guppy-3.4.4-a" or
        constructorString == "guppy-3.6.0-a" or
        constructorString == "r10-guppy-3.4.8-a";
}

// Attempt to construct interpreting the constructor string as
// a built-in configuration name.
// Note: This needs to be kept in sync with isBuiltIn()
bool SimpleBayesianConsensusCaller::constructBuiltin(const string& constructorString)
{
    if(constructorString == "guppy-2.3.1-a"){
        #include "SimpleBayesianConsensusCaller-Builtin-guppy-2.3.1-a.hpp"
        return true;
    }

    if(constructorString == "guppy-2.3.5-a"){
        #include "SimpleBayesianConsensusCaller-Builtin-guppy-2.3.5-a.hpp"
        return true;
    }

    if(constructorString == "guppy-3.0.5-a"){
        #include "SimpleBayesianConsensusCaller-Builtin-guppy-3.0.5-a.hpp"
        return true;
    }

    if(constructorString == "guppy-3.4.4-a"){
        #include "SimpleBayesianConsensusCaller-Builtin-guppy-3.4.4-a.hpp"
        return true;
    }

    if(constructorString == "guppy-3.6.0-a"){
        #include "SimpleBayesianConsensusCaller-Builtin-guppy-3.6.0-a.hpp"
        return true;
    }

    if(constructorString == "r10-guppy-3.4.8-a"){
        #include "SimpleBayesianConsensusCaller-Builtin-r10-guppy-3.4.8-a.hpp"
        return true;
    }
    return false;
}



