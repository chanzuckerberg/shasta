#include "SimpleBayesianConsensusCaller.hpp"
using namespace shasta;



// Attempt to construct interpreting the constructor string as
// a built-in configuration name.
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



