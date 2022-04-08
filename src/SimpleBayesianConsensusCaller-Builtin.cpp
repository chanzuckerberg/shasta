#include "SimpleBayesianConsensusCaller.hpp"
using namespace shasta;


// The list of built-ins Bayesian model.
// This must be kept in sync with constructBuiltin below.
const std::set<string> shasta::SimpleBayesianConsensusCaller::builtIns =
{
    "guppy-2.3.1-a",
    "guppy-3.0.5-a",
    "guppy-3.4.4-a",
    "guppy-3.6.0-a",
    "r10-guppy-3.4.8-a",
    "bonito-0.3.1-a",
    "guppy-5.0.7-a",
    "guppy-5.0.7-b"
};



bool SimpleBayesianConsensusCaller::isBuiltIn(const string& constructorString) {
    return builtIns.find(constructorString) != builtIns.end();
}



// Attempt to construct interpreting the constructor string as
// a built-in configuration name.
// Note: This needs to be kept in sync with isBuiltIn()
bool SimpleBayesianConsensusCaller::constructBuiltin(const string& constructorString)
{
    if(not isBuiltIn(constructorString)) {
        return false;
    }

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

    if(constructorString == "guppy-5.0.7-b"){
        // From SimpleBayesianConsensusCaller-11.csv
        #include "SimpleBayesianConsensusCaller-Builtin-guppy-5.0.7-b.hpp"
        return true;
    }

    // We already checked for a valid built-in before,
    // so if getting here there is missing logic above.
    SHASTA_ASSERT(0);
}



