// Main program for the Shasta static executable.
// The static executable provides
// basic functionality and reduced performance. 
// For full functionality use the shared library built
// under directory src.

// Shasta.
#include "Assembler.hpp"
namespace ChanZuckerberg {
    namespace shasta {
        void shastaMain(int argumentCount, const char** arguments);
        class AssemblyOptions;
    }
}
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/program_options.hpp>

// Standard library.
#include "fstream.hpp"
#include "iostream.hpp"
#include "stdexcept.hpp"



// In class AssemblyOptions, we choose names that are
// consistent with options names in shasta.conf.
class ChanZuckerberg::shasta:: AssemblyOptions {
public:
    class ReadsOptions {
    public:
        // useRunLengthReads; This will be phased out and always be true.
        int minReadLength;
        class PalindromicReadOptions {
        public:
            int maxSkip;
            int maxMarkerFrequency;
            double alignedFractionThreshold;
            double nearDiagonalFractionThreshold;
            int deltaThreshold;
        };
        PalindromicReadOptions palindromicReads;
    };
    ReadsOptions Reads;
};



void ChanZuckerberg::shasta::shastaMain(int argumentCount, const char** arguments)
{
    // Some names in the boost program_options library.
    using boost::program_options::command_line_parser;
    using boost::program_options::options_description;
    using boost::program_options::value;
    using boost::program_options::variables_map;
    
    // Options that are only allowed only on the command line.
    options_description commandLineOnlyOptions(
        "Options allowed only on the command line");
    string configFileName;
    commandLineOnlyOptions.add_options()
        ("help", 
        "Write a help message.")
        ("config", 
        boost::program_options::value<string>(&configFileName),
        "Configuration file name.")
        ;
 
 
        
    // Options allowed on the command line and in the config file.
    // Values specified in the command line take precedence.
    options_description options(
        "Options allowed on the command line and in the config file");
    AssemblyOptions assemblyOptions;
    vector<string> inputFastaFileNames;
    
    options.add_options()
        ("input", 
        value< vector<string> >(&inputFastaFileNames), 
        "Names of input FASTA files. Specify at least one.")
        
        ("Reads.minReadLength", 
        value<int>(&assemblyOptions.Reads.minReadLength)->default_value(10000), 
        "Read length cutoff.")
    
        ("Reads.palindromicReads.maxSkip", 
        value<int>(&assemblyOptions.Reads.palindromicReads.maxSkip)->default_value(100), 
        "Used for palindromic read detection.")
    
        ("Reads.palindromicReads.maxMarkerFrequency", 
        value<int>(&assemblyOptions.Reads.palindromicReads.maxMarkerFrequency)->default_value(10), 
        "Used for palindromic read detection.")
    
        ("Reads.palindromicReads.alignedFractionThreshold", 
        value<double>(&assemblyOptions.Reads.palindromicReads.alignedFractionThreshold)->default_value(0.1), 
        "Used for palindromic read detection.")
    
        ("Reads.palindromicReads.nearDiagonalFractionThreshold", 
        value<double>(&assemblyOptions.Reads.palindromicReads.nearDiagonalFractionThreshold)->default_value(0.10), 
        "Used for palindromic read detection.")
    
       ("Reads.palindromicReads.deltaThreshold", 
        value<int>(&assemblyOptions.Reads.palindromicReads.deltaThreshold)->default_value(100), 
        "Used for palindromic read detection.")
        ;
    

        
    // Get options from the command line.
    // These take precedence over values entered in the config file.
    options_description commandLineOptions;
    commandLineOptions.add(commandLineOnlyOptions);
    commandLineOptions.add(options);
    variables_map variablesMap;
    store(command_line_parser(argumentCount, arguments).
          options(commandLineOptions).run(), variablesMap);
    notify(variablesMap);
    
    if (variablesMap.count("help")) {
        cout << commandLineOptions << endl;
        return;
    }      
      
    // Get options from the config file, if one was specified.
    if(!configFileName.empty()) {
        ifstream configFile(configFileName);
        if (!configFile) {
            throw runtime_error("Unable to open open config file " + configFileName);
        }
        store(parse_config_file(configFile, options), variablesMap);
        notify(variablesMap);
    }
       
    // Check that we have at least one input FASTA file.     
    if (inputFastaFileNames.empty()) {
        cout << commandLineOptions << endl;
        throw runtime_error("Specify at least one input FASTA file.");
    }      
    
    // The rest is not implemented.
    throw runtime_error("The Shasta static executable is not yet functional.");
    // Assembler assembler("Data/", 2*1024*1024, true);
}



int main(int argumentCount, const char** arguments)
{
    try {
    
        shastaMain(argumentCount, arguments);      
          
    } catch(boost::program_options::error_with_option_name e) {
        cout << "Invalid option: " << e.what() << endl;
        return 1;
    } catch (runtime_error e) {
        cout << "Terminated after catching a runtime error exception:" << endl;
        cout << e.what() << endl;
        return 2;
    } catch (exception e) {
        cout << "Terminated after catching a standard exception:" << endl;
        cout << e.what() << endl;
        return 3;
    } catch (...) {
        cout << "Terminated after catching a non-standard exception." << endl;
        return 4;
    }
    return 0;
}
