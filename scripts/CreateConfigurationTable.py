#!/usr/bin/python3

import os

"""
This must be run from the shasta directory.
It expects to see a conf and src directory.

This creates shasta/src/ConfigurationTable.cpp which contains
the definition of shasta::configurationTable,
a map that map configuration names to strings.

The file will be overwritten if it exists.

Each entry inthe table corresponds to a file in shasta/conf.

"""

# Check that we have the conf and src directories.
if not (os.path.isdir('conf') and os.path.isdir('src')):
    raise Exception('Must run from the shasta directory')



# The list of configurations that will included in the table.
# Each of them must have a corresponding file in
# the conf directory, with the same name plus a .conf extension.
# Add them in approximately chronological order.
configurations = [
    'Nanopore-Dec2019',
    'Nanopore-UL-Dec2019',
    'Nanopore-Jun2020',
    'Nanopore-UL-Jun2020',
    'Nanopore-UL-iterative-Sep2020',
    'Nanopore-Sep2020',
    'Nanopore-UL-Sep2020',
    'Nanopore-OldGuppy-Sep2020',
    'Nanopore-Plants-Apr2021',
    'Nanopore-Oct2021',
    'Nanopore-UL-Oct2021',
    'HiFi-Oct2021',
    ]



# Before doing anything, check that we have all the files.    
for configuration in configurations:    
    fileName = 'conf/' + configuration + '.conf'
    if not os.path.isfile(fileName):
        raise Exception(fileName + ' not found.')
        
# Open the output file.
out = open('src/ConfigurationTable.cpp', 'w')   

# Write the initial portion of the file.
out.write("""
#include "ConfigurationTable.hpp"

namespace shasta {
   const std::map<string, string> configurationTable = {
""")

# Write the configurations.
for i in range(len(configurations)):
    configuration = configurations[i]    
    fileName = 'conf/' + configuration + '.conf'
    configFile = open(fileName, 'r')
    out.write('    {"' + configuration + '", R"zzz(')
    out.write(configFile.read())
    out.write(')zzz"}');
    if not i == len(configurations) - 1:
        out.write(',')    
    out.write('\n');

# Write the final portion of the file.
out.write('    };\n}\n')



    



