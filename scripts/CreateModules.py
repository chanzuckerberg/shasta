#!/usr/bin/python3

import os
import pathlib

"""
This creates C++20 gcc compiled modules (CMI)
for some of the system include files used by Shasta.
This creates:
1. The C++ modules, in directory gcm.cache.
2. File modules.txt which contains a list and locations
   of the compiled modules in gcm.cache. 

To use these modules during a build:
1. Add "-fmodules-ts -fmodule-mapper=modules.txt" to the command line
2. C++ code should use "import <iostream>;" instead of "#include iostream>",
   etc. for other modules.
"""

stdDirectory = 'usr/include/c++/11'
stdIncludes = [
    'algorithm', 
    'array', 
    'chrono', 
    'cstddef', 
    'cstdint', 
    'fstream', 
    'iosfwd', 
    'iostream', 
    'iterator', 
    'map', 
    # 'memory', For some reason this errors out while building.
    'numeric', 
    'set', 
    'span', 
    # 'stdexcept', For some reason this errors out while building.
    # 'string', For some reason this errors out while building.
    'tuple', 
    'utility', 
    'vector',
    ]

modulesFile = open('modules.txt', 'w')
modulesFile.write('$root ' + str(pathlib.Path().absolute()) + '/gcm.cache\n')

for stdInclude in stdIncludes:
    command = 'g++ -std=c++20 -fmodules-ts -c -x c++-system-header ' + stdInclude
    os.system(command)
    modulesFile.write('/%s/%s %s/%s.gcm\n' % (stdDirectory, stdInclude, stdDirectory, stdInclude))
    

