#!/usr/bin/python3

"""
This checks that all built-in configurations work with
the current shasta executable.
It must run from the shasta-install/bin directory.
"""

import shasta
import os

badCount = 0;
for p in shasta.configurationTable:
    configurationName = p[0]
    command = './shasta --command listConfiguration --config ' + configurationName
    returnCode = os.system(command + ' > /dev/null')
    if not (returnCode == 0):
    	print(configurationName, 'does not work.')
    	badCount = badCount + 1

if badCount== 0:
    print('All built-in configurations work.')
else:
    print(badCount, 'built-in configurations are not working.')



