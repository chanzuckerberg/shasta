#!/usr/bin/python3

import os
import pwd
import sys

helpMessage = """
This script can be used to set up the run directory for a small
run for which performance is not important.

All binary files are created on disk, not on the huge page filesystem.

This script creates empty directories data, data, and threadLogs directory.
It also creates symbolik links dataOnDist and DataOnDisk pointing
to data and Data respectively.

When using this script, there is not need to use SetuHugePages.py to 
set up the huge pages.

After this runs, a shasta run can be started in the current directory.
"""

if not len(sys.argv) == 1:
    print(helpMessage);
    exit(1)
    

# If any of these is present, don't do anything.
mustNotExist = ['data', 'Data', 'threadLogs']
for name in mustNotExist:
    if os.path.lexists(name):
        print('%s must not exist. Remove it before running this script.' % name)
        exit(1)

# Create the directories.
os.mkdir('data')
os.mkdir('Data')
os.mkdir('threadLogs')

# Create the symbolic links.
os.symlink('data', 'dataOnDisk')
os.symlink('Data', 'DataOnDisk')



