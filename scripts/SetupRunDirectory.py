#!/usr/bin/python3

import os
import pwd
import sys

helpMessage = """
This sets up the symbolic links (data and Data) required
to start a shasta run. It also creates the threadLogs directory.

Before running this, use SetuHugePages.py to 
set up the huge pages.

After this runs, a shasta run can be started in the current directory.
"""

if not len(sys.argv) == 1:
    print(helpMessage);
    exit(1)
    
largePagesMountPoint = '/hugepages';
data = '/dev/shm/data'
Data = '%s/Data' % largePagesMountPoint


# If any of these is not present, don't do anything.
mustExist = [data, Data]
for name in mustExist:
    if not os.path.lexists(name):
        print('%s must exist. You can use SetupHugePages.py to create it.' % name)
        exit(1)

# If any of these is present, don't do anything.
mustNotExist = ['data', 'Data', 'threadLogs']
for name in mustNotExist:
    if os.path.lexists(name):
        print('%s must not exist. Remove it before running this script.' % name)
        exit(1)

# Create the Data and data symbolic links.
os.symlink(Data, 'Data')
os.symlink(data, 'data')

# Crreate the threadLogs directory.
os.system('mkdir threadLogs')


