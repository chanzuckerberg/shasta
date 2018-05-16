#!/usr/bin/python3

import os
import pwd
import sys

helpMessage = """
This sets up huge pages for a Nanopore2 run.

After running this, use SetUpRunDirectory.py
to create the required symbolic links and
the threadLogs directory.

This makes the following assumptions:
- /hugepages does not exist.
- /dev/shm/data does not exist.

Invoke with one argument, the number of GB to allocate to large pages.
"""

if not len(sys.argv) == 2:
    print(helpMessage);
    exit(1)
    
gigaBytes = int(sys.argv[1]) 

largePagesMountPoint = '/hugepages';
data = '/dev/shm/data'
Data = '%s/Data' % largePagesMountPoint



# If any of these is present, don't do anything.
mustNotExist = [data, largePagesMountPoint]
for name in mustNotExist:
    if os.path.lexists(name):
        print('%s must not exist. Remove it before running this script.' % name)
        exit(1)
   
# Get the user name.
userName = pwd.getpwuid(os.getuid()).pw_name

# Make sure the hugepages package is installed.
if not os.path.exists('/usr/bin/hugeadm'):
    os.system('sudo apt update')
    os.system('sudo apt upgrade')
    os.system('sudo apt install hugepages')

# Allocate the requested number of pages.
# Consider adding --add-temp-swap option.
os.system('sudo hugeadm --pool-pages-min=2M:%dG' % gigaBytes)

# Create the mount point for the huge page filesystem.
os.system('sudo mkdir %s' % largePagesMountPoint)

# Mount the huge page filesystem.
os.system('sudo mount -t hugetlbfs -o size=100%% none %s' % largePagesMountPoint)

# Create the Data directory.
os.system('sudo mkdir %s' % Data)
os.system('sudo chown %s %s' % (userName, Data))

# Create the data directory.
os.system('mkdir %s' % data)



