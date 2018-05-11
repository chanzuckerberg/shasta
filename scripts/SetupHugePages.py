#!/usr/bin/python3

import os
import pwd
import sys

helpMessage = """
This sets up huge pages for a Nanopore2 run
and creates data and Data symbolic links in the current directory.
It also creates a theadLogs directory.

This makes the following assumptions, which are true for a fresh AWS EC2 
instance:
- /hugepages does not exist.
- /dev/shm/data does not exist.

After this runs, a Nanopore2 run can be started in the current directory.

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
mustNotExist = [data, largePagesMountPoint, 'data', 'Data', 'threadLogs']
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
os.system('sudo hugeadm --pool-pages-min=2M:%dG --add-temp-swap' % gigaBytes)

# Create the mount point for the huge page filesystem.
os.system('sudo mkdir %s' % largePagesMountPoint)

# Mount the huge page filesystem.
os.system('sudo mount -t hugetlbfs -o size=100%% none %s' % largePagesMountPoint)

# Create the Data directory and the symbolic link to it.
os.system('sudo mkdir %s' % Data)
os.system('sudo chown %s %s' % (userName, Data))
os.symlink(Data, 'Data')

# Create the data directory and the symbolic link to it.
os.system('mkdir %s' % data)
os.symlink(data, 'data')

# Crreate the threadLogs directory.
os.system('mkdir threadLogs')


