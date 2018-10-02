#!/usr/bin/python3

import os
import pwd
import sys

helpMessage = """
This script can be used to clean up huge pages after a shasta run 
was completed and saved to disk.
It must run with root privileges.

Invoke without arguments.
"""

if not len(sys.argv) == 1:
    print(helpMessage);
    exit(1)
    
    

largePagesMountPoint = '/hugepages';
Data = '%s/Data' % largePagesMountPoint

response = input('This will destroy any binary data in memory at %s.\nEnter "Yes" if you really want to do that:\n' % Data)   
if not response == 'Yes':
    print('Nothing done.')
    exit(0)


# Create the Data directory.
os.system('sudo rm -rf %s' % Data)

# Unmount the huge page filesystem.
os.system('sudo umount %s' % largePagesMountPoint)

# Remove the mount point for the huge page filesystem.
os.system('sudo rmdir %s' % largePagesMountPoint)

# Free up the huge pages.
os.system('sudo hugeadm --pool-pages-min=2M:0')


