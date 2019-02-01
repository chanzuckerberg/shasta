#!/usr/bin/python3

import os
import pwd
import sys


def parseArguments():
    helpMessage = """
    This script can be used to clean up huge pages after a shasta run 
    was completed and saved to disk.
    It must run with root privileges.
    
    Invoke without arguments.
    """
    
    if not len(sys.argv) == 1:
        print(helpMessage)
        exit(1)


def cleanUpHugePages(Data, largePagesMountPoint, requireUserInput=True):
    if requireUserInput:
        response = input('This will destroy any binary data in memory at %s.\nEnter (Yes/yes/Y/y) if you really want to do that:\n' % Data)   
        if not response in set(['Yes', 'Y', 'y', 'yes']):
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


def main():
    largePagesMountPoint = '/hugepages'
    Data = os.path.join(largePagesMountPoint, "Data")
    
    parseArguments()
    cleanUpHugePages(Data=Data, largePagesMountPoint=largePagesMountPoint)


if __name__ == "__main__":
    main()
