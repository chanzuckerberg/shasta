#!/usr/bin/python3

import os
import sys


def parseArguments():
    helpMessage = """
    This sets up the Data directory required for a Shasta run.
    The Data directory becomes the mount point of a huge page filesystem.
    
    After this runs, you can use RunAssembly.py to start a Shasta run can be started in the current directory.
    Invoke without arguments.
    
    This script uses sudo to invoke commands that requires root privileges,
    so depending on your sudo set up it may ask for your password.

    """

    if not len(sys.argv) == 1:
        print(helpMessage)
        exit(1)



# This is used to check that we don't overwrite existing data.
def verifyDirectoryFiles(runDirectory = ''):

    # If any of these is present in the run directory, don't do anything.
    mustNotExist = ['Data']
    for name in mustNotExist:
        path = os.path.abspath(os.path.join(runDirectory, name))

        if os.path.lexists(path):
            print('%s must not exist. Remove it before running this script.' % name)
            exit(1)



def setupRunDirectory(runDirectory = ''):

    # Define the minimum and maximum amount of huge page memory.
    # Leaving the minimum at zero has pro's and con's:
    # Good: we don't need to worry about unmounting the
    MB = 1024 * 1024
    GB = 1024 * MB
    totalSystemMemoryBytes =  os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
    minimumHugePageMemoryBytes = 0   # Consider increasing this
    maximumHugePageMemoryBytes = totalSystemMemoryBytes - 8 * GB
    maximumHugePageMemoryBytes = max(maximumHugePageMemoryBytes, minimumHugePageMemoryBytes)
    hugePageSize = 2 * MB
    minimumHugePageMemoryHugePages = int(minimumHugePageMemoryBytes / hugePageSize)
    maximumHugePageMemoryHugePages = int(maximumHugePageMemoryBytes / hugePageSize)
    os.system('sudo sh -c "echo %s > /sys/kernel/mm/hugepages/hugepages-2048kB/nr_hugepages"' 
        % minimumHugePageMemoryHugePages)
    os.system('sudo sh -c "echo %s > /sys/kernel/mm/hugepages/hugepages-2048kB/nr_overcommit_hugepages"' 
        % maximumHugePageMemoryHugePages)
    

    # Create the Data directory.
    dataPath = os.path.abspath(os.path.join(runDirectory, 'Data'))
    os.mkdir(dataPath)
    
    # Mount the huge page filesystem.
    command = 'sudo mount -t hugetlbfs -o uid=%s,gid=%s,pagesize=2M none %s' % (os.getuid(), os.getgid(), dataPath)
    os.system(command)



def main():

    parseArguments()
    
    # Check that we don't overwrite existing data.
    verifyDirectoryFiles()
    
    # Create the Data and threadlogs directories.
    setupRunDirectory()


if __name__ == "__main__":
    main()
