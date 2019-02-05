#!/usr/bin/python3

import os
import pwd
import sys


def parseArguments():
    helpMessage = """
    This sets up huge pages for a shasta run.
    It must run with root privileges.

    After running this, use SetUpRunDirectory.py
    to create the required symbolic links and
    the threadLogs directory.

    This makes the following assumptions:
    - /hugepages does not exist.
    - /dev/shm/data does not exist.
    - The hugepages package is installed.

    Invoke with one argument, the number of GB to allocate to large pages.
    """

    if not len(sys.argv) == 2:
        print(helpMessage)
        exit(1)


def verifyPageMemoryDirectory(largePagesMountPoint):
    # If any of these is present, don't do anything.
    mustNotExist = [largePagesMountPoint]
    for name in mustNotExist:
        if os.path.lexists(name):
            print('%s must not exist. Remove it before running this script.' % name)
            exit(1)


def allocatePages(gigaBytes, largePagesMountPoint, Data):

    # Get the user name.
    userName = pwd.getpwuid(os.getuid()).pw_name

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


def main():
    largePagesMountPoint = '/hugepages'
    Data = os.path.join(largePagesMountPoint, "Data")

    gigaBytes = int(sys.argv[1])

    parseArguments()
    verifyPageMemoryDirectory(largePagesMountPoint)
    allocatePages(gigaBytes=gigaBytes, largePagesMountPoint=largePagesMountPoint, Data=Data)


if __name__ == "__main__":
    main()
