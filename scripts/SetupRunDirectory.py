#!/usr/bin/python3

import os
import pwd
import sys


def parseArguments():
    helpMessage = """
    This sets up the Data symbolic link required
    to start a shasta run. It also creates the threadLogs directory.
    
    Before running this, use SetuHugePages.py to 
    set up the huge pages.
    
    After this runs, a shasta run can be started in the current directory.
    """

    if not len(sys.argv) == 1:
        print(helpMessage)
        exit(1)


def verifyDirectoryFiles(Data, parentDirectory=""):
    # If any of these is not present, don't do anything.
    mustExist = [Data]
    for name in mustExist:
        if not os.path.lexists(name):
            print('%s must exist. You can use SetupHugePages.py to create it.' % name)
            exit(1)

    # If any of these is present, don't do anything.
    mustNotExist = ['Data', 'threadLogs']
    for name in mustNotExist:
        path = os.path.abspath(os.path.join(parentDirectory, name))

        if os.path.lexists(path):
            print('%s must not exist. Remove it before running this script.' % name)
            exit(1)


def setupRunDirectory(Data, parentDirectory=""):
    # Generate absolute paths to the files that will be created
    localDataPath = os.path.abspath(os.path.join(parentDirectory, "Data"))
    threadLogsPath = os.path.abspath(os.path.join(parentDirectory, "threadLogs"))

    # Create the Data and data symbolic links.
    os.symlink(Data, localDataPath)

    # Create the threadLogs directory.
    os.mkdir(threadLogsPath)


def main():
    largePagesMountPoint = '/hugepages'
    # data = '/dev/shm/data'                            # unused?
    Data = os.path.join(largePagesMountPoint, "Data")

    parseArguments()
    verifyDirectoryFiles(Data=Data)
    setupRunDirectory(Data=Data)


if __name__ == "__main__":
    main()
