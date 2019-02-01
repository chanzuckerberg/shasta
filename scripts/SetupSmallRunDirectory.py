#!/usr/bin/python3

import os
import pwd
import sys


def verifyArguments():
    """
    Parse Fasta path arguments
    """
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
        print(helpMessage)
        exit(1)


def verifyDirectoryFiles(parentDirectory=""):
    """
    Make sure the run directory is clean before starting
    """
    # If any of these is present, don't do anything.
    mustNotExist = ['Data', 'threadLogs']
    
    for name in mustNotExist:
        path = os.path.abspath(os.path.join(parentDirectory, name))
        
        if os.path.lexists(path):
            print('%s must not exist. Remove it before running this script.' % path)
            exit(1)


def setupSmallRunDirectory(parentDirectory=""):
    """
    Generate directories and symlink required for the assembler to run
    """
    # Generate absolute paths to the files that will be created
    dataPath = os.path.abspath(os.path.join(parentDirectory, "Data"))
    dataOnDiskPath = os.path.abspath(os.path.join(parentDirectory, "DataOnDisk"))
    threadLogsPath = os.path.abspath(os.path.join(parentDirectory, "threadLogs"))

    # Create the directories.
    os.mkdir(dataPath)
    os.mkdir(threadLogsPath)
    
    # Create the symbolic links.
    os.symlink(dataPath, dataOnDiskPath)


def main():
    verifyArguments()
    verifyDirectoryFiles()
    setupSmallRunDirectory()


if __name__ == "__main__":
    main()