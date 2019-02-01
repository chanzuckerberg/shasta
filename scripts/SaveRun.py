#!/usr/bin/python3

import os
import shutil
import sys


def parseArguments():
    # Get from the arguments the list of input fasta files and check that they all exist.
    helpMessage = "This script copies a run to directory DataOnDisk."
    if not len(sys.argv)==1:
        print(helpMessage)
        exit(1)


def saveRun(parentDirectory=""):
    dataPath = os.path.abspath(os.path.join(parentDirectory, "Data"))
    dataOnDiskPath = os.path.abspath(os.path.join(parentDirectory, "DataOnDisk"))

    if not os.path.lexists(dataPath):
        raise Exception('Data does not exist.')
    
    if os.path.lexists(dataOnDiskPath):
        raise Exception('DataOnDisk already exists. Remove before running this script.')
        
    shutil.copytree('Data', 'DataOnDisk')


def main():
    parseArguments()
    saveRun()


if __name__ == "__main__":
    main()
