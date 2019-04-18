#!/usr/bin/python3

import os
import sys


def parseArguments():
    helpMessage = """
    This unmounts the Data directory for the Shasta run in the current directory.
    Invoke without arguments.
    It uses sudo to invoke a command that requires root privileges,
    so depending on your sudo set up it may ask for your password.
    """

    if not len(sys.argv) == 1:
        print(helpMessage)
        exit(1)


def cleanUpRunDirectory(requireUserInput=True):
    if requireUserInput:
        response = input('This will destroy binary data in memory in the Data directory.\nEnter (Yes/yes/Y/y) if you really want to do that:\n')   
        if not response in set(['Yes', 'Y', 'y', 'yes']):
            print('Nothing done.')
            exit(0)
            
    os.system('sudo umount Data')  
    os.rmdir('Data')          


def main():
    parseArguments()
    cleanUpRunDirectory()
    
if __name__ == "__main__":
    main()
