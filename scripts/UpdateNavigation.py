#!/usr/bin/python3

import glob
import sys

helpMessage = """
This is used to update the navigation section in all html files.
It uses the contents of Navigation.html
to update the corresponding section of each html file.

This should be rerun every time Navigation.html changes.

Invoke from the docs directory, without arguments.
"""

# Check that there are no arguments.
if not len(sys.argv) == 1:
    print(helpMessage)
    exit(1)

# Define the strings marking the portion to be replaced.    
beginString = '    <nav role="navigation">'
endString =   '    </nav>\n'

# Get the replacement string.
file = open('Navigation.html', 'r')
replacementString = file.read()
file.close()
    
# Loop over all html files in the current directory.    
fileNames = glob.glob('*.html')
for fileName in fileNames:

    # Skip Navigation.html.
    if fileName == 'Navigation.html':
        continue

    # Get the contents of this file.
    file = open(fileName, 'r')
    fileContents = file.read()
    file.close()
    
    # Locate the begin and end string.
    beginPosition = fileContents.find(beginString)
    endPosition   = fileContents.find(endString) + len(endString)
    
    # Skip if the begin and end string were not found
    # or are in inconsistent positions.
    if beginPosition==-1 or endPosition==-1 or beginPosition>endPosition:
        print('Skipping', fileName)
        continue
    
    # Replace the portion between the begin string and the end string.
    fileContents = fileContents[:beginPosition] + replacementString + fileContents[endPosition:]

    # Write out the modified file contents.
    file = open(fileName, 'w')
    file.write(fileContents)
    file.close()

