#!/usr/bin/python3

import os
import shasta

# Find the path to the docs directory.
thisScriptPath = os.path.realpath(__file__)
thisScriptDirectory = os.path.dirname(thisScriptPath)
thisScriptParentDirectory = os.path.dirname(thisScriptDirectory)
docsDirectory = thisScriptParentDirectory + "/docs"

# Initialize the assembler and access what we need.
a = shasta.Assembler(
    smallDataFileNamePrefix='dataOnDisk/',
    largeDataFileNamePrefix='DataOnDisk/')
a.accessReadsReadOnly()
a.accessReadNamesReadOnly()
a.accessKmers()
a.accessMarkers()
a.accessOverlaps()
a.accessAlignmentData()
a.accessMarkerGraphVertices()
a.accessMarkerGraphConnectivity()

a.setDocsDirectory(docsDirectory)
a.explore()


