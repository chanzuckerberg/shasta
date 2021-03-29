#!/usr/bin/python3

import shasta
import sys

helpMessage = """
Write all alignment candidates (pairs of ReadId) to a file named AlignmentCandidates.csv. Takes no arguments.
"""

if not len(sys.argv) == 1:
    print(helpMessage)
    exit(1)


a = shasta.Assembler()
print("Accessing data ...")
a.accessMarkers()
a.accessAlignmentCandidates()
print("Writing to AlignmentCandidates.csv ...")
a.writeAlignmentCandidates()
print("Done")
