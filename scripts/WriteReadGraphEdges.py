#!/usr/bin/python3

import shasta
import sys

helpMessage = """
Write all read graph edges (pairs of ReadId) to a file named ReadGraphEdges.csv. Takes no arguments.
"""

if not len(sys.argv) == 1:
    print(helpMessage)
    exit(1)


a = shasta.Assembler()
print("Accessing data ...")
a.accessReadGraph()
print("Writing to ReadGraphEdges.csv ...")
a.writeReadGraphEdges()
print("Done")
