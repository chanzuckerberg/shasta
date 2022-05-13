#!/usr/bin/python3

import shasta
import argparse

parser = argparse.ArgumentParser(description=
    'Create an assembly path in the mode 3 assembly graph.')
parser.add_argument('segmentId', type=int, help='The segment id to start from.')
parser.add_argument('direction', type=int, help='The path direction (0=forward, 1=backward).')
arguments = parser.parse_args()


a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphEdges()
a.accessMode3AssemblyGraph()

path = a.createMode3AssemblyPath(arguments.segmentId, arguments.direction)

print('Found the following path: ', *path)
