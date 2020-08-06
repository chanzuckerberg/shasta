#!/usr/bin/python3

import shasta

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('readId0', type=int)
parser.add_argument('strand0', type=int, choices=range(2))
parser.add_argument('readId1', type=int)
parser.add_argument('strand1', type=int, choices=range(2))
arguments = parser.parse_args()


a = shasta.Assembler()
a.accessMarkers()
a.accessReadGraph()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphEdgeLists()
a.alignPseudoPaths(arguments.readId0, arguments.strand0, arguments.readId1, arguments.strand1)

	





