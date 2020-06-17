#!/usr/bin/python3

import shasta
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('segmentId', type=int)
parser.add_argument('direction', choices=['forward', 'backward'])
arguments = parser.parse_args()



a = shasta.Assembler()
a.accessMarkers()
a.accessReadGraph()
a.accessReadFlags()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphEdgeLists()
a.accessAssemblyGraphSequences()
a.followOrientedReadPaths(arguments.segmentId, arguments.direction == 'forward')

	





