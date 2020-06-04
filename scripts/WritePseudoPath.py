#!/usr/bin/python3

import shasta

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--readId', type=int, required=True)
parser.add_argument('--strand', type=int, choices=range(2), required=True)
arguments = parser.parse_args()

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.writePseudoPath(readId=arguments.readId, strand=arguments.strand)





