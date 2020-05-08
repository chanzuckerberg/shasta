#!/usr/bin/python3

import shasta

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--readId0', type=int, required=True)
parser.add_argument('--strand0', type=int, choices=range(2), required=True)
parser.add_argument('--readId1', type=int, required=True)
parser.add_argument('--strand1', type=int, choices=range(2), required=True)
arguments = parser.parse_args()

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.colorGfaWithTwoReads(
	readId0=arguments.readId0, strand0=arguments.strand0, 
	readId1=arguments.readId1, strand1=arguments.strand1)





