#!/usr/bin/python3

import shasta

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--segmentId', type=int, required=True)
arguments = parser.parse_args()

a = shasta.Assembler()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphOrientedReadsByEdge()
a.colorGfaBySimilarityToSegment(arguments.segmentId)





