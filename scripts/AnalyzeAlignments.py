#!/usr/bin/python3

import shasta
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('readId', type=int)
parser.add_argument('strand', type=int, choices=range(2))
arguments = parser.parse_args()

a = shasta.Assembler()
a.accessMarkers()
a.accessAlignmentData()
a.accessCompressedAlignments()
a.analyzeAlignments(readId = arguments.readId, strand = arguments.strand)





