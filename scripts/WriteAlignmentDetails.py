#!/usr/bin/python3

import shasta
import argparse


parser = argparse.ArgumentParser(description=
    'Write CSVs with details for each alignment')
    
arguments = parser.parse_args()

a = shasta.Assembler()
a.accessMarkers()
a.accessAlignmentCandidates()
a.accessCompressedAlignments()
a.accessAlignmentData()
a.writeAlignmentDetails()


