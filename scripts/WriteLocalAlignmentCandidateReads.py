#!/usr/bin/python3

import shasta
import argparse


parser = argparse.ArgumentParser(description=
    'Write a FASTA file containing all the reads present in a local alignment candidate graph.')
    
parser.add_argument('--readId', type=int, required=True)
parser.add_argument('--strand', type=int, choices=range(2), required=True)
parser.add_argument('--maxDistance', type=int, required=True)
parser.add_argument('--allowChimericReads', action='store_true')
parser.add_argument('--allowCrossStrandEdges', action='store_true')
parser.add_argument('--allowInconsistentAlignmentEdges', action='store_true')
arguments = parser.parse_args()

readId = arguments.readId
strand = arguments.strand
maxDistance = arguments.maxDistance
allowChimericReads = arguments.allowChimericReads
allowCrossStrandEdges = arguments.allowCrossStrandEdges
allowInconsistentAlignmentEdges = arguments.allowInconsistentAlignmentEdges

a = shasta.Assembler()
a.accessMarkers()
a.accessAlignmentCandidates()
a.accessAlignmentCandidateTable()
a.accessAlignmentData()
a.accessReadGraph()
a.writeLocalAlignmentCandidateReads(
    readId=readId, strand=strand,
    maxDistance=maxDistance,
    allowChimericReads=allowChimericReads,
    allowCrossStrandEdges=allowCrossStrandEdges,
    allowInconsistentAlignmentEdges=allowInconsistentAlignmentEdges)


