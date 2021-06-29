#!/usr/bin/python3

import shasta
import argparse


parser = argparse.ArgumentParser(description=
    'Write a FASTA file containing all the reads present in a local read graph.')
    
parser.add_argument('--readId', type=int, required=True)
parser.add_argument('--strand', type=int, choices=range(2), required=True)
parser.add_argument('--maxDistance', type=int, required=True)
parser.add_argument('--useReadName', action='store_true')
parser.add_argument('--allowChimericReads', action='store_true')
parser.add_argument('--allowCrossStrandEdges', action='store_true')
parser.add_argument('--allowInconsistentAlignmentEdges', action='store_true')
arguments = parser.parse_args()

readId = arguments.readId
strand = arguments.strand
maxDistance = arguments.maxDistance
useReadName = arguments.useReadName
allowChimericReads = arguments.allowChimericReads
allowCrossStrandEdges = arguments.allowCrossStrandEdges
allowInconsistentAlignmentEdges = arguments.allowInconsistentAlignmentEdges

a = shasta.Assembler()
a.accessAlignmentData()
a.accessReadGraph()
a.accessMarkers()
a.writeLocalReadGraphReads(
    readId=readId, strand=strand, 
    maxDistance=maxDistance,
    useReadName=useReadName,
    allowChimericReads=allowChimericReads,
    allowCrossStrandEdges=allowCrossStrandEdges,
    allowInconsistentAlignmentEdges=allowInconsistentAlignmentEdges)


