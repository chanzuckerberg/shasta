#!/usr/bin/python3

import shasta
import argparse


parser = argparse.ArgumentParser(description=
    'Set the flags of all marker graph edges to specified values, or specify 2 to leave a flag unchanged.')
    
parser.add_argument('--wasRemovedByTransitiveReduction', type=int, default=2, choices=range(3))
parser.add_argument('--wasPruned', type=int, default=2, choices=range(3))
parser.add_argument('--isSuperBubbleEdge', type=int, default=2, choices=range(3))
parser.add_argument('--isLowCoverageCrossEdge', type=int, default=2, choices=range(3))
parser.add_argument('--wasAssembled', type=int, default=2, choices=range(3))

arguments = parser.parse_args()



a = shasta.Assembler()
a.accessMarkerGraphEdges(accessEdgesReadWrite = True)
a.setMarkerGraphEdgeFlags(
    wasRemovedByTransitiveReduction = arguments.wasRemovedByTransitiveReduction,
    wasPruned = arguments.wasPruned,
    isSuperBubbleEdge = arguments.isSuperBubbleEdge,
    isLowCoverageCrossEdge = arguments.isLowCoverageCrossEdge,
    wasAssembled = arguments.wasAssembled)


