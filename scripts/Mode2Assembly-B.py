#!/usr/bin/python3

"""

This run the final portion of Mode 2 assembly.

"""

import shasta
import GetConfig

config = GetConfig.getConfig()

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphEdges()
a.accessMarkerGraphReverseComplementEdge()
a.accessMarkerGraphConsensus()

a.createAssemblyGraph2(
    bubbleRemovalDiscordantRatioThreshold = float(config['Assembly']['bubbleRemoval.discordantRatioThreshold']),
    bubbleRemovalAmbiguityThreshold = float(config['Assembly']['bubbleRemoval.ambiguityThreshold']),
    bubbleRemovalMaxPeriod = int(config['Assembly']['bubbleRemoval.maxPeriod']),
    superbubbleRemovalEdgeLengthThreshold = int(config['Assembly']['superbubbleRemoval.edgeLengthThreshold']),
    phasingMinReadCount = int(config['Assembly']['phasing.minReadCount']))




