#!/usr/bin/python3

"""

This run the final portion of Mode 2 assembly.

"""

import ast
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
    phasingMinReadCount = int(config['Assembly']['phasing.minReadCount']),
    suppressGfaOutput = ast.literal_eval(config['Assembly']['suppressGfaOutput']),
    suppressFastaOutput = ast.literal_eval(config['Assembly']['suppressFastaOutput']), 
    suppressDetailedOutput = ast.literal_eval(config['Assembly']['suppressDetailedOutput']), 
    suppressPhasedOutput = ast.literal_eval(config['Assembly']['suppressPhasedOutput']), 
    suppressHaploidOutput = ast.literal_eval(config['Assembly']['suppressHaploidOutput']) 
    )




