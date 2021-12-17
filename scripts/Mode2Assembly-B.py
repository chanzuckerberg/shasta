#!/usr/bin/python3

"""

This run the final portion of Mode 2 assembly.

"""

import ast
import shasta
import GetConfig

config = GetConfig.getConfig()

shasta.openPerformanceLog('Mode2Assembly-B.log')

a = shasta.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphEdges()
a.accessMarkerGraphReverseComplementEdge()
a.accessMarkerGraphConsensus()

a.createAssemblyGraph2(
    superbubbleRemovalEdgeLengthThreshold = int(config['Assembly']['superbubbleRemoval.edgeLengthThreshold']),
    pruneLength = int(config['Assembly']['pruneLength']),
    phasingMinReadCount = int(config['Assembly']['phasing.minReadCount']),
    suppressGfaOutput = ast.literal_eval(config['Assembly']['mode2.suppressGfaOutput']),
    suppressFastaOutput = ast.literal_eval(config['Assembly']['mode2.suppressFastaOutput']), 
    suppressDetailedOutput = ast.literal_eval(config['Assembly']['mode2.suppressDetailedOutput']), 
    suppressPhasedOutput = ast.literal_eval(config['Assembly']['mode2.suppressPhasedOutput']), 
    suppressHaploidOutput = ast.literal_eval(config['Assembly']['mode2.suppressHaploidOutput']) 
    )




